# Python implementation that can be compiled into C using numba to make it faster
# based on the Matlab code of jgsa@civil.aau.dk
# As much as it was possible the code and variable naming are kept to allow for the verification of the implementation.#
# Matlab comments:
# black_box(X)
# this function returns a vector with the value of the
# performance function D=1-L/R (L=Load,R=resistance)
# the limit state function g=D
# the failure region is defined as g<0
# the input X must contain 14 columns which represents each random variable
# each row represents a simulation
# email: jgsa@civil.aau.dk
# whatsapp: +56987226332

import numpy as np
import pickle
from numba import jit


# @jit is not speeding this up (tested with moving the loading outside)
def juan_first_passage(X):
    # system properties
    # m_values: the mean value
    # cv_values: coefficient of variation
    m_values = np.array([np.pi, 0.03, 5, 0.05])
    cv_values = np.array([0.1, 0.3, 0.2, 0.2])
    gg, x2p = load_excitation()
    # threshold level
    R = 0.45
    # predefine limit state function
    n_row = X.shape[0]
    g = np.empty(n_row)
    for i in range(n_row):
        # Z corresponds to the random variables of the excitation model
        Z = X[i, 0:10]
        # Q corresponds to the random variables of the system model
        Q = X[i, 10:14]
        # generate the excitation
        dt, t, acc = excitation_model(Z, gg, x2p)
        # calculate the system response
        x_aux = fem_response(Q, acc, t, dt, m_values, cv_values)
        # calculate the limit state function
        g[i] = 1 - x_aux/R
    return g


def load_excitation():
    # load the matrix file which contains the information
    # of the orthogonal expansion
    # with open(path + 'x2p-g.pkl', 'rb') as f:
    x = pickle.load(open('code_rp/gfun/juan_first_passage/x2p-g.pkl', 'rb'))
    gg = x['g']
    x2p = x['x2p']
    return gg, x2p


# excitation model(Z) generate and returns the excitation
# according to the excitation model
# for testing
#   Z = np.ones(10)
#   excitation_model(Z)
@jit(nopython=True, cache=True)
def excitation_model(Z, gg, x2p):
    # time steps
    dt = 0.01
    # define time history
    n = np.ceil(20/dt)
    t = np.arange(n+1)*20./n
    # "gg" and "x2p" are matrices loaded from file 'x2p-g'
    # define the excitation "acc"
    F = np.sum((np.real(x2p).T * Z).T, 0) * gg
    F = F.flatten()
    acc = F * np.sqrt(81.15/10000)
    return dt, t, acc


# fem_response() returns the maximum response of the FEM model
@jit(nopython=True, cache=True)
def fem_response(U, acc, t, dt, m_values, cv_values):
    # general numerical data
    b_ = 0.25
    g_ = 0.50
    tol = 1e-7
    nGL = 1.
    N = len(t)
    gen_data = [dt, g_, b_, tol, nGL]

    # transformation from U-space to real space
    T = np.empty(4)
    T[:] = np.nan
    T[0] = u2log(m_values[0], cv_values[0], U[0])
    T[1] = u2log(m_values[1], cv_values[1], U[1])
    T[2] = u2log(m_values[2], cv_values[2], U[2])
    T[3] = u2log(m_values[3], cv_values[3], U[3])

    # definition of system properties
    ke = T[2]
    db = 0.
    alpha = 0.00
    Uy = T[3]
    alpha_i = 1.
    beta_i = 1.
    gamma_i = 1.
    n = 1.
    ais_data = [ke, db, alpha, Uy, alpha_i, beta_i, gamma_i, n]
    Ms = 1
    Ks = T[0]**2 + ke*alpha
    Cs = 2*np.sqrt(Ks)*T[1]
    b_geo = 1

    # matrices or scalar to use Newmark method
    M_0 = Ms/dt**2 + g_*Cs/dt + b_*Ks
    M_1 = Ms/dt**2 + g_*Cs/dt
    M_2 = Ms/dt + (g_ - b_)*Cs
    M_3 = (0.5 - b_)*Ms + dt/2*(g_ - 2*b_)*Cs
    M_11 = (nGL*b_) / M_0
    M_22 = M_1 / M_0
    M_33 = M_2 / M_0
    M_44 = M_3 / M_0
    M_1a = b_geo

    # generated excitation
    f = acc*Ms

    # dissipation parameters / NOT USED mu=0 / must be defined
    W = 0.1
    mu = 0.
    dis_data = [W, mu]

    # Solve equation of movement | MATLAB version | uncomment next line to use
    x = mdof_rk(M_1a, M_11, M_22, M_33, M_44, f, t, ais_data, dis_data, gen_data, N)

    # returns the maximum response
    x_max = np.abs(x).max()

    return x_max


# transform variables from U space to T space
# testing:
#   m, cv, U = np.pi, 0.1, 1
#   u2log(m, cv, U)
@jit(nopython=True, cache=True)
def u2log(m, cv, U):
    v = (cv*m)**2
    mu = np.log(m**2/np.sqrt(v + m**2))
    s = np.sqrt(np.log(v/m**2 + 1))
    X = np.exp(mu + s*U)
    return X


# mdof_rk solve the equation of movement of the dynamical FEM
@jit(nopython=True, cache=True)
def mdof_rk(M_1a, M_11, M_22, M_33, M_44, f, t, ais_data, dis_data, gen_data, N):
    # collect general and dissipation data
    tol = gen_data[3]
    nGL = gen_data[4]
    # dissipation data
    W = dis_data[0]
    mu = dis_data[1]
    # predefine vectors
    nt = len(t)
    x = np.zeros(nt)
    xp = np.zeros(nt)
    x2p = np.zeros(nt)
    fis1 = np.zeros(nt)
    fis2 = np.zeros(nt)
    z = np.zeros(nt)
    F = np.zeros(nt)
    # load initial load
    F[0] = f[0]
    # auxiliary scalar
    k_aux_1 = (1 - ais_data[2]) * ais_data[0] * ais_data[3]
    for i in range(1, N):
        # update auxiliary variable
        z[i] = z[i-1]
        # counter to 1
        cont = 1
        # defines a large initial error
        error = 1
        # save the previous velocity
        xp_a = xp[i-1]

        while error > tol and cont < 30:
            # define the np.sign of xp
            signxp = 1 - 2 * (xp_a < 0)
            # define external force
            fis1[i] = k_aux_1 * z[i]
            fis2[i] = signxp * W * mu
            # define actual time step state
            F[i] = M_1a * f[i] - (fis1[i] + fis2[i])
            x[i], xp[i], x2p[i] = newmark(F[i], M_11, M_22, M_33, M_44, x[i-1], xp[i-1], x2p[i-1], gen_data)
            z[i] = runge_kutta(ais_data, [xp[i-1], xp[i]], [z[i-1], z[i]], gen_data)
            # verify error
            xp_nl_p = xp[i]
            error = abs((xp_nl_p - xp_a) / xp_nl_p)
            xp_a = xp_nl_p
            # update counter
            cont = cont + 1
    return x


# runge_kutta() solve differential equation of auxiliary variable
# using Runge Kutta 4
@jit(nopython=True, cache=True)
def runge_kutta(ais_data, xp, z, gen_data):
    # Collect general data
    dt = gen_data[0]
    Uy = ais_data[3]
    alpha_i = ais_data[4]
    beta_i = ais_data[5]
    gamma_i = ais_data[6]
    n = ais_data[7]
    # define intermediate time step
    xp_m = (xp[0] + xp[1])/2
    # defines next time step variables
    K1 = xp[0]*(alpha_i-(z[1])**n *        (beta_i*np.sign(z[1])         + gamma_i*np.sign(xp[0])))/Uy
    K2 = xp_m*(alpha_i-(z[1]+K1*dt/2)**n * (beta_i*np.sign(z[1]+K1*dt/2) + gamma_i*np.sign(xp_m))) /Uy
    K3 = xp_m*(alpha_i-(z[1]+K2*dt/2)**n * (beta_i*np.sign(z[1]+K2*dt/2) + gamma_i*np.sign(xp_m))) /Uy
    K4 = xp[1]*(alpha_i-(z[1]+K3*dt)**n *  (beta_i*np.sign(z[1]+K3*dt)   + gamma_i*np.sign(xp[1])))/Uy
    # return next time step
    z_ = z[0] + dt/6*(K1 + 2*K2 + 2*K3 + K4)
    return z_


# Newmark method to solve equation of movement
@jit(nopython=True, cache=True)
def newmark(F, M_11, M_22, M_33, M_44, x, xp, x2p, gen_data):
    # Collect general data
    x = list([x])
    xp = list([xp])
    x2p = list([x2p])
    dt = gen_data[0]
    g_ = gen_data[1]
    b_ = gen_data[2]
    # define next time step
    x_i = M_11*F + M_22*x[0] + M_33*xp[0] + M_44*x2p[0]
    x.append(x_i)
    x2p_i = 1/(b_*dt**2)*(x[1] - x[0]) - 1/(b_*dt)*xp[0] - (1/(2*b_) - 1)*x2p[0]
    x2p.append(x2p_i)
    xp_i = xp[0] + ((1 - g_)*x2p[0] + g_*x2p[1])*dt
    xp.append(xp_i)
    # return next time step
    x_ = x[1]
    xp_ = xp[1]
    x2p_ = x2p[1]
    return x_, xp_, x2p_

# X = np.ones(14)*-1
# g = juan_first_passage(X)
# print(g)