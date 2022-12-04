import numpy as np
from datetime import datetime
from os import path, makedirs
from sklearn.gaussian_process import GaussianProcessRegressor # is it needed?
from sklearn.gaussian_process.kernels import RBF # is it needed?
import scipy.special
from scipy.special import eval_legendre, eval_hermitenorm
from numpy.linalg import inv   #if needed

class PC_Kriging():
    def __init__(self, config=None):
        if config is None:
            config = {"pol_type": ['hermite', 'legendre'],
                      "order": [2, 1]}
        assert "pol_type" in config and \
               "order" in config and \
            "Missing config"

        self.date_record = datetime.now().strftime("%Y_%m_%d_%H%M%S")
        self.M1 = self.hermite
        self.M2 = self.hermite

    def train (self, xn, y, p, theta):
        #xn, inputs training points
        #y, observations
        #p, degree of polynomial expansion
        #theta [ l, hyperparameter (length scaled); v, matern coefficient ]
        #---------------------------------------------------
        
        # information matrix based on observations and a chosen polynomial expansion
        
        dim = xn.shape[1]
        indices = np.arange(0, p+1)     #polynomials indices to be combined (1 colum ver variable)
        comb = np.zeros((p**dim, dim))  #to store all possible combinations
        
        # CHECK DIMENSIONS (INCREASE 1 indices vector per dimension)
        #------------------------------------------------------------------------- #all possible combinations of indices
        comb = np.array(np.meshgrid(indices,indices)).T.reshape(-1,dim)     #(indices, indices, indices...) per dimension
        alpha = []                 #list of product combinations with degree <= "p" (truncation term)

        for i in range (0,len(comb)): 
            if (np.sum(comb[i])) <= p:  
                alpha.append(comb[i])

        alpha = np.array(alpha) #alpha for multivariate combination
        
        phi = self.info_matrix(xn, alpha)     #check the number of variables included 
        
        B, sig2 = self.coefficients(phi, xn, y, theta[0], theta[1])
        
        return B, sig2, phi, alpha

    def predict (self, XN, xn, y, theta, modelpar1):
        # B, sig2, phi, alpha 
        B    = modelpar1[0]
        sig2 = modelpar1[1]
        phi  = modelpar1[2]
        alpha = modelpar1[3]
        
        #allocanting storage ------------------------------------------
        fx = np.zeros((len(XN),len(alpha[1])))
        rx = np.zeros((len(XN),len(xn)))
        Rn = np.zeros((len(xn),len(xn)))
        
        mean1 = np.zeros(len(XN))
        mean2 = np.zeros(len(XN))
        PCKmean = np.zeros(len(XN))
            
        ux = np.zeros((len(alpha[1]),len(XN)))
        term1 = np.zeros((len(XN),len(XN)))
        term2 = np.zeros((len(XN),len(XN)))
        variance = np.zeros((len(XN),len(XN)))
        
        #XN, inputs  predictions ---------------------------------------
        fx = self.info_matrix(XN, alpha)                # f(x) information matrix about the predictions
        rx = self.matern(XN, xn, theta[0], theta[1])       # r(x) correlation matrix between predictions and observations
        Rn = self.matern(xn, xn, theta[0], theta[1])       # R    correlation matrix between predictions
        Rn_inv = np.linalg.inv(Rn)                     # R inverse

        #---------------------------------------------- Mean prediction
        mean1 = fx @ B 
        mean2 = rx @ Rn_inv @ (y - (phi @ B))
        PCKmean = mean1 + mean2
        #---------------------------------------------- Variance prediction
        ux = ( phi.T @ Rn_inv @ rx.T) - fx.T

        term1 = rx @ Rn_inv @ rx.T
        term2 = ux.T @ np.linalg.inv(phi.T @ Rn_inv @ phi) @ ux

        variance = sig2 * ( 1 - term1 + term2)

        variaDiag = np.diagonal(variance) 
        
        return PCKmean, variaDiag


    def info_matrix (self, X, alpha):    # X must be normalized for polynomial evaluations
    
        #alpha = polynomial indices to be combined
        
        fx = np.zeros ((len(X), len(alpha))) 

        for i in range (0,len(X)):
            for j in range (0, len(alpha)):
                fx[i,j] = self.M1(X[i,0], alpha[j,0]) * self.M2(X[i,1], alpha[j,1])  #it has to increase per dimension 

        return fx

    def coefficients(self, F, xn, y, l, v):   
    
        R = self.matern(xn , xn, l, v)    # Correlation matrix R between observations
        
        R_inv = np.linalg.inv(R)
        left_r =  np.linalg.inv((F.T @ R_inv) @ F)
        right_r = F.T @ R_inv @ y
        B_hat = left_r @ right_r
        #----------------------------------------------------------
        ins = y - (F @ B_hat)              # Scale sigma^2
        left_s = (1/len(xn))*((ins).T)
        right_s = R_inv @ ins
        sig2 = left_s @ right_s
        
        return  B_hat, sig2

    def matern(self, xr, xn, theta, v):     
        #theta, l, hyperparameter (length scaled)
        #matern parameter, v can be 3/2, 5/2
        
        d = self.distance(xr,xn)                 #eucledian distance
        R = np.zeros((len(xr),len(xn)))
        
        if v == 3/2:
            R = (1+ np.sqrt(3)*d/theta) * np.exp(-(np.sqrt(3)*d/theta))
        elif v == 5/2:
            R = ((1+ (np.sqrt(5)*d/theta) + (5/3)*(d/theta)**2 )* np.exp(-(np.sqrt(5)*d/theta)))
        return R

    def distance(self, x, xk):                        #multidimensional distance between 2 samples
        d = np.zeros((len(x),len(xk)))
        
        for j in range(0,len(xk)):
            for i in range(0,len(x)):
                d[i,j] = d[i,j] + np.sqrt(np.sum((x[i]-xk[j])**2))
        return d

    def hermite(self, x, n):
        return eval_hermitenorm(n, x) / np.sqrt(scipy.special.factorial(n))

    def legendre(self, x, n):
        return eval_legendre(n,x)*np.sqrt(2*n+1)

    def scalelegendre(self, x, new_min, new_max): 
        return ((new_min+new_max)+((new_max-new_min)*x))/2

    def scalehermite(self, x , mean, sigma):
        return mean+sigma*x