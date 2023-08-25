import numpy as np
from datetime import datetime
from os import path, makedirs
import scipy.special
from scipy.special import eval_legendre, eval_hermitenorm
from numpy.linalg import inv
from scipy import optimize
import time
import math

class PC_Kriging():
    def __init__(self, config=None):
        if config is None:
            config = {"pol_type": ['hermite', 'legendre']}
        assert "pol_type" in config and \
            "Missing config"

        self.date_record = datetime.now().strftime("%Y_%m_%d_%H%M%S")
        
        self.polynomials = []
        for polyn in config["pol_type"]:
            if polyn == 'hermite':
                self.polynomials.append(self.hermite)
            else:
                self.polynomials.append(self.legendre)

    def train (self, xn, y, p, l, v):
        self.doe = xn          #xn, inputs training points
        self.N_samples = xn.shape[0]   #number of samples
        self.observ = y        #y, observations
        self.degree = p        #p, degree of polynomial expansion
        self.maternCoef = v    #v, matern coefficient 
        self.theta = l         #theta, hyperparameter (length scale, l)
        
        # information matrix based on observations and a chosen polynomial expansion
        
        dim = xn.shape[1]
        
        indices = np.arange(0, p+1)     #polynomials indices to be combined (1 colum ver variable)
        comb = np.zeros((p**dim, dim))  #to store all possible combinations
        
        # CHECK DIMENSIONS (INCREASE 1 indices vector per dimension)
        #------------------------------------------------------------------------- #all possible combinations of indices
        indices_ = []
        for i in range(dim):
            indices_.append(indices)
        comb = np.array(np.meshgrid(*indices_)).T.reshape(-1,dim)
        
        alpha = []                 #list of product combinations with degree <= "p" (truncation term)

        for i in range (0,len(comb)): 
            if (np.sum(comb[i])) <= p:  
                alpha.append(comb[i])

        alpha = np.array(alpha) #alpha for multivariate combination
        
        phi = self.info_matrix(xn, alpha)     #check the number of variables included
        
        #coefficients, calibrated weights
        #Sigma Squared
        self.coeff, self.sigmaSQ = self.coefficients(phi, xn, y, l, v)
        #--------------------------------------------------- 
        self.Poly_ind = alpha      #polynomial indices
        self.InfoMat = phi         #information matrix with training points
        
        return self.coeff, self.sigmaSQ

    def predict(self, XN):
        Rn = self.matern(self.doe, self.doe, self.theta, self.maternCoef)
        Rn_inv = np.linalg.inv(Rn)
        
        size_XN = XN.shape
        mean_predict = np.zeros(size_XN[0])
        var_predict = np.zeros(size_XN[0])
    
        for i in range(size_XN[0]):
            fx = self.info_matrix(XN[i].reshape(1,size_XN[1]), self.Poly_ind)
            rx = self.matern_fast(XN[i].reshape(1,size_XN[1]), self.doe, self.theta, self.maternCoef)

            mean1 = fx @ self.coeff
            mean2 = rx @ Rn_inv @ (self.observ - (self.InfoMat @ self.coeff))
            mean_total = mean1 + mean2
                
            ux = ( self.InfoMat.T @ Rn_inv @ rx.T) - fx.T
            term1 = rx @ Rn_inv @ rx.T
            term2 = ux.T @ np.linalg.inv(self.InfoMat.T @ Rn_inv @ self.InfoMat) @ ux
            variance = self.sigmaSQ * ( 1 - term1 + term2)
        
            mean_predict[i] = mean_total
            var_predict[i] = variance[0]
        
        return mean_predict , var_predict
        
    def matern_fast(self, xr, xn, l, v):
        #l, hyperparameter (length scaled)
            #matern parameter, v can be 3/2, 5/2
        size_xn = xn.shape
        d = np.linalg.norm(xr - xn, axis=1).reshape(1, size_xn[0])
      
        if v == 3/2:
            R = (1+ np.sqrt(3)*d/l) * np.exp(-(np.sqrt(3)*d/l))
        elif v == 5/2:
            R = ((1+ (np.sqrt(5)*d/l) + (5/3)*(d/l)**2 )* np.exp(-(np.sqrt(5)*d/l)))
            
        return R

    def info_matrix (self, X, alpha):    # X must be normalized for polynomial evaluations
    
        #alpha = polynomial indices to be combined
        
        fx = np.zeros ((len(X), len(alpha))) 

        for i in range (0,len(X)):
            for j in range (0, len(alpha)):
                indx = 0
                fx_ = 1
                for polyn in self.polynomials:
                    fx_ *= polyn(X[i,indx], alpha[j,indx])
                    indx += 1
                fx[i,j] = fx_  #it has to increase per dimension
        return fx

    def coefficients(self, F, xn, y, l, v):   
        vnug= 0 #1e-9                  #adding nuget to main diagonal in case of numerical instability
        nug = np.eye(len(xn))*vnug

        R = self.matern(xn, xn, l, v)+ nug     # Correlation matrix R between observations
        
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

    def matern(self, xr, xn, l, v):     
        #l, hyperparameter (length scaled)
        #matern parameter, v can be 3/2, 5/2
        
        d = self.distance(xr, xn)                 #eucledian distance
        
        if v == 3/2:
            R = (1+ np.sqrt(3)*d/l) * np.exp(-(np.sqrt(3)*d/l))
        elif v == 5/2:
            R = ((1+ (np.sqrt(5)*d/l) + (5/3)*(d/l)**2 )* np.exp(-(np.sqrt(5)*d/l)))
        return R
    
    def distance(self, x, xk):                        #multidimensional distance between 2 samples with numpy      
        d = np.zeros((len(x),len(xk)))
        
        for i in range(0,len(x)):
            
            d[i] = [np.linalg.norm(x[i] - xk[j]) for j in range(0, len(xk))]
            
        return d

    def hermite(self, x, n):
        return eval_hermitenorm(n, x) / np.sqrt(scipy.special.factorial(n))

    def legendre(self, x, n):
        return eval_legendre(n,x)*np.sqrt(2*n+1)

    def scalelegendre(self, x, new_min, new_max):
        #to go from (-1,1) to (new_min, new_max) 
        return ((new_min+new_max)+((new_max-new_min)*x))/2

    def scalehermite(self, x , mean, sigma):
        return mean+sigma*x
    
    def LinearNorm(self, x, oldmin, oldmax, newmin, newmax):    # scaling linearly X to new domain limits
        return newmin + ((x-oldmin)*(newmax-newmin)/(oldmax-oldmin))
    
    def objective_F(self, l):

        corr = self.matern(self.doe , self.doe, l, self.maternCoef)
        detR = np.linalg.det(corr)
        modelpar = self.train(self.doe, self.observ, self.degree, l, self.maternCoef)    # returns B, sig2
        ### -------Theta_ by UQLab User Manual PCK(C. Lataniotis, D. Wicaksono, S. Marelli, B. Sudret)
        sig2 = modelpar[1].reshape(-1)
        # return 0.5*(np.log(detR)+ N*np.log(2*np.pi*sig2)+ N)

        ### ------------------Theta_ by MLE PCK(Schobi,Sudret,Wiart)------------------------------
        FB = self.InfoMat @ modelpar[0]
        ins = (self.observ - FB).reshape(-1)
        R_1 = np.linalg.inv(corr)
        return ((ins.T) @ R_1 @ ins) * ( 1 /self.N_samples) * (detR**(1/self.N_samples))
    
    def optimize(self, method):
        # available methods
        # 'shgo' stands for “simplicial homology global optimization”
        # 'DA' stands for "dual_annealing"
        # 'DE' differential evolution

        if method is None:
            method = 'shgo'

        assert method and \
            "Optimizer not specified, 'shgo', 'DA', 'DE' "
        
        d = self.distance(self.doe, self.doe)
        
        lmax = np.max(d)
        lmin = np.min(d[d!=0])
        bounds = [(lmin, lmax)]

        if method == 'shgo':
            return optimize.shgo(self.objective_F, bounds)['x'][0]

        elif method ==  'DA':
            return optimize.dual_annealing(self.objective_F, bounds)['x'][0]
        
        elif method ==  'DE':
            return optimize.differential_evolution(self.objective_F, bounds)['x'][0]