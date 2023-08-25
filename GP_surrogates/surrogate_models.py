import numpy as np
from pc_kriging import PC_Kriging
from datetime import datetime
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from scipy.stats import norm
from doepy import build
from scipy import optimize
import pickle

class surrogate_models():

    def __init__(self, settings=None):
        
        self.metamodeltype = settings["metamodel"] if settings.get("metamodel") is not None else 'PCK'
        self.p_max = settings["max_pol"] if settings.get("max_pol") is not None else 10

        # if self.metamodeltype is 'PCK': assert "max_pol" in settings and  \
        #         'Missing max polynomial degree, max_pol'
       
        # self.date_record = datetime.now().strftime("%Y_%m_%d_%H%M%S")

        #----------------------------------------------------------------------
        
        self.doe_limits = {}              #dict with marginals limits (for initial sampling)
        self.config, poly_type = {}, []
        self.marginals = settings["marginals"][0]
        
        self.dim = len(settings["marginals"][0])      #Input dimensions
        self.passive_strat = settings["passive_sampling"][1]
        self.n_o = settings["passive_sampling"][0]

        # self.MaternCoef = settings["model"][0] 
        self.MaternCoef = 5/2

        # self.ModelType = model["metamodel"]

        self.targetF = settings["active_sampling"][0]
        self.learningF = settings["active_sampling"][1]
        #----------------------------------------------------------------------
        for margin in range (0, self.dim):
            self.doe_limits['x' + str (margin + 1)] = [-1, 1]  
            poly_type.append(settings["marginals"][0]['x' + str (margin + 1)][2])

        self.config["pol_type"] = poly_type

        #----------------------------------------------------------------------
        if self.metamodeltype == 'PCK':
            self.surrogate = PC_Kriging(self.config)

        elif self.metamodeltype == 'GP':
            kernel_matern = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0), nu=self.MaternCoef)
            self.surrogate = GaussianProcessRegressor(kernel=kernel_matern, random_state=0)           
            
    def passive_sampling (self, n):
        
        xn_o = np.zeros((n, self.dim))      #normalized samples
        xr_o = np.zeros((n, self.dim))      #scaled samples
        yn_o = np.zeros((n))           #observations
        #----------------------------------------------------------------------
        if  self.passive_strat == 'LHS':   #Latin hypercube sampling

            Xdoe_N = build.space_filling_lhs( self.doe_limits , num_samples = n )
        
        # else WE SHOULD ADD 'RND' in case of random sampling on the marginals  
        
        #----------------------------------------------------------------------
        # evaluation and isotransformation

        for margin in range (0, self.dim):
            
            xn_o[:, margin] = Xdoe_N['x' + str (margin + 1)]
            
            if self.config["pol_type"][margin] == 'hermite':
                mean_ = self.marginals['x' + str (margin + 1)][0]
                var_ = self.marginals['x' + str (margin + 1)][1]
                xr_o[:, margin] = self.scalehermite( xn_o[:, margin], mean_, var_ )
            else:
                min_ = self.marginals['x' + str (margin + 1)][0]
                max_ = self.marginals['x' + str (margin + 1)][1]
                xr_o[:, margin] = self.scalelegendre( xn_o[:, margin], min_, max_ )

        yn_o = self.targetF(xr_o)
        
        return xn_o, xr_o, yn_o  #returns initial sampling (normalized, scaled and R_evaluations)
        

    def active_sampling(self, n_exp, n_act, MCpool, GroundT):
    
        for experiments in range(n_exp):
            
            print('Experiment: ', experiments+1 , '########################################################' )
            
            self.training_results = {}   #results file
            
            xn, xr, yn = self.passive_sampling(self.n_o)
            
            for points in range(n_act + 1):
                
                print('Training model with', len(xn),'samples..')
                mean_loo = np.zeros(len(xn))
                var_loo = np.zeros(len(xn))
                
                if self.metamodeltype == 'PCK':
                    ModelName = 'PCK_' + str(len(xn))

                    opt_length_it = np.zeros(self.p_max - 1)
                    eloo_results = np.zeros(self.p_max - 1)
                    sumat = np.zeros((len(xn),(self.p_max - 1)))
                    # optimal polymonial degree search -----------------------------------
                
                    for p in range(1, self.p_max):
                        
                        d = self.surrogate.distance(xn, xn)
                        lmin = np.min(d[d!=0])
                        
                        ModelParam_temp = self.surrogate.train(xn, yn, p, lmin, self.MaternCoef)
                        
                        opt_length = self.surrogate.optimize('shgo')
                        
                        #getting LOO-CV error--------------------------------
                        
                        e_loo, var_loo = self.eloo(xn, yn, p, opt_length)
                        
                        mean_eloo = np.mean(e_loo)
            
                        sumat[:,p-1] = np.divide(e_loo, var_loo)
                        
                        eloo_results[p-1] = mean_eloo
                        opt_length_it[p-1] = opt_length 

                    ## training optimal PCK model ----------------------------

                    opt = np.argmin(eloo_results)    #selected based 'mean eloo'
                    opt_length = opt_length_it[opt]
                    opt_degree = int(opt+1)

                    ModelParam = self.surrogate.train (xn, yn, opt_degree, opt_length , self.MaternCoef)
                
                elif self.metamodeltype == 'GP':
                    ModelName = 'GP_' + str(len(xr))
                    ## training GP model ----------------------------
                    self.surrogate.fit(xr, yn)
                
                MCinputs_norm = np.zeros((int(MCpool), self.dim))
                MCinputs = np.zeros((int(MCpool), self.dim))
                LOOCV = np.zeros(int(MCpool))
                
                ## Generating pool of samples - MCS -----------------------------------
                for margin in range (0, self.dim):

                    if self.config["pol_type"][margin] == 'hermite':
                        
                        MCinputs_norm[:, margin] = np.random.normal(0, 1, size=int(MCpool))
                        
                        mean_ = self.marginals['x' + str (margin + 1)][0]
                        var_ = self.marginals['x' + str (margin + 1)][1]
                        
                        MCinputs[:, margin] = self.scalehermite( MCinputs_norm[:, margin], mean_, var_ )
                        
                    else:
                        
                        MCinputs_norm[:, margin] = np.random.uniform(-1, 1, size=int(MCpool))
                        
                        min_ = self.marginals['x' + str (margin + 1)][0]
                        max_ = self.marginals['x' + str (margin + 1)][1]
                        
                        MCinputs[:, margin] = self.scalelegendre( MCinputs_norm[:, margin], min_, max_ )
                
                # Pf estimation ----------------------------------------------
                print('Model predictions with MC population...')
                if self.metamodeltype == 'PCK':
                    meanMC, stdMC = self.surrogate.predict(MCinputs_norm)    # mean, std
                    
                elif self.metamodeltype == 'GP':
                    meanMC, stdMC = self.surrogate.predict(MCinputs, return_std=True)

                fail_samples_SUMO = np.sum(np.asarray(meanMC) < 0 )
                Pf_SUMO = fail_samples_SUMO / MCpool
                
                if GroundT == True :
                    yMC_ref = self.targetF(MCinputs)  
                    fail_samples_ref = np.sum(yMC_ref < 0 )
                    Pf_Ref = fail_samples_ref / MCpool
                    
                else:
                    Pf_Ref = 'None'
                
                cov_pf = np.sqrt((1 - Pf_SUMO ) / (Pf_SUMO * MCpool) )

                B_sumo = - norm.ppf( Pf_SUMO ) 

                if self.metamodeltype == 'PCK':
                    print('Degree:', opt_degree, 'e_LOO:', "%.4f" % round(np.min(eloo_results), 4), 'Pf_ref:',
                        Pf_Ref ,'Pf_SuMo:', Pf_SUMO , 'B_sumo:', "%.3f" % round(B_sumo, 4), 'CoV_SuMo:', "%.3f" % round(cov_pf, 4))
                    self.training_results[ModelName] = Pf_SUMO, cov_pf, opt_length, opt_degree, np.min(eloo_results)
                    filename = 'PCK_Batch_'+ str(experiments+1)+'.sav'
                    pickle.dump(self.training_results, open(filename, 'wb'))

                elif self.metamodeltype == 'GP':
                    print('Pf_ref:', Pf_Ref ,'Pf_SuMo:', Pf_SUMO , 'B_sumo:', "%.3f" % round(B_sumo, 4), 'CoV_SuMo:', "%.3f" % round(cov_pf, 4))
                    self.training_results[ModelName] = Pf_SUMO
                    filename = 'GP_Batch_'+ str(experiments+1)+'.sav'
                    pickle.dump(self.training_results, open(filename, 'wb'))

                print(' ')
                            
                # Learning Function ----------------------------
        
                if self.learningF == 'U' :
                    U_f = self.U_function(meanMC.reshape(-1), stdMC.reshape(-1))
                    xr_new = MCinputs[np.argmin(U_f)]
                    xn_new = MCinputs_norm[np.argmin(U_f)]
                    
                elif self.learningF == 'EFF':
                    eff = self.EFF(meanMC.reshape(-1),stdMC.reshape(-1), 0)
                    xr_new = MCinputs[np.argmax(eff)]
                    xn_new = MCinputs_norm[np.argmax(eff)]
                
                # LOO CV errors ###################################################
                #variance modification based LOO CV erros around voronoi cells
                
                elif self.learningF == 'ULOO' :
                    for k in range (0, MCpool):               
                        voro = self.VoronoiCell(MCinputs[k], xr)
                        LOOCV[k]= stdMC[k] * (1 + sumat[voro, opt])
                    
                    U_f = self.U_function(meanMC.reshape(-1), LOOCV.reshape(-1))
                    xr_new = MCinputs[np.argmin(U_f)]
                    xn_new = MCinputs_norm[np.argmin(U_f)]
                    
                elif self.learningF == 'EFFLOO' :
                    for k in range (0, MCpool):               
                        voro = self.VoronoiCell(MCinputs[k], xr)
                        LOOCV[k]= stdMC[k]*(1 + sumat[voro, opt])

                    eff = self.EFF(meanMC.reshape(-1), LOOCV.reshape(-1), 0)
                    xr_new = MCinputs[np.argmax(eff)]
                    xn_new = MCinputs_norm[np.argmax(eff)]
                
                else:
                    print ('No learning metric speficied')
                
                xr = np.concatenate((xr, xr_new.reshape(-1, self.dim)), axis=0)
                xn = np.concatenate((xn, xn_new.reshape(-1, self.dim)), axis=0)

                y_new = self.targetF(xr_new.reshape(-1, self.dim))
                yn = np.concatenate((yn, y_new), axis=0)

            Model_last = 'last'
            self.training_results[Model_last] = self.surrogate, xn, xr, yn 
            pickle.dump(self.training_results, open(filename, 'wb'))
                
    def eloo (self, xn, yn, p, length):   #LOO CV squared errors
    
        mean_loo = np.zeros(len(xn))
        std_loo = np.zeros(len(xn))
        
        if self.metamodeltype == 'PCK':

            self.surrogate_loo = PC_Kriging(self.config)    # for LOOCV with same 'config' as specified in the original model

        for out in range (0, len(xn)):

            out_el = [ (out* self.dim + i) for i in range (self.dim)]

            yn_loo = np.delete(yn, out)                             #y_n-i      leaving element i out the observations 
            xn_loo = np.delete(xn, out_el).reshape(-1, self.dim)    #x_n-i     leaving element i out the nomalized inputs (xn)

            #training LOO
            modelpar_loo = self.surrogate_loo.train(xn_loo , yn_loo , p , length, self.MaternCoef)

            #predicting LOO over each removed sample

            mean_loo[out], std_loo[out] = self.surrogate_loo.predict(xn[out].reshape(1, -1))
        
        eloo = (yn - mean_loo)**2  #vector of LOO CV Squared errors
        
        return eloo, std_loo


    def EFF(self, u, v, z):
        zl=-2*v
        zh=2*v
        return ((u-z)*( 2*norm.cdf((z-u)/v) - norm.cdf((zl-u)/v) - norm.cdf((zh-u)/v)) 
                -(v)*( 2*norm.pdf((z-u)/v) - norm.pdf((zl-u)/v) - norm.pdf((zh-u)/v))  
                +(2*v)*(norm.cdf((zh-u)/v) - norm.cdf((zl-u)/v)))
    
    def U_function(self, u, v):
        return np.abs(u)/v
    
    def VoronoiCell(self, x,xn):   #given x [single value] return the index of the closest xn [array]
        dist = self.surrogate.distance(x.reshape(1,-1),xn)
        return np.argmin(dist)
    
    def scalelegendre(self, x, new_min, new_max):
        #to go from (-1,1) to (new_min, new_max) 
        return ((new_min+new_max)+((new_max-new_min)*x))/2

    def scalehermite(self, x , mean, sigma):
        return mean+sigma*x