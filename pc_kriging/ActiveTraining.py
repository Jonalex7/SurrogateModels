import numpy as np
from pc_kriging import PC_Kriging
from datetime import datetime
from scipy.stats import norm
from doepy import build
from scipy import optimize

class Active_Training():

    def __init__(self, model=None, settings=None):
        if model is None:
            model = {"metamodel" : 'PCK'}
        assert "metamodel" in model and \
            "Missing model type"
        
        self.date_record = datetime.now().strftime("%Y_%m_%d_%H%M%S")

        #----------------------------------------------------------------------
        self.doe_limits = {}              #dict with marginals limits (for initial sampling)
        self.config, poly_type = {}, []
        self.marginals = settings["marginals_R"][0]
        
        self.dim = len(settings["marginals_R"][0])      #Input dimensions
        self.passive_strat = settings["passive_sampling"][1]
        self.n_o = settings["passive_sampling"][0]

        self.MaternCoef = settings["model"][0]
        self.p_max = settings["model"][1]
        self.ModelType = model["metamodel"]
        self.targetF = settings["active_sampling"][0]
        self.learningF = settings["active_sampling"][1]
        #----------------------------------------------------------------------
        for margin in range (0, self.dim):
            self.doe_limits['x' + str (margin + 1)] = [-1, 1]  
            poly_type.append(settings["marginals_R"][0]['x' + str (margin + 1)][2])

        self.config["pol_type"] = poly_type

        #----------------------------------------------------------------------
        if self.ModelType == 'PCK':
            
            self.surrogate = PC_Kriging(self.config)
            
        # else we may add 'PCE' or 'Kriging' later
            

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
                xr_o[:, margin] = self.surrogate.scalehermite( xn_o[:, margin], mean_, var_ )
            else:
                min_ = self.marginals['x' + str (margin + 1)][0]
                max_ = self.marginals['x' + str (margin + 1)][1]
                xr_o[:, margin] = self.surrogate.scalelegendre( xn_o[:, margin], min_, max_ )
        
        yn_o = self.targetF(xr_o)
        
        return xn_o, xr_o, yn_o  #returns initial sampling (normalized, scaled and R_evaluations)
        

    def active_sampling(self, n_act, MCpool, GroundT):
        
        self.training_results = {}   #results file { Surrogate , Pf, CoV_Pf , eLoo , mse }
        
        xn, xr, yn = self.passive_sampling(self.n_o)
        
        for points in range(n_act + 1):
            
            print('N° Training Samples: ', len(xn))

            mse_results = np.zeros(self.p_max - 1)
            opt_length_it = np.zeros(self.p_max - 1)
            eloo_results = np.zeros(self.p_max - 1)

            mean_loo = np.zeros(len(xn))
            var_loo = np.zeros(len(xn))
            sumat = np.zeros((len(xn),(self.p_max - 1)))

            # OPTIMAL SURROGATE MODEL -----------------------------------
        
            for p in range(1, self.p_max):
                
                d = self.surrogate.distance(xn, xn)
                lmin = np.min(d[d!=0])
                
                ModelParam_temp = self.surrogate.train(xn, yn, p, lmin, self.MaternCoef)
                
                opt_length = self.surrogate.optimize('shgo')
                
                #getting LOO-CV error
                
                e_loo, var_loo = self.eloo(xn, yn, p, opt_length)
                
                mean_eloo = np.mean(e_loo)
    
                sumat[:,p-1] = np.divide(e_loo, var_loo)
                
                eloo_results[p-1] = mean_eloo
                opt_length_it[p-1] = opt_length 
                
            #--------------------------------------- Gen. error over a set of test points
#                 mean0, var0 = self.surrogate.predict(XN)    # test points predictions mean, variance
#                 mse = np.mean ((YN - mean0)**2)
#                 mse_results[p-1] = mse

            ## training optimal model ----------------------------

            opt = np.argmin(eloo_results)    #selected based 'mean eloo'
            opt_length = opt_length_it[opt]
            opt_degree = int(opt+1)

            ModelParam = self.surrogate.train (xn, yn, opt_degree, opt_length , self.MaternCoef)
            
            MCinputs_norm = np.zeros((int(MCpool), self.dim))
            MCinputs = np.zeros((int(MCpool), self.dim))
            LOOCV = np.zeros(int(MCpool))
            
            ## Generating pool of samples - MCS -----------------------------------

            for margin in range (0, self.dim):

                if self.config["pol_type"][margin] == 'hermite':
                    
                    MCinputs_norm[:, margin] = np.random.normal(0, 1, size=int(MCpool))
                    
                    mean_ = self.marginals['x' + str (margin + 1)][0]
                    var_ = self.marginals['x' + str (margin + 1)][1]
                    
                    MCinputs[:, margin] = self.surrogate.scalehermite( MCinputs_norm[:, margin], mean_, var_ )
                    
                else:
                    
                    MCinputs_norm[:, margin] = np.random.uniform(-1, 1, size=int(MCpool))
                    
                    min_ = self.marginals['x' + str (margin + 1)][0]
                    max_ = self.marginals['x' + str (margin + 1)][1]
                    
                    MCinputs[:, margin] = self.surrogate.scalelegendre( MCinputs_norm[:, margin], min_, max_ )
            
            # Pf estimation ----------------------------------------------
            
            meanMC, varMC = self.surrogate.predict(MCinputs_norm)    # mean, variance
            fail_samples_SUMO = np.sum(np.asarray(meanMC) < 0 )
            Pf_SUMO = fail_samples_SUMO / MCpool
            
            if GroundT == True :
                yMC_ref = self.targetF(MCinputs)  
                fail_samples_ref = np.sum(yMC_ref < 0 )
                Pf_Ref = fail_samples_ref / MCpool
                
            else:
                Pf_Ref = 'None'
            
            cov_pf = np.sqrt((1 - Pf_SUMO ) / (Pf_SUMO * MCpool) )

            print('Degree', opt_degree, 'e_LOO', "%.5f" % round(np.min(eloo_results), 4), 'Pf_ref',
                  Pf_Ref ,'Pf_SuMo', Pf_SUMO , 'CoV_SuMo', "%.5f" % round(cov_pf, 4))
            
            # Saving results ----------------------------

            self.training_results[str(len(xn))+'_points'] = opt_degree, opt_length, np.min(eloo_results), Pf_SUMO, cov_pf 
            
            # Learning Function ----------------------------
    
            if self.learningF == 'U' :
                U_f = self.U_function(meanMC.reshape(-1), varMC.reshape(-1))
                xr = np.append(xr, MCinputs[np.argmin(U_f)]).reshape(-1, self.dim)
                xn = np.append(xn, MCinputs_norm[np.argmin(U_f)]).reshape(-1, self.dim)
                
            elif self.learningF == 'EFF':
                eff = self.EFF(meanMC.reshape(-1),varMC.reshape(-1), 0)
                xr = np.append(xr, MCinputs[np.argmax(eff)]).reshape(-1, self.dim)
                xn = np.append(xn, MCinputs_norm[np.argmax(eff)]).reshape(-1, self.dim)
            
            # LOO CV errors ###################################################
            #variance modification based LOO CV erros around voronoi cells
            
            elif self.learningF == 'ULOO' :
                for k in range (0, MCpool):               
                    voro = self.VoronoiCell(MCinputs[k], xr)
                    LOOCV[k]= varMC[k] * (1 + sumat[voro, opt])
                
                U_f = self.U_function(meanMC.reshape(-1), LOOCV.reshape(-1))
                xr = np.append(xr, MCinputs[np.argmin(U_f)]).reshape(-1, self.dim)
                xn = np.append(xn, MCinputs_norm[np.argmin(U_f)]).reshape(-1, self.dim)
                
            elif self.learningF == 'EFFLOO' :
                for k in range (0, MCpool):               
                    voro = self.VoronoiCell(MCinputs[k], xr)
                    LOOCV[k]= varMC[k]*(1 + sumat[voro, opt])
                
                eff = self.EFF(meanMC.reshape(-1), LOOCV.reshape(-1), 0)
                xr = np.append(xr, MCinputs[np.argmax(eff)]).reshape(-1, self.dim)
                xn = np.append(xn, MCinputs_norm[np.argmax(eff)]).reshape(-1, self.dim)
            
            else:
                print ('No learning metric speficied')
            
            yn = self.targetF(xr)
            
                
    def eloo (self, xn, yn, p, length):   #LOO CV squared errors
    
        mean_loo = np.zeros(len(xn))
        var_loo = np.zeros(len(xn))
        
        if self.ModelType == 'PCK':

            self.surrogate_loo = PC_Kriging(self.config)    # for LOOCV with same 'config' as specified in the original model

        for out in range (0, len(xn)):

            out_el = [ (out* self.dim + i) for i in range (self.dim)]

            yn_loo = np.delete(yn, out)                             #y_n-i      leaving element i out the observations 
            xn_loo = np.delete(xn, out_el).reshape(-1, self.dim)    #x_n-i     leaving element i out the nomalized inputs (xn)

            #training LOO
            modelpar_loo = self.surrogate_loo.train(xn_loo , yn_loo , p , length, self.MaternCoef)

            #predicting LOO over each removed sample

            mean_loo[out], var_loo[out] = self.surrogate_loo.predict(xn[out].reshape(1, -1))
        
        eloo = (yn - mean_loo)**2  #vector of LOO CV Squared errors
        
        return eloo, var_loo


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