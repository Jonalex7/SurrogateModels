import numpy as np
#Example 2: Parabolic/Concave limit-state function
#ref. (HMCMC High Dim, Prof. Kostas) PF= 5.90E-8
# 2 independent standard normal random variables
def example_2(x):   
    r = 6.0
    k = 0.3
    e = 0.1
    g = r - (x[:,1]) - k*(x[:,0] - e)**2
    return g   
   
#Example 4: The Himmelblau Function
#ref. (HMCMC High Dim, Prof. Kostas)
# 2 independent standard normal random variables
def example_4 (x):
    beta = 95   # (beta = 95 for Ref. PF=1.65E-4)   (beta = 50 for Ref. PF=2.77E-7)
    term1 = (((0.75*x[:,0] - 0.5)**2 / 1.81) + ((0.75*x[:,1] - 0.5) / 1.81) - 11)**2
    term2 = (((0.75*x[:,0] - 1.0)/ 1.81) + ((0.75*x[:,1] - 0.5)**2 / 1.81) - 7)**2
    g = term1 + term2 - beta
    return g

# Example 8: High-dimensional highly nonlinear problem - 100 D
#ref. (HMCMC High Dim, Prof. Kostas)
# 100 independent standard normal random variables

def example_8(X):
    y_0 = 2.5   # 2.5 for PF_ref = 3.40E-5  //  3.5 for PF_ref = 7.96E-7 // 4.5 for PF_ref = 6.75E-9
    d = len(X[0])
    g = np.zeros((len(X)))

    j_ind = np.arange(1,10)
    k_ind = np.arange(11,14)
    l_ind = np.arange(15,17)

    for sample in range(0, len(X)):
        sum_i = np.sum(X[sample])
        sum_j = X[sample][j_ind].sum()
        sum_k = X[sample][k_ind].sum()
        sum_l = X[sample][l_ind].sum()

        g[sample] =  y_0 - ((1/np.sqrt(d))*(sum_i)) + 2.5*(X[sample][0] - sum_j)**2 + (X[sample][10] - sum_k)**4  + (X[sample][14] - sum_l)**8
    
    return g