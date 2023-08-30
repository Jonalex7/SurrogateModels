#Parabolic/Concave limit-state function
#ref. (HMCMC High Dim, Prof. Kostas) PF= 5.90E-8
# 2 independent standard normal random variables
def example_2(x):   
    r = 6.0
    k = 0.3
    e = 0.1
    g = r - (x[:,1]) - k*(x[:,0] - e)**2
    return g   
   
#4: The Himmelblau Function
#ref. (HMCMC High Dim, Prof. Kostas)
# 2 independent standard normal random variables
def example_4 (x):
    beta = 95   # (beta = 95 for Ref. PF=1.65E-4)   (beta = 50 for Ref. PF=2.77E-7)
    term1 = (((0.75*x[:,0] - 0.5)**2 / 1.81) + ((0.75*x[:,1] - 0.5) / 1.81) - 11)**2
    term2 = (((0.75*x[:,0] - 1.0)/ 1.81) + ((0.75*x[:,1] - 0.5)**2 / 1.81) - 7)**2
    g = term1 + term2 - beta
    return g