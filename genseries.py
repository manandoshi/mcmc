import numpy as np
import pandas as pd
import sys

#arguments are phi, d, sigma, length

length = sys[4]

def genbinomcoeff(n, r):
    
    if r==0:
        return 1
    res = 1
    for i in range(r):
        res = res*(n-i)/(r-i)
    return res

def genSeries(phi, d, sigma):
    
    a = range(length)
    a[0] = np.random.normal(0, sigma, 1)
    
    for i in range(1,length):
        ARI = 0
        for j in range(i):
            if j==i-1:
                ARI += -phi*genbinomcoeff(d, j)*(-1**j)*a[0]
            else:
                ARI -= a[i-j-1]*(-1**(j+1))*(genbinomcoeff(d, j+1) + phi*(genbinomcoeff(d, j)))
        a[i] = np.random.normal(0, sigma, 1) - ARI
    
    return a

phi = float(sys.argv[1])
d = float(sys.argv[2])
sigma = float(sys.argv[3])
series = genSeries(phi, d, sigma)
df = pd.DataFrame(series)
df.to_csv("series_" + str(phi) + "_" + str(d) + "_" + str(sigma) + "_"".csv")