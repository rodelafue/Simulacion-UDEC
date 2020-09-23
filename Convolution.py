# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 20:29:13 2016

@author: rodrigodelafuente
"""

import numpy as np
import scipy.special
from matplotlib import pyplot as plt

"""
# ex-Gaussian distribution
def convolve_exp_norm(alpha, mu, sigma, x):
    co = alpha/2.0 * np.exp( alpha*mu+ alpha*alpha*sigma*sigma/2.0)
    x_erf = (mu + alpha*sigma*sigma - x)/(np.sqrt(2.0)*sigma)
    y = co * np.exp(-alpha*x) * (1.0 - scipy.special.erf(x_erf))
    return y

# input parameters
x = np.arange(-5.0,5.0,0.1)
alpha = 1.0 # index of the exponential function
sigma = 0.5 # dispersion of the normal distribution
mu = -0.5 # center of the normal distribution
exponential = []
for i in range(len(x)):
    value = alpha * np.exp(-alpha*x[i]) if (x[i]>=0.0) else 0.0
    exponential.append(value)
np.array(exponential)
normal = 1.0/(np.sqrt(2.0*np.pi)*sigma) * np.exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))

# convolve
convolution = convolve_exp_norm(alpha, mu, sigma, x)

# plot result
plt.figure(figsize=(7,7))
plt.plot(x,convolution,color='red',linewidth=2.0)
plt.plot(x,exponential,color='blue')
plt.plot(x,normal,color='green')
plt.legend(('Exponentially modified Gaussian','Exponential Distribution (input)','Normal Distribution (input)'),frameon=False)
plt.xlim(-2,5)
plt.text(3.0,0.7,r'$\lambda=1.5$',fontsize=20,verticalalignment='center')
plt.text(3.0,0.6,r'$\mu=-0.5$',fontsize=20,verticalalignment='center')
plt.text(3.0,0.5,r'$\sigma=0.5$',fontsize=20,verticalalignment='center')
plt.savefig('convolution.pdf')
"""

class RandomVariable(object):
    """Parent class for all random variables."""

class Exponential(RandomVariable):
    def __init__(self, lam):
        self.lam = lam

    def generate(self):
        return np.random.exponential(self.lam)

class Uniforme(RandomVariable):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        
    def generate(self):
        return np.random.uniform(self.a,self.b)

class Erlang(RandomVariable):
    def __init__(self, lam, k):
        self.lam = lam
        self.k = k
        self.expo = Exponential(lam)

    def generate(self):
        total = 0
        for i in range(self.k):
            total += self.expo.generate()
        return total

class SumDist(RandomVariable):
    def __init__(self,lista_dist):
        self.list_dist = lista_dist

    def generate(self,n):
        n_list = []
        for i in range(n):
            n_list.append(sum([x.generate() for x in self.list_dist]))
        return n_list

### Plots
lista=[Uniforme(0,1) for i in range(100)]
conv = SumDist(lista)
out = conv.generate(50000)
plt.hist(out)

lista=[Exponential(2) for i in range(2000)]
conv = SumDist(lista)
out = conv.generate(50000)
plt.hist(out)


