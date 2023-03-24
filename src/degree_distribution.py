#%%
import numpy as np
import scipy as sp
from scipy.stats import poisson, uniform, powerlaw, expon
from functions import prob_list_poisson,prob_list_powerlaw,discrete_inverse_trans,discrete_samples,config_model


#Example of definition of degree distribution 
#and configuration model
#for a Poisson distribution with mean 4
#%%
n = 1000
k = 4
prob_vec = prob_list_poisson(k,n)
#n=1000
#gamma = 2.5
#min = 3
#prob_vec = prob_list_powerlaw(gamma,min,n)
degree_list = discrete_samples(prob_vec,n)
graph = config_model(degree_list)

# %%
degree_list
# %%
