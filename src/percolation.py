
#THIS SCRIPT PERFORMS DIFFERENT PERCOLATION STRATEGIES#
#%%
import random as rand
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.stats import poisson, uniform, powerlaw, expon 
from scipy.optimize import fsolve
import pandas as pd
import seaborn as sns
import collections  

#%%
#functions to perform percolation on a graph
from functions import remove_edges_uniform,remove_nodes_degree,remove_nodes_uniform 
#functions to obtain evolution of size of giant cluster on a graph
from functions import get_S_synth_site,get_S_synth_bond,get_S_synth_non_uniform
#functions to obtain theoretical evolution of size of giant cluster on a graph
from functions import get_S_theory_site,get_S_theory_bond,get_S_theory_non_uniform
#functions to solve the equation of the critical point
from functions import g_0,g_1,f_0,f_1
#function to obtain the size of the giant cluster 
from functions import get_size_biggest_component
#funtions to obtain critical point
from functions import get_k_s,get_k


#%%
#initial definition of the graph
from functions import prob_list_poisson
from functions import discrete_samples
from functions import discrete_inverse_trans
from functions import config_model
mean_deg = 4
n = 1000
prob_vec = prob_list_poisson(mean_deg,n)
degree_list = discrete_samples(prob_vec,n)
graph = config_model(degree_list)


#%%
#SITE PERCOLATION
size_biggest_component, k, k_s = get_S_synth_site(degree_list, 50) #50 is the number of realizations
S,k1,k_s1 = get_S_theory_site(prob_vec) 
p = np.linspace(0,1,50)
plt.scatter(p,[s for s in size_biggest_component], color = "blue", facecolor = "none")
plt.plot(np.linspace(0,1,50),S, color = "red")
plt.axvline(x=k/(k_s-k), color = "blue", linestyle = "dotted")
plt.axvline(x=k1/(k_s1-k1), color = "red", linestyle = "dotted")
plt.xlabel(r'$\phi$')
plt.ylabel("Size of Giant Cluster S")
plt.legend([r'Simulation',r'Theory',r'$\phi_{c}$ - Simulation',r'$\phi_{c}$ - Theory'])
plt.title("Site Percolation on ER Network" +"\n"+ "Poisson Distribution with <k> = 4")
plt.show()

#%%
#BOND PERCOLATION
size_biggest_component, k, k_s = get_S_synth_bond(degree_list, 50) #50 is the number of realizations
S,k1,k_s1 = get_S_theory_bond(prob_vec)
p = np.linspace(0,1,50)
plt.scatter(p,[s for s in size_biggest_component], color = "blue", facecolor = "none")
plt.plot(np.linspace(0,1,len(S)),S, color = "red")
plt.axvline(x=k/(k_s-k), color = "blue", linestyle = "dotted")
plt.axvline(x=k1/(k_s1-k1), color = "red", linestyle = "dotted")
plt.xlabel(r'$\phi$')
plt.ylabel("Size of Giant Cluster S")
plt.legend([r'Simulation',r'Theory',r'$\phi_{c}$ - Simulation',r'$\phi_{c}$ - Theory'])
plt.title("Bond Percolation on ER Network" +"\n"+ "Poisson Distribution with <k> = 4")
plt.show()

# %%
# NON UNIFORM PERCOLATION

S = get_S_theory_non_uniform(degree_list,prob_vec,mean_deg)
components = get_S_synth_non_uniform(graph)

phi = [i for i in range(max(degree_list))]
plt.scatter(phi,[s for s in components], color = "blue", facecolor = "none")
plt.plot(phi,S, color = "red")
plt.xlabel("Maximum degree $k_{0}$")
plt.ylabel("Size of Giant Cluster S")
plt.legend([r'Simulation',r'Theory',r'$\phi_{c}$ - Simulation',r'$\phi_{c}$ - Theory'])
plt.title("Non uniform Percolation on ER Network" +"\n"+ "Poisson distribution with " r'k = 4')
plt.show()