#%%
from functions import *
#%%
#definition of the model
t = 1000
prob_vec = prob_list_poisson(8,t)
degree_list = discrete_samples(prob_vec,t)
alpha = 0.3
graph = config_model(degree_list)
#%%
#compute SATRW theory
p,q,d = prob_matrix(prob_vec,alpha,1)
s = SATRW_theory(p,q)
#%%
#compute SATRW synthetic over 35 simulations
nsim = 35
S = np.zeros(t-1)
P = np.zeros((t,max(degree_list)))
for i in range(nsim):
    s_1,p_1 = SATRW_synth(graph,degree_list,1,alpha)
    S+=np.array(s_1)
    print(i+1)
    P+=p_1
S = S/nsim
P = P/nsim
#%%
#plot probability distribution
colors = plt.cm.inferno(np.linspace(0,1,15))
for i in range(0,t,int(t/10)):
    plt.plot(p[i],color = colors[i//int(t/10)])
    plt.plot(P[i],color = colors[i//int(t/10)],linestyle = '--')
    #set xlim 
    plt.xlim(0,15)
    plt.ylim(0,0.5)
plt.legend(['theory','simul'],loc = 'upper right')
plt.title('Degree distribution - SF graph, $\\alpha$ = 0.3') 
#save the figure
#plt.savefig('degree_dist_SF.png',transparent = False,facecolor = "white", dpi = 200, bbox_inches = "tight")
# %%
#plot giant component size
#neg values of s to 0
s = [0 if i < 0 else i for i in s]
T = np.linspace(0,1,t-10)
t1 = np.linspace(0,t-3,15)
plt.plot(T,s[:990],color = "blue")
plt.scatter([t_1/1000 for t_1 in t1] ,[S[int(i)] for i in t1], color = 'blue', marker = 'o')
plt.legend(['theory','simul'],loc = 'upper right')
plt.title('Giant component size - ER graph, $\\alpha$ = 0.3')
#plt.savefig('giant_comp_ER.png',transparent = False,facecolor = "white", dpi = 200, bbox_inches = "tight")
# %%
