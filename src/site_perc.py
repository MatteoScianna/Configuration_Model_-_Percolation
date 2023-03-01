#%%
prob_vec = prob_list_poisson(4,1) #customizable values of lambda and probability distribution
degree_list = discrete_samples(prob_vec,10000) #customizable number of nodes in the network
graph = config_model(degree_list)
#%%
size_biggest_component, k, k_s = get_S_synth_site(degree_list, 1) #site percolation
S,k1,k_s1 = get_S_theory_site(prob_vec)

#%%
p = np.linspace(0,1,100)
plt.scatter(p,[s for s in size_biggest_component], color = "blue", facecolor = "none")
plt.plot(np.linspace(0,1,10),S, color = "red")
plt.axvline(x=k/(k_s-k), color = "blue", linestyle = "dotted")
plt.axvline(x=k1/(k_s1-k1), color = "red", linestyle = "dotted")
plt.xlabel(r'$\phi$')
plt.ylabel("Size of Giant Cluster S")
plt.legend([r'Simulation',r'Theory',r'$\phi_{c}$ - Simulation',r'$\phi_{c}$ - Theory'])
plt.title("Site Percolation on ER Network" +"\n"+ "Poisson Distribution with <k> = 4")
plt.show()

#%%
