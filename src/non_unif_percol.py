prob_vec = prob_list_powerlaw(2.5,2)
degree_list = discrete_samples(prob_vec,10000)
mean_deg = np.mean(degree_list)
graph = config_model(degree_list)

S = get_S_theory_non_uniform(degree_list,prob_vec)
components = get_S_synth_non_uniform(graph)

phi = [i for i in range(max(degree_list))]
plt.scatter(phi,[s for s in components], color = "blue", facecolor = "none")
plt.plot(phi,S, color = "red")
plt.xlabel("Maximum degree $k_{0}$")
plt.ylabel("Size of Giant Cluster S")
plt.legend([r'Simulation',r'Theory',r'$\phi_{c}$ - Simulation',r'$\phi_{c}$ - Theory'])
plt.title("Non uniform Percolation on SF Network" +"\n"+ "Power Law Distribution with " r'$\alpha$ = 2.5')
plt.show()
