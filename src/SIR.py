#%%
from functions import *
from scipy.integrate import odeint


#%%

def SIR(graph, beta, gamma):
    
    S = []
    I = []
    R = []
    S.append(len(graph.nodes)-1)
    I.append(1)
    R.append(0)
    list_nodes = list(graph.nodes)
    #choose a random node to be infected
    infected = [rand.choice(list_nodes)]
    #choose a random node to be recovered
    recovered = []
    #choose a random node to be susceptible
    susceptible = [list_nodes[i] for i in range(len(list_nodes)) if i not in infected]
    #iterate
    #while no infected
    t = 0
    while len(infected)>0:
        print("time: "+str(t))
        #add new infected nodes
        inf1 = infected.copy() #to avoid changing the list while iterating
        for node in inf1:
            for neighbor in graph.neighbors(node):
                if neighbor in susceptible:
                    if rand.random() < beta:
                        infected.append(neighbor)
                        susceptible.remove(neighbor)
        for node in inf1:
            if rand.random() < gamma:
                recovered.append(node)
                infected.remove(node)
        #update lists
        S.append(len(susceptible))
        I.append(len(infected))
        R.append(len(recovered))
        print("Susceptible: "+str(len(susceptible)))
        print("Infected: "+str(len(infected)))
        print("Recovered: "+str(len(recovered)))
        t+=1
    return S,I,R,t

def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

def SIR1(beta,tau,graph):
    S = []
    I = []
    R = []
    S.append(len(graph.nodes)-1)
    I.append(1)
    R.append(0)
    list_nodes = list(graph.nodes)
    if len(list_nodes) == 0:
        return S,I,R
    #choose a random node to be infected
    patient0 = rand.choice(list_nodes)
    infected = [patient0]
    #choose a random node to be recovered
    recovered = []
    #choose a random node to be susceptible
    susceptible = [i for i in list_nodes if i not in infected]
    dict_recovering_nodes = {patient0:0}
    t=0
    while len(infected) > 0:
        #print("time: "+str(t))
        #add new infected nodes
        inf1 = infected.copy() #to avoid changing the list while iterating
        for node in inf1:
            for neighbor in graph.neighbors(node):
                if neighbor in susceptible:
                    k = rand.random()
                    if k < beta:
                        infected.append(neighbor)
                        susceptible.remove(neighbor)
                        dict_recovering_nodes[neighbor] = 0
        for node in infected:
            #print("node: "+str(node))
            if dict_recovering_nodes[node] == tau:
                recovered.append(node)
                infected.remove(node)
                del dict_recovering_nodes[node]
            else:
                dict_recovering_nodes[node] += 1
                #print("dict_recovering_nodes[node]: "+str(dict_recovering_nodes[node]))
        t += 1
        #update lists
        S.append(len(susceptible))
        I.append(len(infected))
        R.append(len(recovered))
        #print("Susceptible: "+str(len(susceptible)))
        #print("Infected: "+str(len(infected)))
        #print("Recovered: "+str(len(recovered)))
    return S,I,R

#%%
#Evolution SIR Model (Theory)
N = 10000
I0, R0 = N*0.01, 0
S0 = N*0.99
beta, gamma = 1, 0.4
days = 30
t = np.linspace(0,days,days)
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
S_t, I_t, R_t = odeint(deriv, y0, t, args=(N, beta, gamma)).T
r = sp.optimize.fsolve(lambda x: 1-(S0/N)*np.exp(-(beta*x)/gamma)-x,0.5)

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S_t/N, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I_t/N, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R_t/N, 'g', alpha=0.5, lw=2, label='Recovered')
#add horizontal dotted line for r
ax.plot([0,days],[r,r], 'k--', alpha=0.5, lw=2, label='total outbreak size',color = "green")
#add y value for r
ax.text(days,r+0.01, str(round(r[0],2)), color = "green")
ax.set_xlabel('Time /days')
ax.set_ylabel('Number ('+str(N)+')')
ax.set_ylim(0,1.1)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.savefig('SIR_theory.png', transparent = True, dpi = 200, bbox_inches = "tight")

prob_vec = prob_list_poisson(0.5, 7)

# %%
prob_vec = prob_list_poisson(7)
#Percolation vs Recovery
T = np.linspace(0, 1, 100)
R = []
U = []
for t in T:
    #solve equation for x = g_1(1+(x-1)*t)
    u = sp.optimize.fsolve(lambda x: g_1(1+(x-1)*t,prob_vec,4)-x,0.5)[0]
    U.append(u)
    r = 1-g_0(1+(u-1)*t,prob_vec)
    R.append(r)
#%%
PHI = np.linspace(0, 1, 100)
S = []
for phi in PHI:
    u = sp.optimize.fsolve(lambda x: 1-phi + phi*g_1(x,prob_vec,4) - x,0.5)[0]
    s = 1-g_0(u,prob_vec)
    S.append(s)

# %%
plt.plot(T, R, 'r', alpha=0.5, lw=2, label='R')
plt.scatter(PHI,S, color = "blue", label = "S",facecolor = "none")
plt.legend()
#modify legend names
plt.legend(["Total oubreak size - SIR Model","Size giant cluster - Bond Percolation"])
plt.xlabel(r'$\phi$ - T')
plt.ylabel(r'R - S')
#title
plt.title("Total outbreak size vs Size giant cluster - theory")
plt.show()
#save
plt.savefig("SIR_vs_Bond_Percolation_theory.png")

#%%
prob_vec = prob_list_poisson(7)
degree_list = discrete_samples(prob_vec,1000)
graph = config_model(degree_list)
# %%
S_synth = get_S_synth_bond(degree_list, 100)
p = np.linspace(0, 1, 100)
plt.plot(p, S_synth, 'r', alpha=0.5, lw=2, label='S')
# %%
n = 100
BETA = np.linspace(0, 1, n)
R0 = [0]*n
#calculate R0 for 100 realizations of the SIR model
real = 0
while real < 25:
    i = 0
    prob_vec = prob_list_poisson(4)
    degree_list = discrete_samples(prob_vec,1000)
    graph = config_model(degree_list)
    for beta in BETA:
        S,I,R = SIR1(beta,4,graph)
        R0[i] += (R[-1])
        i += 1
        print(i)
    real += 1

# %%
plt.scatter([1-(1-beta)**(4) for beta in BETA], [r0/25000 for r0 in R0],color = "b",facecolor = "none")
plt.plot(p, S_synth, 'r', alpha=0.5, lw=2, label='S')
plt.legend(["Outbreak Size - SIR Model","S - Bond Percolation"])
plt.xlabel(r'T - $\phi$')
plt.ylabel(r'R - S')
plt.title("SIR Model vs Bond Percolation")
#save figure
#plt.savefig("SIR_vs_Bond_Percolation-simulation.png")



# %%
