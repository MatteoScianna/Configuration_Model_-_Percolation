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



def get_k(graph):
    list_degrees = []
    for node in graph.nodes:
        list_degrees.append(graph.degree(node))
    k = sum(list_degrees)/len(graph.nodes)
    return k

def get_k_s(graph):
    list_degrees = []
    for node in graph.nodes:
        list_degrees.append(graph.degree(node))
    k_s = sum([el**(2) for el in list_degrees])/len(graph.nodes)
    return k_s


def g_1(x,prob_list,mean_degree):
    sol = 0
    for i in range(0,len(prob_list)-1):
        if prob_list[i] == 0:
            continue
        else:
            sol += ((i+1)*prob_list[i+1]*x**i)/mean_degree
    return sol

def g_0(x,prob_list):
    sol = 0
    for i in range(len(prob_list)-1):
        sol += (prob_list[i]*(x**(i)))
    return sol

def newton_method(f,df,x0 = 0.5, epsilon=0.0001):
    x = x0
    while abs(f(x)) > epsilon:
        x = x - f(x)/df(x)
    return x


def factorial(n):
    if n == 0:
        return 1
    else:
        return n*factorial(n-1)

#%%
#### get probability distributions
def prob_list_poisson(lambda1,n):
    prob_vec = []
    k = 0
    while k<n:
        prob_vec.append((poisson.pmf(k,lambda1)))
        k+=1
    prob_vec = [i/sum(prob_vec) for i in prob_vec]
    return prob_vec

def prob_list_powerlaw(a,min,n):
    prob_vec = [0]*min
    for k in range(min,int(np.sqrt(n))):
        prob_vec.append(k**(-a))
    prob_vec = [i/sum(prob_vec) for i in prob_vec]
    return prob_vec



### get degree lists

def discrete_inverse_trans(prob_vec):
    U=uniform.rvs(size=1)
    if U<=prob_vec[0]:
        return 1
    else:
        for i in range(1,len(prob_vec)+1):
            if sum(prob_vec[0:i])<U and sum(prob_vec[0:i+1])>U:
                return i

def discrete_samples(prob_vec,n):
    sample=[]
    for i in range(0,n):
        sample.append(discrete_inverse_trans(prob_vec))
    #check if sum is even   
    if sum(sample)%2 != 0:
        sample[0] = sample[0]+1
    return sample

##### get configuration model 

def config_model(degree_list):
    #check if sum is even
    if sum(degree_list)%2 != 0:
        degree_list[0] = degree_list[0]+1
    graph = nx.configuration_model(degree_list)
    graph1= graph.copy()
    for (u,v) in graph1.edges():
        if u == v:
            graph.remove_edge(u,v)
    #remove parallel edges
    graph = nx.Graph(graph)
    return graph
    

##### percolation

def remove_edges_uniform(graph,p):
    lista_random = []
    for i in range(len(graph.edges)):
        lista_random.append(rand.random())
    graph1 = graph.copy()
    i = 0
    for edge in graph.edges:
        if lista_random[i] > p:
            graph1.remove_edge(edge[0],edge[1])
        i = i+1
    return graph1


def remove_nodes_degree(graph,p):
    graph1 = graph.copy()
    i = 0
    for node in graph.nodes:
        if graph.degree(node) <p:
            graph1.remove_node(node)
        i = i+1
    return graph1

def remove_nodes_uniform(graph,p):
    lista_random = []
    for i in range(len(graph.nodes)):
        lista_random.append(rand.random())
    graph1 = graph.copy()
    i = 0
    for node in graph.nodes:
        if lista_random[i] > p:
            graph1.remove_node(node)
        i = i+1
    return graph1

##### Get size of biggest component after percolation 

def get_size_biggest_component(g):
    if len(g.nodes) == 0:
        return 0
    else:
        return len(sorted(nx.connected_components(g), key = len, reverse=True)[0])
    
def f_0(z,prob_vec,max_phi):
    sol = 0
    for k in range(max_phi):
        sol = sol + prob_vec[k]*(z**k)
    return sol

def f_1(z,prob_vec,max_phi,mean_deg):
    sol = 0
    for k in range(max_phi):
        if k != 0:
            sol = sol + k*prob_vec[k]*(z**(k-1))
    sol = sol/mean_deg
    return sol

def get_S_synth_site(degree_list, n):
    size_biggest_component = [0]*n
    p = np.linspace(0,1,n)
    k = 0
    k_s = 0
    for j in range(n):
        graph = config_model(degree_list)
        for i in range(len(p)):
            graph1 = remove_nodes_uniform(graph,p[i])
            size_biggest_component[i] = size_biggest_component[i] + get_size_biggest_component(graph1)
        k +=get_k(graph1)
        k_s += get_k_s(graph1)
        print(j)
    size_biggest_component = [elem/(n*len(degree_list)) for elem in size_biggest_component]
    k = k/n
    k_s = k_s/n
    return size_biggest_component, k, k_s

def get_S_theory_site(prob_vec):
    p = np.linspace(0,1,50)
    degree_list = discrete_samples(prob_vec,1000)
    k1 = sum([degree_list[i] for i in range(len(degree_list))])/len(degree_list)
    k_s1 = sum([degree_list[i]**(2) for i in range(len(degree_list))])/len(degree_list)
    S = []
    for i in range(len(p)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-p[i]+p[i]*g_1(x,prob_vec,np.mean(degree_list)) - x,0.05)
        S.append(p[i]*(1-g_0(x_intersect,prob_vec)))
    return S, k1, k_s1


def get_S_synth_bond(degree_list, n):
    size_biggest_component = [0]*n
    p = np.linspace(0,1,n)
    k = 0
    k_s = 0
    graph = config_model(degree_list)
    for j in range(n):
        for i in range(len(p)):
            graph1 = remove_edges_uniform(graph,p[i])
            size_biggest_component[i] = size_biggest_component[i] + get_size_biggest_component(graph1)
        k +=get_k(graph1)
        k_s += get_k_s(graph1)
        print(j)
    size_biggest_component = [elem/(n*len(degree_list)) for elem in size_biggest_component]
    k = k/n
    k_s = k_s/n
    return size_biggest_component, k, k_s

def get_S_theory_bond(prob_vec):
    p = np.linspace(0,1,100)
    degree_list = discrete_samples(prob_vec,len(prob_vec))
    k1 = sum([degree_list[i] for i in range(len(degree_list))])/len(degree_list)
    k_s1 = sum([degree_list[i]**(2) for i in range(len(degree_list))])/len(degree_list)
    S = []
    for i in range(len(p)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-p[i]+p[i]*g_1(x,prob_vec,np.mean(degree_list)) - x,0.05)
        S.append(abs((1-g_0(x_intersect,prob_vec))))
    return S, k1, k_s1


def get_S_theory_non_uniform(degree_list,prob_vec,mean_deg):
    S = []
    for max_phi in range(max(degree_list)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-f_1(1,prob_vec,max_phi,mean_deg)+f_1(x,prob_vec,max_phi,mean_deg) - x,0.05)
        S.append(f_0(1,prob_vec,max_phi)-f_0(x_intersect,prob_vec,max_phi))
    return S

def get_S_synth_non_uniform(graph):
    components = []
    for max_phi in range(max([deg[1] for deg in graph.degree])):
        graph1 = graph.copy()
        for node in graph.nodes:
            if graph.degree[node] >= max_phi:
                graph1.remove_node(node)
        components.append((get_size_biggest_component(graph1))/len(graph.nodes))
    return components



####### SIR model

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



##### SATRW

def g(x,prob_matrix,t):
    sol = 0
    for i in range(0,len(prob_matrix) -t -1):
        if i >= len(prob_matrix[t]):
            sol+=0
        else:
            sol += (prob_matrix[t][i]*(x**i))
    return sol

def h(x,excess_matrix,t):
    sol = 0
    for i in range(0,len(excess_matrix) -t -2):
        if i >= len(excess_matrix[t]):
            sol+=0
        else:
            sol += (excess_matrix[t][i]*(x**i))
    return sol


def prob_matrix(prob_vec,alpha,tel_rule):
    
    n = len(prob_vec)  # Number of elements in p_0
    # Initialize arrays
    p = np.zeros((n+1, n))
    d = np.zeros((n+1, n))
    q = np.zeros((n+1, n))
    k_avg = np.zeros(n+1)
    r_avg = np.zeros(n+1)
    # Set initial values
    p[0] = prob_vec
    q[-1] = p[0]
    q[0] = p[0]
    if tel_rule == 1:
        d[1] = p[0]
    if tel_rule == 2:
        vec = [alpha*p[0,k-1]+(1-alpha)*p[0,k] for k in range(1,n)]
        #add first element of d[1] equal to 0
        vec.insert(0,0)
        d[1] = vec

    r_avg[1] = np.sum(np.arange(n)*p[0])
    # Compute p, q, k_avg, r_avg at each time step
    for t in range(1, n):
        for k in range(n):
            p[t, k] = (1 / (n-t)) * (p[t-1, k] * (n-t+1) - d[t, k] +
                                        r_avg[t] * (q[t-1, k] - (q[t-1, k-1] if k > 0 else 0)))
        k_avg[t] = np.sum(np.arange(n) * p[t])
        for k in range(n):
            if tel_rule == 1:
                tel_factor = 1
            if tel_rule == 2:
                tel_factor = (k)/k_avg[t]
            q[t, k] = ((k+1)*p[t,k+1]) / k_avg[t] if k < n-1 else 0
            d[t+1, k] = d[t, 0] * tel_factor*p[t, k] + (1 - d[t, 0]) * (alpha *tel_factor*p[t, k] + (1-alpha) * q[t-1, k])
        r_avg[t+1] = np.sum(np.arange(n) * d[t+1])
    return p,q,d

def prob_list(degree_list,n):
    prob_vec = []
    for i in range(n):
        prob_vec.append(degree_list.count(i)/len(degree_list))
    return prob_vec

def SATRW_synth(graph,degree_list,tel_rule,alpha,steps):
    t = steps
    S = []
    P = np.zeros((t,max(degree_list)))
    #get a random node
    graph1 = graph.copy() 
    node= rand.choice(list(graph1.nodes))
    time = 0
    while time <t-1:
        #if the degree of the node is 0, select another node
        if graph1.degree(node) == 0:
            graph1.remove_node(node)
            if tel_rule == 1:
                    node1 = rand.choice(list(graph1.nodes))
            if tel_rule == 2:
                if len(list(graph1.edges)) == 0:
                    node1 = rand.choice(list(graph1.nodes))
                else:
                    edge = rand.choice(list(graph1.edges))
                    node1 = rand.choice(edge)
            #remove the node from the graph
            #else, with probability alpha select a neighbor of the node, with probability 1-alpha select a random node
        else:
            if rand.random() < alpha:
                node1 = rand.choice(list(graph1.neighbors(node)))
                graph1.remove_node(node)
            else:
                graph1.remove_node(node)
                if tel_rule == 1:
                    node1 = rand.choice(list(graph1.nodes))
                if tel_rule == 2:
                    if len(list(graph1.edges)) == 0:
                        node1 = rand.choice(list(graph1.nodes))
                    else:
                        edge = rand.choice(list(graph1.edges))
                        node1 = rand.choice(edge)
        #compute the size of the giant component
        giant = max(nx.connected_components(graph1), key=len)
        S.append(len(giant)/t)
        #compute the degree probability distribution of the new graph
        degree_list = [graph1.degree(i) for i in graph1.nodes]
        prob_vec1 = prob_list(degree_list,P.shape[1])
        P[time] = prob_vec1
        #update the node
        node = node1
        time+=1
    return S,P,graph1

def SATRW_theory(p,q):
    S = []
    for t in range(0,len(p)-3):
        u = sp.optimize.fsolve(lambda x: h(x,q,t) -x, 0.5)[0]
        S1 = (1-g(u,p,t))*(len(p)-t)
        S.append(S1/len(p))
    return S
