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


def g_1(x,prob_list,mean_deg):
    sol = 0
    for i in range(10000-1):
        q = (prob_list[i+1]*(i+1))/mean_deg
        sol += q*(x**(i))
    return sol

def g_0(x,prob_list):
    sol = 0
    for i in range(10000):
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


def prob_list_poisson(lambda1,min):
    prob_vec = [0]*min
    for k in range(min,50000):
        prob_vec.append((poisson.pmf(k,lambda1)))
    prob_vec = prob_vec/sum(prob_vec)
    return prob_vec

def prob_list_powerlaw(a,min):
    prob_vec = [0]*min
    for k in range(min,10000):
        prob_vec.append(k**(-a))
    prob_vec = [i/sum(prob_vec) for i in prob_vec]
    return prob_vec

def discrete_inverse_trans(prob_vec):
    U=uniform.rvs(size=1)
    if U<=prob_vec[0]:
        return 1
    else:
        for i in range(1,len(prob_vec)+1):
            if sum(prob_vec[0:i])<U and sum(prob_vec[0:i+1])>U:
                return i

def discrete_samples(prob_vec,n=1):
    sample=[]
    for i in range(0,n):
        sample.append(discrete_inverse_trans(prob_vec))
    #check if sum is even   
    if sum(sample)%2 != 0:
        sample[0] = sample[0]+1
    return sample

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

def remove_edges_uniform(graph,p):
    graph1 = graph.copy()
    for edge in graph.edges:
        if rand.random() > p:
            graph1.remove_edge(edge[0],edge[1])
    return graph1


def remove_nodes_uniform(graph,p):
    graph1 = graph.copy()
    for node in graph.nodes:
        if rand.random() > p:
            graph1.remove_node(node)
    return graph1

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
    size_biggest_component = [0]*100
    p = np.linspace(0,1,100)
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
    p = np.linspace(0,1,10)
    degree_list = discrete_samples(prob_vec,10000)
    k1 = sum([degree_list[i] for i in range(len(degree_list))])/len(degree_list)
    k_s1 = sum([degree_list[i]**(2) for i in range(len(degree_list))])/len(degree_list)
    S = []
    for i in range(len(p)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-p[i]+p[i]*g_1(x,prob_vec,np.mean(degree_list)) - x,0.05)
        S.append(p[i]*(1-g_0(x_intersect,prob_vec)))
    return S, k1, k_s1


def get_S_synth_bond(degree_list, n):
    size_biggest_component = [0]*100
    p = np.linspace(0,1,100)
    k = 0
    k_s = 0
    for j in range(n):
        graph = config_model(degree_list)
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
    p = np.linspace(0,1,10)
    degree_list = discrete_samples(prob_vec,10000)
    k1 = sum([degree_list[i] for i in range(len(degree_list))])/len(degree_list)
    k_s1 = sum([degree_list[i]**(2) for i in range(len(degree_list))])/len(degree_list)
    S = []
    for i in range(len(p)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-p[i]+p[i]*g_1(x,prob_vec,np.mean(degree_list)) - x,0.05)
        S.append((1-g_0(x_intersect,prob_vec)))
    return S, k1, k_s1


def get_S_theory_non_uniform(degree_list,prob_vec):
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


