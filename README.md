# Vaccination Strategies and Epidemic Scenarios
In this repositories you can find different tools in order to 
- create graphs from a given degree distribution (configuration model);
- apply different percolation methods to a given graph;
- simulate the evolution of a SIR model on a diluited network.

File functions.py in src folder contains all the functions to do so. Available possibilities for the probaiblity distribution of the graph are Poisson and Power Law, in order to obtain ER and Scale Free networks. 

## Percolation

Once Configuration Model is created, it is possible to perform different kind of percolation processes on the network. Available percolation models are site percolation, bond percolation and degree-based percolation. 
File percolation.py contains an example for all these scenarios. Plots below show the evolution of the biggest component of a Configuration Model comparing synthetic and theoretical results. 

<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/site_percolation.png" width="400" height="280" /> <img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/bond_percolation.png" width="400" height="280"/>
<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/non_uniform_percolation.jpg" width="680" height="280">

## SATRW 

Another possible way to diluite a network is using the Self Avoiding Teleport Random Walk strategy. It is a process that mixes together standard random walk and different possible teleporting rules. At each step a new node is selected either within the neighbors of the previous node or by the selected teleporting rule and hence it is removed from the network. 
File SATRW.py contains scripts to perform this technique on a given network. Below, plots show the evolution of the time dependent probability distribution and the evolution of the size of the giant component of the network (comparison between theoretical and simulated results). 


<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/degree_dist_ER.png" width="500" height="400" /> <img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/degree_dist_SF.png" width="500" height="400" />
<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/giant_comp_ER.png" width="500" height="400" /> <img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/giant_comp_SF.png" width="500" height="400" />

## SIR Model

File SIR.py contains the code for the implementation of an SIR model for a given graph. 

## Epidemic Scenarios on Diluited Networks

Once we have defined different possible vaccination strategies and implemented the evolution of a possible diseas spread, it is possible to investigate whether the final size of the epidemic changes according to different vaccination strategies. File VAX+SIR.py considers all the different vaccination strategies mentioned before (uniform site percolation, degree-based site percolation, SATRW with uniform teleportation and SATRW with biased teleportaiton) and run an SIR on the corresponding diluited network. 
Heatmaps below show the impact of an epidemic given a certain value of vaccinated individuals and a different value of the transmissibility parameter of the epidemic. 

<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/heatmaps-ER.png" width="500" height="400" /> <img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/heatmaps-SF.png" width="500" height="400" />


