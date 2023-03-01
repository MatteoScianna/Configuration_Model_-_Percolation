# Configuration_Model_and_Percolation
Enables to create a Configuration Model from a costumizable probability distribution (Poisson or Power Law). 
From this Configuration Model, uniform percolation - bond and site - is performed. 

File functions.py in src folder contains all the functions to create a Configuration Model starting from a given probability distribution. Available possibilities are Poisson distribution and Power Law Distribution. 

Once Configuration Model is created, it is possible to perform both site and bond percolation for the given newtwork. 

File site_perc.py contains an example for the creation of a Configuration Model for a Poisson distribution with lambda = 4, together with a plot of both site and bond percolation comparing the resultant curve with the theoretical one. 

<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/site_percol.png" width="400" height="280">
<img src="https://github.com/MatteoScianna/Configuration_Model_-_Percolation/blob/main/img/bond_percol.png" width="400" height="280">

