# Implementation of the Maximum Weight Minimization method for Rare Event Simulation


## Examples:
- s-t unreliability: dodecahedron network, lattice 20x20 network (undirected graph)
- Loss tail probability: campus network with k=5 and k=9 user computers (directed graph)
- Maximum flow: lattice 4x4 and 6x6 networks (directed graph)
- Dynamic faule tree: hypothetical computer system example x4 (fault tree)


## Requirement:
- Python3 with networkx, numpy, scipy, matplotlib 
- Optional: nauty (for generating the equivalence group partitioning)


## Simulation pipeline:

For easy of access, each example is stored in a separate folder. For example, to run the **maximum flow** simulation, go to the **maximum-flow** folder and follow these steps:
- Step 1: From the **main.py** file, import the available model (i.e., either **lattice4x4** or **lattice6x6**) and run the **APPX_RE_SET()** function to generate the approximation of the rare event set. The data is stored in a .csv file where each line is a vector that corresponds to a rare event sample. For convenience, the csv files are also provided (i.e., **Y_lattice4x4.csv** and **Y_lattice6x6.csv**).
- Step 2: Run **optimize.py** to obtain the parameters of the proposal distribution. This step uses the .csv data created in the previous step. For convenience, the results have been saved in the python files that store the corresponding model information (i.e., either **lattice4x4.py** or **lattice6x6.py**).
- Step 3: Go back to **main.py**, run the **IMPORTANCE_SAMPLING()** function using the proposal distribution computed in the previous step to estimate the probability of the rare event.
