# Implementation of the Maximum Weight Minimization method for Importance Sampling & Rare Event Simulation

@TODO: add the Dynamic fault tree analysis examples.

----------
Examples:
- s-t unreliability: dodecahedron network, lattice 20x20 network (undirected graph)
- Loss tail probability: campus network with k=5 and k=9 user computers (directed graph)
- Maximum flow: lattice 4x4 and 6x6 networks (directed graph)
- Dynamic faule tree: hypothetical computer system example x4 (fault tree)
----------
Require: Python3
Optional: nauty (for equivalence group partitioning)

----------
Simulation pipeline:

Each example is stored in a separate folder. To run the simulation:
- Step 1: From main.py, import the model and run the APPX_RARE_EVENT_SET() function to generate the approximation of the rare event set. The data is stored in a csv file where each line is a vector that corresponds to a boundary rare event sample.
- Step 2: Run optimize.py to obtain the parameters of the proposal distribution. This step uses the csv data created in the previous step. For convenient, the results have been saved in the python file of the corresponding model.
- Step 3: Go back to main.py, run IMPORTANCE_SAMPLING() using the proposal distribution computed in the previous step to estimate the probability of the rare event.
