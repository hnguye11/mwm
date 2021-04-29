# Implementation of the Maximum Weight Minimization method for Importance Sampling & Rare Event Simulation

@TODO: add the Dynamic fault tree analysis examples.

----------
Current examples:
- Dodecahedron, single and x3 in series (s-t unreliability, undirected graph)
- Campus network with k=5 and k=9 user computers (s-rooted spanning arborescence, directed graph)
- Lattice 4x4 and 6x6 (s-t maximum flow, directed graph)

----------
Require: Python 2.7, Python 3.6
Optional: nauty (for equivalence group partitioning)

----------
Simulation pipeline (dodecahedron and campus network):
- Step 1: From mwm.py, run the APPX_RARE_EVENT_SET() function to generate the approximation of the rare event set. The data is stored in a csv file where each line is a rare event sample
- Step 2: Run optimize.py to obtain the parameters of the Importance Sampling (IS) distribution. This step uses the csv data created in the previous step (NOTE: require Python 3.6 to run scipy.optimize)
- Step 3: Go back to mwm.py, run IMPORTANCE_SAMPLING() using the IS distribution computed in the previous step to perform IS-based Monte Carlo simulation.

The lattice examples use exponential tilting for proposal distributions and are implemented in separated Python scripts mwm_lattice.py and optimize_lattice.py.

----------
If you have trouble running the examples, please feel free to reach out to me at hnguye11(at)illinois(dot)edu or nhh311(at)gmail(dot)com. 
