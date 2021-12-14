# Implementation of the Maximum Weight Minimization method for Importance Sampling & Rare Event Simulation

@TODO: add the Dynamic fault tree analysis examples.

----------
Current examples:
- Dodecahedron, single and x3 in series (s-t unreliability, undirected graph)
- Campus network with k=5 and k=9 user computers (loss tail probability, directed graph)
- Lattice 4x4 and 6x6 (s-t maximum flow, directed graph)

----------
Require: Python 2.7, Python 3.6

Optional: nauty (for equivalence group partitioning)

----------
Simulation pipeline (dodecahedron and campus network examples):
- Step 1: From mwm.py, run the APPX_RARE_EVENT_SET() function to generate an approximation to the rare event set. The data is stored in a csv file where each line is a rare event sample.
- Step 2: Run optimize.py to obtain the parameters of the Importance Sampling (IS) distribution. This step uses the csv data created in the previous step (NOTE: require Python 3.6 to import scipy.optimize)
- Step 3: Go back to mwm.py, run IMPORTANCE_SAMPLING() using the IS distribution computed in the previous step to perform IS-based rare event simulation.

The lattice examples use exponential tilting for proposal distributions and are implemented in separated Python scripts mwm_lattice.py and optimize_lattice.py.
