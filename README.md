# pathman
Path Matrix Algorithm for Networks

This set of scripts calculates statistical properties of the path ensemble for continuous-time random walks (CTRWs) on networks. 

pathman.py - Main algorithm for computing path ensemble properties
Generate_random_barrier_model.py - Creates the input files of pathman.py 

Input for pathmen.py:

  boundary conditions file; file extension .bc
  contents:
  
    line 1 -> names of initial states and probabilties for starting in each state
    line 2 -> 
