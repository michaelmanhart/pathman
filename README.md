# pathman
Path Matrix Algorithm for Networks

This set of scripts calculates statistical properties of the path ensemble for continuous-time random walks (CTRWs) on networks. 

# General Overview
* pathman.py                       - Main algorithm for computing path ensemble properties
*  Generate_random_barrier_model.py - Creates the input files of pathman.py for a mutlidimensional RBM

* Generate_lattice.py              - Creates the input files of pathman.py for a general lattice

* RBM_plot_states_prop.py          - Uses the output of pathman.py for a 2-dimensional RBM and displays various properties

Examples of each the input and output files for these four scripts are found above with the title "2D_ran_barr_lattice" 
# Input for pathmen.py:

boundary conditions file (file extension .bc), file contents:
  
    line 1 -> Names of initial states and probabilties for starting in each state. Corresponding names and probabilites are joined with a comma, and each etry is separated by a tab or space.
    line 2 -> Names of final states. Each entry is separated by a tab or space
  

