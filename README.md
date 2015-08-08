# PathMAN: Path Matrix Algorithm for Networks

This set of scripts calculates statistical properties of the path ensemble for continuous-time random walks (CTRWs) on networks. 

# General Overview
* pathman.py                       - Main algorithm for computing path ensemble properties
* Generate_random_barrier_model.py - Creates the input files of pathman.py for a mutlidimensional RBM

* Generate_lattice.py              - Creates the input files of pathman.py for a general lattice

* RBM_plot_states_prop.py          - Uses the output of pathman.py for a 2-dimensional RBM and displays various properties

Examples of each of the input and output files for these four scripts are found above with the title "2D_ran_barr_lattice" 
# Input for pathmen.py:

* boundary conditions file (file extension .bc), file contents:
  
The boundary conditions file provides a list of the initial states of the system, the probability to start in each state, and a list of final states. The files are to be written in following format.

    line 1 -> Names of initial state and probabilties for starting in each state. (e.g. state_name_1,probability_1  state_name_2,probability_2  etc.)
    line 2 -> Names of final states. (e.g. state_name_1 state_name_2  etc.)

* network file (file extension  .network), file contents:

The network file consists of a line for each of the states of the system. 
The name of the state is followed by a list of neighbors to that state and 
the transition rate to each of those neighbors. The third column consists of the
waiting time moments for transitions out of the current state. The list begins and
requires at least the first moment. The fourth and final column is populated by any 
general state-depenent quantity. An example entry is as follows:

    state_1 neighbor_1,rate_to_neighbor_1;neighbor_2,rate_to_neighbor_2;... time_moment_1,time_moment_2,... state_function_1,state_function_2,...
    
# Output of pathmen.py:
    
* moments file (file extension .moments), file contents:
 
The moments file provides a list of the final states and all of the desired 
moments of the FPT distribution. This file may also contain the path-action distribution 
moments as well if specified. The first line of the file always provides a decriptor for each
column. 

* spatial file (file extension .spatial), file contents:
 
Each line in the spatial file denotes a particular state. The first column specified the name of the 
state. The second and third columns provide the average visits and the average fraction of time spent in 
each state, respectively. 

* length distribution file (file extension .lengths), file contents:

The length distribution file provides the user with each portion of the total time moments 
absorbed at each step $\ell$

