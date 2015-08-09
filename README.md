# PathMAN: Path Matrix Algorithm for Networks

Path Matrix Algorithm for Networks (c) 2015 Michael Manhart and Willow Kion-Crosby

This set of scripts calculates statistical properties of the path ensemble for continuous-time random walks (CTRWs) on networks. 

# General Overview
* pathman.py                       - Main algorithm for computing path ensemble properties
* Generate_random_barrier_model.py - Creates the input files of pathman.py for a mutlidimensional RBM

* Generate_lattice.py              - Creates the input files of pathman.py for a general lattice

* RBM_plot_states_prop.py          - Uses the output of pathman.py for a 2-dimensional RBM and displays various properties

Examples of each of the input and output files for these four scripts are found above with the title "2D_ran_barr_lattice" 
# Input for pathmen.py:

Note: A description of all commandline options are provided in the help menu:
python pathman.py --help

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
    
* moments file (file extension .moments),
 
The moments file provides a list of the final states and all of the desired 
moments of the FPT distribution. This file may also contain the path-action distribution 
moments as well if specified. The first line of the file always provides a decriptor for each
column. 

* spatial file (file extension .spatial),
 
Each line in the spatial file denotes a particular state. The first column specified the name of the 
state. The second and third columns provide the average visits and the average fraction of time spent in 
each state, respectively. 

* length distribution file (file extension .lengths),

The length distribution file provides the user with each portion of the total time moments 
absorbed at each step. If the path-action is also being calculated, the portion of these
moments at each step are also provided. The first column in the file specifies the current
path step. If any general state-dependent quantity is provided and specified, the average
of that quantity over all states at each given step is also found.

# Input for Generate_random_barrier_model.py:

Note: A description of all commandline options are provided in the help menu:
python Generate_random_barrier_model.py --help

* energy file (file extension .energy),

The energy file of the random barrier model generator must take the form
of a list of each possible transition between states and the corresponding
symmetric energy. The reverse of a given transition does need to be provided
as symmetry rates are assumed in this model. If the file does not exist, this script
will randomly generate the file before producing the final output. An example of 
several lines are as follows:

    line 1 -> state_1->state_2 energy_1
    line 2 -> state_2->state_1 energy_1
    line 3 -> state_2->state_3 energy_2
    ...

# Output of Generate_random_barrier_model.py:

This script generates both necessary files to run a full analysis of the RBM
with the pathman.py script. Examples of both of these files can be found above
under the names "2D_ran_barr_lattice.bc" and "2D_ran_barr_lattice.network". The
.energy file for the system can also be found in the above list.

# Input for Generate_lattice.py:

Note: A description of all commandline options are provided in the help menu:
python Generate_lattice.py --help

* energy file (file extension .energy),

The energy file for the general lattice model generator must take the form
of a list of the energy at each state. If the file does not exist, this script
will randomly generate the file before producing the final output. An example of 
several lines are as follows:

    line 1 -> state_1 energy_1
    line 2 -> state_2 energy_2
    line 3 -> state_3 energy_3
    ...

# Output of Generate_lattice.py:

This script generates both necessary files to run a full analysis for random 
walks on a multidimensional lattice with a general energy landscape
with the pathman.py script.

# Input for RBM_plot_states_prop.py:

On the commandline the name of the file which represents the current RBM to 
be visualized must be included:

python RBM_plot_states_prop.py file_name

The script will then attempt to access files of this name with file extensions .energy, 
.spatial, .lengths, and .network.
