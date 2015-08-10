# PathMAN: Path Matrix Algorithm for Networks

(c) 2015 Michael Manhart and Willow Kion-Crosby

This set of scripts calculates statistical properties of the first-passage path ensemble for continuous-time random walks (CTRWs) on networks.  PathMAN requires Python 2 or 3, NumPy, and SciPy.  Full details of the algorithm and some simple examples are in our paper:

Manhart M, Kion-Crosby W, Morozov AV.  (2015)  "Path statistics, memory, and coarse-graining of continuous-time random walks on networks."  arXiv:[1508.01578](http://arxiv.org/abs/1508.01578).

PathMAN is licensed under GNU General Public License version 3 (see LICENSE for details).  If you use PathMAN in your research, please cite the paper listed above.

# General Overview

Here is an overview of the files included in the repository.  

* pathman.py                       - Main script for calculating path statistics in any CTRW

* Generate_lattice.py              - Generates input files for a simple N-dimensional square lattice lattice to be analyzed by pathman.py

* Generate_random_barrier_model.py - Generates input files for an N-dimensional random barrier model (RBM) on a square lattice to be analyzed by pathman.py

* RBM_plot_states_prop.py          - Makes a few simple plots of path statistics for the 2-dimensional RBM

* RBM_example                      - Contains example input and output files for a 2-dimensional RBM

# Usage for pathman.py:

This is the main script that calculates path statistics for a CTRW on any discrete network of states.  To bring up a list of all command line options for pathman.py, type

    python pathman.py --help

This script requires two input files.  The first defines the network of states and the CTRW.  In the example, this has the extension ".network".  It is easiest to explain by example.  Consider a 1D lattice of length 5, where rates of jumping between neighboring states are all equal.  The network file would be

    # State name    Outgoing jump weights   Waiting time moments    State functions
    1               2,1.0                   1.0,2.0,6.0,24.0        1
    2               1,1.0;3,1.0             0.5,0.5,0.75,1.5        2
    3               2,1.0;4,1.0             0.5,0.5,0.75,1.5        3
    4               3,1.0;5,1.0             0.5,0.5,0.75,1.5        4
    5               4,1.0;6,1.0             0.5,0.5,0.75,1.5        5
    6               5,1.0;7,1.0             0.5,0.5,0.75,1.5        6
    7               6,1.0;8,1.0             0.5,0.5,0.75,1.5        7
    8               7,1.0;9,1.0             0.5,0.5,0.75,1.5        8
    9               8,1.0;10,1.0            0.5,0.5,0.75,1.5        9
    10              9,1.0                   1.0,2.0,6.0,24.0        10

Each line corresponds to a single state in the network.  The first entry in the line shows the state's name (represented by a string without whitespace).  The second entry shows the outgoing jump weights, listing each jump according to the destination state and its corresponding weight separated by a comma.  Different jumps are separated by semicolons.  Note that the weights do not have to be normalized to 1 (the script will do this automatically), but they must be nonnegative.  The third entry lists the waiting time moments for that state, separated by columns.  The first waiting time moment is assumed to be the mean (since the zeroth moment is always 1).  Finally, the fourth column, which is optional, lists any state functions to average over path lengths.  You can list as many values of the state function as you want for each state (separated by commas), but you must list the same number for each state.



* boundary conditions file (file extension .bc), file contents:
  
The boundary conditions file provides a list of the initial states of the system, the probability to start in each state, and a list of final states. The files are to be written in following format.

    line 1 -> Names of initial state and probabilties for starting in each state. (e.g. state_name_1,probability_1  state_name_2,probability_2  etc.)
    line 2 -> Names of final states. (e.g. state_name_1 state_name_2  etc.)

* network file (file extension  .network), file contents:


    
# Output of pathman.py:
    
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
