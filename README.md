# PathMAN: Path Matrix Algorithm for Networks

(c) 2015 Michael Manhart and Willow Kion-Crosby

This set of scripts calculates statistical properties of the first-passage path ensemble for continuous-time random walks (CTRWs) on networks.  PathMAN requires Python 2 or 3, NumPy, and SciPy.  Full details of the algorithm and some simple examples are in the following paper:

[Manhart M, Kion-Crosby W, Morozov AV.  (2015)  "Path statistics, memory, and coarse-graining of continuous-time random walks on networks."  *J Chem Phys*, in press.  arXiv:1508.01578](http://arxiv.org/abs/1508.01578).

PathMAN is licensed under GNU General Public License version 3 (see LICENSE for details).  If you use PathMAN in your research, please cite the paper listed above.

# Overview

Here is an overview of the files included in the repository:

* `pathman.py`: main script for calculating path statistics for a CTRW
* `Generate_lattice.py`: generates input files for a CTRW on an N-dimensional square lattice to be analyzed by `pathman.py`
* `Generate_RBM.py`: generates input files for a CTRW on a random barrier model (RBM, on an N-dimensional square lattice) to be analyzed by `pathman.py`
* `Plot_RBM_output.py`: makes a few simple plots of path statistics for the 2D RBM
* `Example_1Dlattice`: contains input and output files for an example on a 1D lattice (referred to in the documentation below)
* `Example_RBM`: contains input and output files for an example of the 2D RBM (Fig. 5 in the paper)

# Usage for `pathman.py`

This is the main script that calculates path statistics for a CTRW on any discrete network of states.

## Input files

This script requires two input files.  The **network file** defines the network of states and the jump probabilities and waiting time moments of the CTRW.  The file name is specified with the command line argument `--network`, or if the whole run is named using `--name`, the network file will be assumed to be `NAME.network`.  As an example, consider an ordinary Markov process on a 1D lattice of length 10, where rates of jumping between neighboring states are all equal.  (All input and output files for this example are included in `Example_1Dlattice`.)  The network file would be:

```
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
```

Each line corresponds to a state in the network, with four columns (separated by whitespace) containing properties of that state:

1. The first column shows the state's name (a string without whitespace).  For the 1D lattice example, we just use the integer position on the lattice.
2. The second column shows the outgoing jump weights, listing each jump as a pair of parameters (name of destination state and weight) separated by a comma.  Different pairs of jump parameters are separated by semicolons.  The jump weights, when normalized over all outgoing jumps for a state, will be the jump probabilities; the script will normalize them automatically.  Hence, for an ordinary Markov process the jump weights can be given as the transition rates.  In the example each site has jumps to its nearest neighbors on the lattice with equal rates (arbitrarily set to 1.0).
3. The third column lists the waiting time moments for that state, separated by commas.  Since the zeroth moment is always 1, the first waiting time moment in the list is assumed to be the mean.  Each state must have at least the first moment and at least as many moments as specified by `--max-moment`, which determines the maximum moment of path statistics `pathman.py` will calculate.  In the example we include the first four moments, assuming an ordinary Markov process where the mean waiting time is the inverse of the sum of outgoing jump rates, and the waiting time distribution is exponential.
4. The fourth column, which is optional, lists values of any state functions (separated by commas) to average over, i.e., the script will calculate the average value for each state function at each jump along the path.  The file must include the same number of state function values for all states.  For the example we consider a single state function that is just the position on the lattice.

The second input is the **boundary conditions file**, which defines the path boundary conditions: the initial probability distribution over states and the set of final states.  For example, with the 1D lattice we might have

```
# Initial states and weights
5,0.5 6,0.5

# Final states
1 10
```

This file always has two lines:

1. The first line lists state names and their initial probabilities separated by commas, with spaces separating the different states.  The script assumes any states in the network unlisted here have zero initial probability.  The initial probabilities do not need to be normalized, as the script will do this automatically.  For the 1D lattice example, the initial distribution is localized in the middle of the lattice, half at position 5 and half at position 6.
2. The second line lists the final states separated by spaces.  For the example, there are two final states, one at each edge of the lattice (positions 1 and 10).

## Output

The script will print basic output to stdout.  For the 1D lattice example, the output should look like this:

```
Command:                          python pathman.py --name Example_1Dlattice/1Dlattice --max-moment 4
Reading network from:             Example_1Dlattice/1Dlattice.network
Reading boundary conditions from: Example_1Dlattice/1Dlattice.bc
Writing output to:                Example_1Dlattice/1Dlattice

Initializing...
	Reading input files...done.                                     (0.0 seconds)
	Processing jump probabilities and waiting time moments...done.  (0.002 seconds)
	Initializing linear algebra...done.                             (0.005 seconds)
	Total states:                     10
	Final states:                     1
	Average (outward) connectivity:   1.7777777777777777
	Initial unnormalized probability: 1.0
Running...
	100 jumps...
	Done: converged after 360 jumps.                                (0.01 seconds)
Writing data...done.                                                (0.003 seconds)
    
Total path moments:
	Length             Raw                 Cumulant            Standardized        
	lbar0              0.999999990767      0.999999990767      0.999999990767      
	lbar1              20.499996492        20.499996492        1.09231990307e-08   
	lbar2              720.498663503       300.248807329       1.0                 
	lbar3              38635.9892758       11555.5702076       2.22110687057       
	lbar4              2849660.68724       697969.003263       10.742363703        

	Time               Raw                 Cumulant            Standardized        
	tbar0              0.999999990767      0.999999990767      0.999999990767      
	tbar1              12.4999978952       12.4999978952       1.05908454071e-08   
	tbar2              274.999517268       118.749569888       1.0                 
	tbar3              9179.88858872       2773.65645437       2.14340683405       
	tbar4              419014.118459       102285.207712       10.2535190307       

Run started at:		2015-11-06 21:08:34.439408
Run ended at:		2015-11-06 21:08:34.461347
Total run time:		0.022 seconds
Peak memory usage:	0.17212 GB
```

Besides printing some basic information about the run (command used, input and output file names), `pathman.py` prints a few basic statistics on the network (number of states, connectivity).  After the numerical calculation has converged, `pathman.py` prints the total moments (averaged over final states) for both path length and time (as well as path action if requested; see options).

A few output files contain more detailed path statistics.  These are all labeled using the argument `--name`, but with different extensions:

* Moments file (`NAME.moments`): this contains a list of final states with path time moments for each.  For the 1D example:
```
# Final_state       tbar0               tbar1               tbar2               tbar3               tbar4   
5                   0.999999990767      12.4999978952       274.999517268       9179.88858872       419014.118459 
```
* Spatial distribution file (`NAME.spatial`): this contains a list of all states with the average number of visits and the average fraction of time at each.  For final states, the average number of visits corresponds to the commitment probability.  For the example:
```
# State             Average_visits      Average_fraction_of_time
1                   2.0                 0.160000026941
2                   3.0                 0.120000020206
3                   2.0                 0.0800000134706
4                   0.999999999999      0.0400000067353
5                   0.999999990767      0.0
6                   0.999999983298      0.0400000060672
7                   1.9999999666        0.0800000121345
8                   2.99999995627       0.120000018457
9                   3.99999994595       0.160000024779
10                  2.49999997298       0.200000031515
```
* Length distribution file (`NAME.lengths`): this gives the total time moments absorbed at final states at each jump along the path.  For the example, 
```
# l                 tbar0               tbar1               tbar2               tbar3               tbar4           
0                   0.0                 0.0                 0.0                 0.0                 0.0         
1                   0.0                 0.0                 0.0                 0.0                 0.0             
2                   0.0                 0.0                 0.0                 0.0                 0.0             
3                   0.0                 0.0                 0.0                 0.0                 0.0             
4                   0.0625              0.15625             0.5                 1.96875             9.28125         
5                   0.03125             0.09375             0.34375             1.5                 7.640625        
6                   0.0625              0.234375            1.046875            5.4609375           32.765625       
7                   0.0390625           0.1640625           0.80078125          4.4765625           28.3359375      
8                   0.0546875           0.2734375           1.56640625          10.16015625         73.8984375      
9                   0.0390625           0.2109375           1.28515625          8.75390625          66.1640625      
10                  0.046875            0.29296875          2.046875            15.861328125        135.41015625
...
360                 4.33200405896e-14   9.74700913266e-12   2.20057993487e-09   4.9852056295e-07    0.000113319869237 
```

## Options

To bring up a list of all command line options and their descriptions for PathMAN, type

```
python pathman.py --help
```


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

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

