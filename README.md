# PathMAN: Path Matrix Algorithm for Networks

(c) 2015 Michael Manhart and Willow Kion-Crosby

This set of scripts calculates statistical properties of the first-passage path ensemble for continuous-time random walks (CTRWs) on networks.  PathMAN requires Python 2 or 3, NumPy, and SciPy.  Full details of the algorithm and some simple examples are in the following paper:

[Manhart M, Kion-Crosby W, Morozov AV.  (2015)  "Path statistics, memory, and coarse-graining of continuous-time random walks on networks."  *J Chem Phys* **143**:214106.](https://dx.doi.org/10.1063/1.4935968).

PathMAN is licensed under GNU General Public License version 3 (see LICENSE for details).  If you use PathMAN in your research, please cite the paper listed above.

# Overview

Here is an overview of the files included in the repository:

* `pathman.py`: main script that calculates path statistics for a CTRW
* `Generate_lattice.py`: generates input files for a CTRW on an multidimensional square lattice to be analyzed by `pathman.py`
* `Generate_RBM.py`: generates input files for a CTRW on a random barrier model (RBM, on an multidimensional square lattice) to be analyzed by `pathman.py`
* `Plot_RBM_output.py`: makes a few simple plots of path statistics for the 2D RBM (similar to Fig. 5 in the aforementioned paper)
* `Example_1Dlattice`: contains input and output files for an example on a 1D lattice (referred to in the documentation below)
* `Example_RBM`: contains input and output files for an example of the 2D RBM (same as Fig. 5 in the paper)

# Usage for `pathman.py`

This is the main script that calculates path statistics for a CTRW on any discrete network of states.

## Input files

This script requires two input files.  The **network file** defines the network of states and their jump probabilities and waiting time moments in the CTRW.  The file name is specified with the command line argument `--network`, or if the whole run is named using `--name`, the network file will be assumed to be `NAME.network`.  As an example, consider an ordinary Markov process on a 1D lattice of length 10, where rates of jumping between neighboring states are all equal.  (All input and output files for this example are included in `Example_1Dlattice`.)  The network file would be:

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
4. The fourth column, which is optional, lists values of any state functions (separated by commas) to average over, i.e., the script will calculate the average value of each state function at each jump along the path.  The file must include the same number of state function values for all states.  For the example we consider a single state function that is just the position on the lattice.

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

Empty lines or lines beginning with `#` are ignored in all input files.

## Output

The script will print basic output to stdout.  For the 1D lattice example, the output should look like this:

```
Command:                          python pathman.py --name Example_1Dlattice/1Dlattice --max-moment 4 --num-state-funcs 1
Reading network from:             Example_1Dlattice/1Dlattice.network
Reading boundary conditions from: Example_1Dlattice/1Dlattice.bc
Writing output to:                Example_1Dlattice/1Dlattice

Initializing...
	Reading input files...done.                                     (0.0 seconds)
	Processing jump probabilities and waiting time moments...done.  (0.003 seconds)
	Initializing linear algebra...done.                             (0.006 seconds)
	Total states:                     10
	Final states:                     2
	Average (outward) connectivity:   2.0
	Initial unnormalized probability: 1.0
Running...
	100 jumps...
	Done: converged after 407 jumps.                                (0.016 seconds)
Writing data...done.                                            	(0.004 seconds)

Total path moments:
	Length             Raw                 Cumulant            Standardized        
	lbar0              0.999999999987      0.999999999987      0.999999999987      
	lbar1              19.9999999947       19.9999999947       1.55815137326e-11   
	lbar2              659.999997743       259.999997956       1.0                 
	lbar3              31915.9990411       8315.99917428       1.9836014666        
	lbar4              2052899.59194       400819.666768       8.92928510491       

	Time               Raw                 Cumulant            Standardized        
	tbar0              0.999999999987      0.999999999987      0.999999999987      
	tbar1              9.99999999734       9.99999999734       1.50147223248e-11   
	tbar2              169.999999434       69.9999994876       1.0                 
	tbar3              4241.99987929       1141.99989602       1.94992994659       
	tbar4              140735.974135       28355.9788619       8.78693454663       

Run started at:		2015-11-07 00:56:06.612505
Run ended at:		2015-11-07 00:56:06.643924
Total run time:		0.031 seconds
Peak memory usage:	0.172104 GB
```

Besides printing some basic information about the run (command used, input and output file names), `pathman.py` prints a few basic statistics on the network (number of states, connectivity).  After the numerical calculation has converged, `pathman.py` prints the total moments (averaged over final states) for both path length and time (as well as path action if requested; see options).

A few output files contain more detailed path statistics.  These are all labeled using the argument `--name`, but with different extensions:

* Moments file (`NAME.moments`): this contains a list of final states with path time moments (and path action moments if requested) for each.  The moments are unconditional, i.e., they are weighted by the total probability of reaching that final state from the initial distribution (commitment probability).  For the 1D example, this will be:
```
# Final_state       tbar0               tbar1               tbar2               tbar3               tbar4         
1                   0.499999999994      4.99999999867       84.9999997172       2120.99993964       70367.9870674 
10                  0.499999999994      4.99999999867       84.9999997172       2120.99993964       70367.9870674
```
* Spatial distribution file (`NAME.spatial`): this contains a list of all states with the average number of visits and the average fraction of time at each.  For final states, the average number of visits corresponds to the commitment probability.  For the example this will look like:
```
# State             Average_visits      Average_fraction_of_time
1                   0.499999999994      0.0
2                   0.999999999988      0.0500000000127
3                   1.99999999998       0.100000000025
4                   2.99999999997       0.150000000038
5                   3.99999999997       0.200000000052
6                   3.99999999997       0.200000000052
7                   2.99999999997       0.150000000038
8                   1.99999999998       0.100000000025
9                   0.999999999988      0.0500000000127
10                  0.499999999994      0.0
```
* Length distribution file (`NAME.lengths`): this gives the total time moments absorbed at final states at each jump (`l`) along the path, along with the average value of any state functions at that jump.  For the example this will be:
```
# l                 tbar0               tbar1               tbar2               tbar3               tbar4               state_func0         
0                   0.0                 0.0                 0.0                 0.0                 0.0                 5.5                 
1                   0.0                 0.0                 0.0                 0.0                 0.0                 5.5                 
2                   0.0                 0.0                 0.0                 0.0                 0.0                 5.5                 
3                   0.0                 0.0                 0.0                 0.0                 0.0                 5.5                 
4                   0.0625              0.125               0.3125              0.9375              3.28125             5.5                 
5                   0.03125             0.078125            0.234375            0.8203125           3.28125             5.5                 
6                   0.0625              0.1875              0.65625             2.625               11.8125             5.5                 
7                   0.0390625           0.13671875          0.546875            2.4609375           12.3046875          5.5                 
8                   0.0546875           0.21875             0.984375            4.921875            27.0703125          5.5                 
9                   0.0390625           0.17578125          0.87890625          4.833984375         29.00390625         5.5                 
10                  0.046875            0.234375            1.2890625           7.734375            50.2734375          5.5 
...
407                 8.0618015937e-13    1.64057662432e-10   3.34677631361e-08   6.84415756133e-06   0.00140305230007    5.5    
```

## Options

To bring up a list of all command line options and their descriptions for PathMAN, type

```
python pathman.py --help
```

Besides options for specifying the input and output files, there are options to tell `pathman.py` to calculate statistics of path action and state functions and to modify the convergence conditions.

# Usage for `Generate_lattice.py`

This script generates input files (network and boundary condition files) that can be read by `pathman.py` for a CTRW on an N-dimensional square lattice.  The help menu contains descriptions of all command line options:

```
python Generate_lattice.py --help
```

In particular, the user can specify the lattice dimension and number of points, as well as initial and final states for the CTRW.  The CTRW is assumed to be a Markov process with Metropolis transition rates on an energy landscape over the lattice.  The user can specify the energy landscape in an input file using `--energy`.  For example, a linear energy ramp on the 1D lattice would be

```
# State     Energy
1           1.0
2           0.8888888888888889
3           0.7777777777777778
4           0.6666666666666667
5           0.5555555555555556
6           0.4444444444444444
7           0.3333333333333333
8           0.2222222222222222
9           0.1111111111111111
10          0.0
```

If no input energy file is given, the script assumes a linear energy function.

# Usage for `Generate_RBM.py`

This script is similar to `Generate_lattice.py`, but for the random barrier model.  As before the help menu contains descriptions of all command line options:

```
python Generate_RBM.py --help
```

The user can specify the energy barriers in an input file using `--energy`.  Since an energy barrier lies between states, each line specifies a pair of neighboring states and the height of the energy barrier between them.  The energy barriers must be symmetric.  For example, on a 1D lattice we might have

```
# State_pair     Energy
1-2              1.0
2-1              1.0
2-3              1.0
3-2              1.0
3-4              1.0
4-3              1.0
4-5              1.0
5-4              1.0
5-6              1.0
6-5              1.0
6-7              1.0
7-6              1.0
7-8              1.0
8-7              1.0
8-9              1.0
9-8              1.0
9-10             1.0
10-9             1.0
```

# Usage for `Plot_RBM_output.py`:

This script makes a few simple plots of path statistics for the 2D RBM, similar to each row of Fig. 5 in the paper.  The usage for this script is

```
python Plot_RBM_output.py NAME
```

where the RBM energy file and input/output files from `pathman.py` all have root `NAME` (e.g., `NAME.energy`, `NAME.network`, `NAME.spatial`, `NAME.lengths`).
