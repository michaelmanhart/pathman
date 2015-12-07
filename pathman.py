################################################################################
#    PathMAN: Path Matrix Algorithm for Networks
#    Copyright (c) 2015 Michael Manhart and Willow Kion-Crosby
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	E-mail: mmanhart@fas.harvard.edu
#	Mail: 
#		Harvard University
#		Department of Chemistry and Chemical Biology
#		12 Oxford Street
#		Cambridge, MA 02138, USA
################################################################################


import sys
import argparse
import datetime
from time import time
import numpy
import scipy.sparse
import scipy.special
from bisect import bisect_left


################################################################################
# This class contains the parameters for the CTRW model and carries out all 
# calculations.  Initializing the class reads the input files and defines all
# matrix and vector objects.  The Run function iterates over the recursion
# relations to perform sums over all paths.  The Output_Data function then 
# writes all data to files.
################################################################################
class Pathman:

	###########################################################################
	# This initializes an instance of the CTRW.
	###########################################################################
	def __init__(self, args):

		# Set run parameters
		self.name = args.name
		self.max_moment = args.max_moment
		self.calc_action = args.calc_action
		self.num_state_funcs = args.num_state_funcs
		self.length_res = args.length_res
		self.length_min = args.length_min
		self.no_length_dists = args.no_length_dists
		self.debug = args.debug
		self.converge = args.converge
		self.max_jumps = args.max_jumps
		self.no_converge = args.no_converge

		# Set input file names
		network_file = args.network
		bc_file = args.bc

		######################################################################
		# Read the network file and store data, create list of all states, and
		# count total number
		######################################################################

		mytime = time()
		Print("\tReading input files...")

		# Store all network data in a list
		try:
			myFile = open(network_file, "r")
		except IOError:
			Print("\nERROR: Cannot read network file \"" + str(network_file) + "\".  Exiting.\n")
			exit()
		network_data = [ [entry for entry in line.split()] for line in myFile.readlines() if not (line.startswith('#') or line.isspace()) ]
		myFile.close()

		# Extract states from first entry of each line, sort them into
		# lexicographical order, and count total number
		self.all_states = sorted([row[0] for row in network_data])
		self.num_states = len(self.all_states)

		# Check for duplicates
		if len(set(self.all_states)) != self.num_states:
			Print("\nERROR: Duplicate states in network file.  Exiting.\n")
			exit()

		######################################################################
		# Initialize boundary conditions: create vector of initial state
		# probabilities and set of final states
		######################################################################

		# Set initial state probability distribution to all zeros
		initial_distribution = numpy.zeros(self.num_states)

		# Read boundary condition file
		try:
			myFile = open(bc_file, "r")
		except IOError:
			Print("\nERROR: Cannot read boundary condition file \"" + str(bc_file) + "\".  Exiting.\n")
			exit()
		bc_data = [ [entry for entry in line.split()] for line in myFile.readlines() if not (line.startswith('#') or line.isspace()) ]
		myFile.close()

		# Extract initial states and their weights from first line of file
		for pair in bc_data[0]:
			pair = pair.split(",")
			s, weight = pair[0], float(pair[1])
			if s not in self.all_states:
				Print("\nERROR: Initial state " + s + " not in network file.  Exiting.\n")
				exit()
			initial_distribution[self.Get_State_Index(s)] = weight

		# Normalize initial distribution
		initial_normalization = initial_distribution.sum()
		initial_distribution /= initial_normalization

		# Extract set of all final states from second line of file
		self.final_states = set(bc_data[1])
		if not (self.final_states <= set(self.all_states)):
			Print("\nERROR: At least one final state is not in the network file.  Exiting.\n")
			exit()

		Print("done.".ljust(42) + "(" + str(round(time() - mytime, 3)) + " seconds)\n")
	
		######################################################################
		# Initialize matrices for jump probabilities and waiting time moments
		######################################################################

		mytime = time()
		Print("\tProcessing jump probabilities and waiting time moments...")

		# Jump matrix (Q), dimensions N x N
		jump_matrix = scipy.sparse.dok_matrix((self.num_states, self.num_states), dtype=float)

		# Jump action moment matrices (Qtilde), dimensions N x N
		jam_matrices = [scipy.sparse.dok_matrix((self.num_states, self.num_states), dtype=float) for i in range(self.max_moment + 1)]

		# Waiting time moment matrices (Theta), list of nmax+1 N x N diagonal 
		# matrices with the first one being an identity
		self.wtm_matrices = [scipy.sparse.identity(self.num_states, format="dok")] + [scipy.sparse.dok_matrix((self.num_states, self.num_states), dtype=float) for i in range(self.max_moment)]

		# Functions of state (B), list of N x N diagonal matrices, one for 
		# each state function in network file
		state_func_matrices = [scipy.sparse.dok_matrix((self.num_states, self.num_states), dtype=float) for i in range(self.num_state_funcs)]

		# Loop over states in network data
		edges = 0
		for row in network_data:

			# Get state name and index
			s = row[0]
			ind = self.Get_State_Index(s)

			# Extract state functions
			if self.num_state_funcs > 0:
				sf_data = row[3].split(",")
				if len(sf_data) < self.num_state_funcs:
					Print("\nERROR: Not enough state function data given for state " + s + ".  Exiting.\n")
					exit()
				for i in range(self.num_state_funcs):
					state_func_matrices[i][ind, ind] = float(sf_data[i])

			# Leave zero jump probabilities and WTMs for final states
			if s in self.final_states: continue

			# Extract jump weights and normalize
			jump_data = [(jump.split(",")[0], float(jump.split(",")[1])) for jump in row[1].split(";")]
			edges += len(jump_data)
			total_weight = sum(entry[1] for entry in jump_data)
			for (neighbor, weight) in jump_data:
				q = weight/total_weight			
				#if q > args.jump_cutoff:
				ind_neighbor = self.Get_State_Index(neighbor)
				jump_matrix[ind_neighbor, ind] = q
				if self.calc_action:
					for i in range(self.max_moment + 1):
						jam_matrices[i][ind_neighbor, ind] = q*(-numpy.log(q))**i

			# Extract waiting time moments
			moment_data = row[2].split(",")
			if self.max_moment > len(moment_data):
				Print("\nERROR: Not enough waiting time moments given for state " + s + ".  Exiting.\n")
				exit()
			for i in range(1, self.max_moment + 1):
				self.wtm_matrices[i][ind, ind] = float(moment_data[i-1])

		# Convert matrices to CSR for linear algebra
		jump_matrix = jump_matrix.tocsr()
		self.wtm_matrices = [m.tocsr() for m in self.wtm_matrices]
		jam_matrices = [m.tocsr() for m in jam_matrices]
		state_func_matrices = [m.tocsr() for m in state_func_matrices]

		# List of waiting time moment matrices weighted by the jump matrix: [Q.Theta0, Q.Theta1, ... ]
		weighted_jump_matrices = [jump_matrix.dot(self.wtm_matrices[i]) for i in range(self.max_moment + 1)]

		Print("done.".ljust(7) + "(" + str(round(time() - mytime, 3)) + " seconds)\n")

		######################################################################
		# Initialize sparse matrices for linear algebra
		######################################################################

		mytime = time()
		Print("\tInitializing linear algebra...")

		# Construct final state "bra" <final|: 1 x N matrix with ones for each final state
		# Construct intermediate state "bra" <int|: 1 x N matrix with ones for each intermediate (non-absorbing) state
		final_bra = scipy.sparse.dok_matrix((1, self.num_states), dtype="float")
		int_bra = scipy.sparse.dok_matrix((1, self.num_states), dtype="float")
		for s in self.all_states: 
			if s in self.final_states:
				final_bra[0, self.Get_State_Index(s)] = 1.
			else:
				int_bra[0, self.Get_State_Index(s)] = 1.

		# Construct F matrix ( (nmax+1) x N(nmax+1) ): nmax+1 copies of <final| along the diagonal
		self.final_matrix = scipy.sparse.block_diag([final_bra for i in range(self.max_moment + 1)], format="csr", dtype="float")

		if self.num_state_funcs > 0:
			# Construct B matrices (number of state functions x N(nmax+1)) for intermediate and final states 
			self.int_B_matrix = scipy.sparse.vstack( [int_bra.dot(state_func_matrices[i]) for i in range(self.num_state_funcs)] )
			self.int_B_matrix = scipy.sparse.hstack( [self.int_B_matrix, scipy.sparse.csr_matrix((self.num_state_funcs, self.num_states*self.max_moment), dtype="float")], format="csr")

			self.final_B_matrix = scipy.sparse.vstack( [final_bra.dot(state_func_matrices[i]) for i in range(self.num_state_funcs)] )
			self.final_B_matrix = scipy.sparse.hstack( [self.final_B_matrix, scipy.sparse.csr_matrix((self.num_state_funcs, self.num_states*self.max_moment), dtype="float")], format="csr")

		# Calculate binomial coefficients
		binomials = numpy.zeros( (self.max_moment + 1, self.max_moment + 1) )
		for n in range(self.max_moment + 1):
			for j in range(n + 1):
				binomials[n][j] = scipy.special.gamma(n + 1)/scipy.special.gamma(j + 1)/scipy.special.gamma(n - j + 1)

		# Construct the time transfer matrix K ( N(nmax+1) x N(nmax+1) )
		blocks = [ [binomials[i][i-j]*weighted_jump_matrices[i - j] for j in range(i + 1)] + [None for j in range(i + 1, self.max_moment + 1)] for i in range(self.max_moment + 1) ]
		self.time_transfer_matrix = scipy.sparse.bmat(blocks, format="csr", dtype="float")

		# Construct the action transfer matrix G ( N(nmax+1) x N(nmax+1) )
		blocks = [ [binomials[i][i-j]*jam_matrices[i - j] for j in range(i + 1)] + [None for j in range(i + 1, self.max_moment + 1)] for i in range(self.max_moment + 1) ]
		self.action_transfer_matrix = scipy.sparse.bmat(blocks, format="csr", dtype="float")

		# Construct initial transfer vectors |tau(0)> and |eta(0)>
		self.time_transfer_ket = numpy.concatenate( [initial_distribution, numpy.zeros(self.max_moment*self.num_states)] )
		self.action_transfer_ket = numpy.concatenate( [initial_distribution, numpy.zeros(self.max_moment*self.num_states)] )
		
		# Initialize vectors of cumulative moments at each state (running sums 
		# of transfer_vectors after each jump) |tau> and |eta>
		self.time_cumulative_ket = self.time_transfer_ket
		self.action_cumulative_ket = self.action_transfer_ket

		self.lbars = numpy.zeros(self.max_moment + 1)

		# Initialize tbar(l), a list of kets where the first index is the 
		# path length l and the second index (in the ket) is the moment.  It 
		# gives the total amount of time moment of all complete paths of 
		# length l
		self.tbarl = [] #[self.final_matrix.dot(self.time_transfer_ket)]

		if self.num_state_funcs > 0:
			# Initialize Bbar(l); it is the average of state function h at the
			# lth jump
			self.Bbarl = [] #[self.state_func_matrix.dot(self.time_transfer_ket)]

		Print("done.".ljust(34) + "(" + str(round(time() - mytime, 3)) + " seconds)\n")

		Print("\tTotal states:".ljust(35) + str(jump_matrix.shape[0]) + "\n")
		Print("\tFinal states:".ljust(35) + str(len(self.final_states)) + "\n")
		Print("\tAverage (outward) connectivity:".ljust(35) + str(edges/float(self.num_states - len(self.final_states))) + "\n")
		Print("\tInitial unnormalized probability:".ljust(35) + str(initial_normalization) + "\n")


	###########################################################################
	# This returns the index of state s in the sorted list of all states
	###########################################################################
	def Get_State_Index(self, s):
		ind = bisect_left(self.all_states, s)
		if ind != len(self.all_states) and self.all_states[ind] == s:
			return ind
		raise ValueError


	###########################################################################
	# This returns a ket with the total cumulative path length moments over all 
	# final states.
	###########################################################################
	def Total_Length_Moments(self):
		return self.lbars


	###########################################################################
	# This returns a ket with the total cumulative path time moments over all 
	# final states.
	###########################################################################
	def Total_Time_Moments(self):
		return self.final_matrix.dot(self.time_cumulative_ket)


	###########################################################################
	# This returns a ket with the total cumulative path action moments over all 
	# final states.
	###########################################################################
	def Total_Action_Moments(self):
		return self.final_matrix.dot(self.action_cumulative_ket)

		
	###########################################################################
	# This updates the transfer matrices for the lth jump.  It also absorbs 
	# moments from the final states and updates the running totals.
	###########################################################################
	def Update(self, l):

		# Update the time transfer ket and increment the cumulative sum
		self.time_transfer_ket = self.time_transfer_matrix.dot(self.time_transfer_ket)
		self.time_cumulative_ket += self.time_transfer_ket

		# Update the length moments
		rho = self.final_matrix.dot(self.time_transfer_ket)[0]
		self.lbars += rho*numpy.array([float(l)**n for n in range(self.max_moment + 1)])

		# Update and increment action kets if indicated
		if self.calc_action:
			self.action_transfer_ket = self.action_transfer_matrix.dot(self.action_transfer_ket)
			self.action_cumulative_ket += self.action_transfer_ket


	###########################################################################
	# This stores length distributions
	###########################################################################
	def Store_Length_Distributions(self, l):

		# Append total time moments to tbar(l)
		self.tbarl.append(self.final_matrix.dot(self.time_transfer_ket))

		# Append state functions to Bbar(l): note the difference in final vs.
		# intermediate states
		if self.num_state_funcs > 0:
			self.Bbarl.append(self.final_B_matrix.dot(self.time_cumulative_ket) + self.int_B_matrix.dot(self.time_transfer_ket))

		# If the debug flag is set, write the entire time transfer ket to file
		if self.debug:
			self.debug_file.write(str(l) + "\t")
			for s in self.all_states:
				ind = self.Get_State_Index(s)
				self.debug_file.write(s + ",")
				moments = [str(self.time_transfer_ket[ind + i*self.num_states]) for i in range(self.max_moment + 1)]
				self.debug_file.write(",".join(moments) + "\t")
			self.debug_file.write("\n")


	###########################################################################
	# This returns true or false depending on whether the paths have converged.
	###########################################################################
	def Has_Converged(self, cur_prob, cur_max_moment, prev_max_moment):			

		# First test that total probability has converged
		if abs(1. - cur_prob) < self.converge:
			# If satisfied, check last nonzero contribution to maximum moment
			if (cur_max_moment - prev_max_moment) > 0 and (cur_max_moment - prev_max_moment)/cur_max_moment < self.converge:
				return True
	
		return False


	###########################################################################
	# This iterates the recursion relations until they reach convergence or the
	# maximum number of jumps.
	###########################################################################
	def Run(self):

		# List of jump checkpoints to output: 100, 500, 1000, 5000, ...
		checkpoints = sorted([10**i for i in range(2, 15)] + [5*(10**i) for i in range(2, 15)])
		checkpoints.reverse()

		# If in debug mode, open debug file for writing
		if self.debug:	self.debug_file = open(self.name + ".debug", "w")

		# Initialize length distributions
		if not self.no_length_dists:
			if self.length_min == 0:
				self.Store_Length_Distributions(0)

		# Keep track of the previous and current value of the maximum time
		# moment for testing convergence
		total_time_moments = self.Total_Time_Moments()
		prev_max_moment = total_time_moments[self.max_moment]

		converged = False

		# Loop over path lengths up to max_jumps
		if sys.version_info.major == 3:
			l_range = range(1, self.max_jumps + 1)
		else:
			l_range = xrange(1, self.max_jumps + 1)
		for l in l_range:

			# Print jump number l if it is a checkpoint
			if l%checkpoints[-1] == 0:
				checkpoints.pop()
				Print("\n\t" + str(l) + " jumps...")

			# Update the transfer vectors
			self.Update(l)

			# Add to the length distributions
			if not self.no_length_dists:
				if l >= self.length_min and (l - self.length_min)%self.length_res == 0:
					self.Store_Length_Distributions(l)

			total_time_moments = self.Total_Time_Moments()
			cur_prob = total_time_moments[0]
			cur_max_moment = total_time_moments[self.max_moment]

			# Check convergence
			if self.Has_Converged(cur_prob, cur_max_moment, prev_max_moment):
				converged = True
			if (not self.no_converge) and converged:
				Print("\n\t" + ("Done: converged after " + str(l) + " jumps.").ljust(64))
				return							

			prev_max_moment = cur_max_moment

		if converged:
			Print("\n\tDone: reached MAX_JUMPS and also converged.".ljust(66))
		else:
			Print("\n\tWARNING: Reached MAX_JUMPS without converging.".ljust(66))


	###########################################################################
	# This writes all data to the output files.
	###########################################################################
	def Output_Data(self):

		if self.debug:	self.debug_file.close()

		######################################################################
		# Output total moments for each final state
		######################################################################

		output_file = open(self.name + ".moments", "w")

		# Print headers
		output_file.write("# Final_state".ljust(20))
		for i in range(self.max_moment + 1): 
			output_file.write(("tbar" + str(i)).ljust(20))
		if self.calc_action:
			for i in range(self.max_moment + 1): 
				output_file.write(("sbar" + str(i)).ljust(20))
		output_file.write("\n")

		# For each final state, print time and action moments
		for s in self.final_states:
			output_file.write(s.ljust(20))
			ind = self.Get_State_Index(s)
			for i in range(self.max_moment + 1): 
				moment = self.time_cumulative_ket[ind + i*self.num_states]
				output_file.write(str(moment).ljust(20))
			if self.calc_action:
				for i in range(self.max_moment + 1): 
					moment = self.action_cumulative_ket[ind + i*self.num_states]
					output_file.write(str(moment).ljust(20))
			output_file.write("\n")

		output_file.close()


		######################################################################
		# Output distribution of time moments over path lengths
		######################################################################

		if not self.no_length_dists:
			output_file = open(self.name + ".lengths", "w")

			# Print headers
			output_file.write("# l".ljust(20))
			for i in range(self.max_moment + 1): 
				output_file.write(("tbar" + str(i)).ljust(20))
			for i in range(self.num_state_funcs): 
				output_file.write(("state_func" + str(i)).ljust(20))
			output_file.write("\n")

			# Print tbar(l) for each l
			for i in range(len(self.tbarl)):
				output_file.write(str(i*self.length_res + self.length_min).ljust(20))
				for x in self.tbarl[i]:
					output_file.write(str(x).ljust(20))
				if self.num_state_funcs > 0:
					for x in self.Bbarl[i]:
						output_file.write(str(x).ljust(20))
				output_file.write("\n")

			output_file.close()


		######################################################################
		# Output path densities over all states
		######################################################################

		total_time_moments = self.Total_Time_Moments()
		mean_time = total_time_moments[1]
		output_file = open(self.name + ".spatial", "w")
		
		# Print headers
		output_file.write("# State".ljust(20) + "Average_visits".ljust(20) + "Average_fraction_of_time" + "\n")

		# For each intermediate and final state, print average number of visits and average fraction of time
		for s in self.all_states:
			output_file.write(s.ljust(20))
			ind = self.Get_State_Index(s)

			average_visits = self.time_cumulative_ket[ind + 0*self.num_states]
			average_fraction_of_time = average_visits*self.wtm_matrices[1][ind, ind]/mean_time

			output_file.write(str(average_visits).ljust(20))
			output_file.write(str(average_fraction_of_time) + "\n")

		output_file.close()


################################################################################
# This function checks all arguments to ensure values are within range and are
# self-consistent.
################################################################################
def Check_Args(args):
	if args.name == None:
		Print("ERROR: you must specify NAME.  Exiting.\n")
		exit()
	if args.network == None:
		args.network = args.name + ".network"
	if args.bc == None:
		args.bc = args.name + ".bc"
	if args.max_moment < 1:
		Print("ERROR: Maximum moment must be at least 1.  Exiting.\n")
		exit()
	if args.max_moment > 6:
		Print("WARNING: Maximum moment is > 6, which may be memory-intensive and/or numerically unreliable.\n")
	if (args.converge < 0) or (args.converge > 1):
		Print("ERROR: Convergence criterion CONVERGE is out of range (must be between 0 and 1).  Exiting.\n")
		exit()


################################################################################
# This function checks the peak memory usage (in kb) of the current process.
################################################################################
def Get_Memory_Usage():
	status = None
	result = {'peak': 0, 'rss': 0}
	try:
		# This will only work on systems with a /proc file system (like Linux)
		status = open('/proc/self/status')
		for line in status:
			parts = line.split()
			key = parts[0][2:-1].lower()
			if key in result:
				result[key] = int(parts[1])
	finally:
		if status is not None:
			status.close()
	return result["peak"]


################################################################################
# A custom print function that does not add trailing whitespace.
################################################################################
def Print(s):
	sys.stdout.write(s)
	sys.stdout.flush()


################################################################################
# This function calculates cumulant moments from raw moments.
################################################################################
def Raw_to_Cumulant(raw_moments):
	cumulants = numpy.zeros(len(raw_moments))
	for n in range(len(cumulants)):
		total = 0.
		for i in range(1, n):
			b = scipy.special.gamma(n)/scipy.special.gamma(i)/scipy.special.gamma(n - i + 1)
			total += b*cumulants[i]*raw_moments[n - i]
		cumulants[n] = raw_moments[n] - total
	return cumulants


################################################################################
# This function calculates standardized moments from raw moments.
################################################################################
def Raw_to_Standardized(raw_moments):
	if len(raw_moments) == 1: return raw_moments
	if len(raw_moments) == 2: return numpy.array([raw_moments[0], 0.0])

	central_moments = numpy.zeros(len(raw_moments))
	for n in range(len(central_moments)):
		total = 0.
		for i in range(0, n + 1):
			b = scipy.special.gamma(n + 1)/scipy.special.gamma(i + 1)/scipy.special.gamma(n - i + 1)
			if (n-i)%2 == 0:
				total += b*raw_moments[i]*(raw_moments[1]**(n-i))
			else:
				total -= b*raw_moments[i]*(raw_moments[1]**(n-i))
		central_moments[n] = total

	sd = numpy.sqrt(central_moments[2])
	return numpy.array([central_moments[i]/sd**i for i in range(len(central_moments))])


################################################################################
# Main function that runs the algorithm.
################################################################################
def main():

	# Initialize timing
	start_date = str(datetime.datetime.now())
	start_time = time()

	# Parse command line arguments
	parser = argparse.ArgumentParser(description="PathMAN (Path Matrix Algorithm for Networks) calculates path statistics of first-passage CTRWs on networks using a transfer matrix algorithm.")
	
	# Command line arguments for run parameters: name, type, default value, and description	
	parser.add_argument("--name", 		type=str,		default=None,		action="store",		help="Name for all output files (required, no default), and default name for input files")
	parser.add_argument("--network",		type=str,		default=None,		action="store",		help="Name of network file (list of states, jump weights, waiting time moments, and state functions; default is NAME.network)")
	parser.add_argument("--bc",			type=str,		default=None,		action="store",		help="Name of boundary conditions file (initial distribution and final states; default is NAME.bc)")
	parser.add_argument("--max-moment",	type=int,		default=1,		action="store",		help="Maximum moment to calculate (default is 1)")
	parser.add_argument("--calc-action", 								action="store_true",	help="If set, the script will calculate moments of path action")	
	parser.add_argument("--num-state-funcs", type=int,	default=0,		action="store",		help="Number of state functions included in network file for which to calculate their average values at each jump (default is 0)")	
	parser.add_argument("--length-res",	type=int,		default=1,		action="store",		help="Store path length distributions at this resolution (default is 1); this can significantly reduce runtime and output file size when calculating very long paths")
	parser.add_argument("--length-min",	type=int,		default=0,		action="store",		help="Store path length distributions starting from this value (default is 0)")
	parser.add_argument("--no-length-dists",							action="store_true",	help="If set, the script will not record any distributions over path length; this can significantly reduce runtime and output file size")
	parser.add_argument("--debug",									action="store_true",	help="If set, the script will output the entire distribution of moments over path lengths (according to LENGTH_RES and LENGTH_MIN)")
	parser.add_argument("--converge", 		type=float,	default=1.e-8,		action="store",		help="Convergence criterion (default is 1e-8)")
	parser.add_argument("--max-jumps",		type=int,		default=10000000,	action="store",		help="Maximum number of jumps to calculate (default is 10000000)")
	parser.add_argument("--no-converge",								action="store_true",	help="If set, the script will calculate up to MAX_JUMPS without testing convergence")	
	#parser.add_argument("--jump-cutoff",	type=float,	default=1.e-10,	action="store",		help="The script will ignore jump probabilities below this cutoff to improve efficiency (default 1e-10)")	

	# Parse arguments and perform basic sanity checks	
	args = parser.parse_args()
	Check_Args(args)

	# Print startup information
	Print("################################################################################\n\n")
	Print("Command:".ljust(34) + "python " + " ".join(sys.argv) + "\n")
	Print("Reading network from:".ljust(34) + args.network + "\n")
	Print("Reading boundary conditions from:".ljust(34) + args.bc + "\n")
	Print("Writing output to:".ljust(34) + str(args.name) + "\n\n")

	Print("Initializing...\n")

	# Initialize network
	mypaths = Pathman(args)

	# Run path algorithm
	Print("Running...")
	mytime = time()
	mypaths.Run()
	Print("(" + str(round(time() - mytime, 3)) + " seconds)\n")

	# Output results
	Print("Writing data...")
	mytime = time()
	mypaths.Output_Data()
	Print("done.".ljust(57) + "(" + str(round(time() - mytime, 3)) + " seconds)\n\n")

	# Print total moments
	Print("Total path moments:\n")

	# Length moments
	Print("\tLength".ljust(20) + "Raw".ljust(20) + "Cumulant".ljust(20) + "Standardized".ljust(20) + "\n")
	lbars_raw = mypaths.Total_Length_Moments()
	lbars_cum = Raw_to_Cumulant(lbars_raw)
	lbars_std = Raw_to_Standardized(lbars_raw)
	for i in range(len(lbars_raw)):
		Print(("\tlbar" + str(i)).ljust(20))
		Print(str(lbars_raw[i]).ljust(20))
		Print(str(lbars_cum[i]).ljust(20))
		Print(str(lbars_std[i]).ljust(20))
		Print("\n")
	Print("\n") 

	# Time moments
	Print("\tTime".ljust(20) + "Raw".ljust(20) + "Cumulant".ljust(20) + "Standardized".ljust(20) + "\n")
	tbars_raw = mypaths.Total_Time_Moments()
	tbars_cum = Raw_to_Cumulant(tbars_raw)
	tbars_std = Raw_to_Standardized(tbars_raw)
	for i in range(len(tbars_raw)):
		Print(("\ttbar" + str(i)).ljust(20))
		Print(str(tbars_raw[i]).ljust(20))
		Print(str(tbars_cum[i]).ljust(20))
		Print(str(tbars_std[i]).ljust(20))
		Print("\n")
	Print("\n") 

	# Action moments (if calculated)
	if args.calc_action:
		Print("\tAction".ljust(20) + "Raw".ljust(20) + "Cumulant".ljust(20) + "Standardized".ljust(20) + "\n")
		sbars_raw = mypaths.Total_Action_Moments()
		sbars_cum = Raw_to_Cumulant(sbars_raw)
		sbars_std = Raw_to_Standardized(sbars_raw)
		for i in range(len(sbars_raw)):
			Print(("\tsbar" + str(i)).ljust(20))
			Print(str(sbars_raw[i]).ljust(20))
			Print(str(sbars_cum[i]).ljust(20))
			Print(str(sbars_std[i]).ljust(20))
			Print("\n")
		Print("\n")  

	# Print running time and memory usage
	Print("Run started at:\t\t" + start_date + "\n")
	Print("Run ended at:\t\t" + str(datetime.datetime.now()) + "\n")
	Print("Total run time:\t\t" + str(round(time() - start_time, 3)) + " seconds\n")
	Print("Peak memory usage:\t" + str(Get_Memory_Usage()/1.e6) + " GB\n\n")
	Print("################################################################################\n")
	

if __name__ == '__main__':
	main()
