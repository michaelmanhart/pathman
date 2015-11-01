################################################################################
#    Multidimensional Lattice Model Generator 
#    Copyright (C) 2015 Michael Manhart and Willow Kion-Crosby
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
#	Contact information:
#	E-mail: mmanhart@fas.harvard.edu
#	Mail: 
#		Harvard University
#		Department of Chemistry and Chemical Biology
#		12 Oxford Street
#		Cambridge, MA 02138, USA
################################################################################


import numpy
import scipy.misc
import random
import argparse
import math


###########################################################################
# The following function calculates the rate between 
# two locations on the lattice given the energies at those 
# sites, an effective inverse-temperature and maximum rate.
###########################################################################
def Find_Rate(beta, Ei, Ef, Gamma_0):
	
	# The Monte-Carlo transition rates are calculated 
	# using the difference in energies between the 
	# initial and final sites. beta is understood 
	# as the inverse-temperature of the system.
	deltaE = Ef-Ei
	Gamma = Gamma_0*min([1, math.exp(-deltaE*beta)])
	return Gamma
	
###########################################################################
# This functions provides a basic energy landscape given 
# the coordinate on the lattice.
###########################################################################
def Find_Energy(E_0, position):
	
	# If no energy file is found, a landscape
	# will be developed using the following function.
	
	return E_0*sum(position)

###########################################################################
# The following function calculates the nth moment of an exponential 
# distribution given the mean, tau.
###########################################################################
def Exponential_Moment(n, tau): 

	return scipy.misc.factorial(n)*(tau**n)


###########################################################################
# The following function develops a dictionary to represent 
# the energy landscape. 
###########################################################################
def Build_Energy_Landscape(L, D, E_0):

	# This function finds the energy landscape by
	# iterating across the entire lattice starting 
	# from the site 1_1_1...

	site_coordinates = numpy.ones(D)
	num_of_sites = L**D

	energy_landscape = {}

	# The step parameter keeps track of 
	# whether the function should step 
	# further in the current dimension
	# or if an edge has been reached.
	
	step = False
	while not step:

		# Label site based on coordinates.
		site = "_".join([str(int(x)) for x in site_coordinates])
		
		# Calculate energy
		energy_landscape[site] = Find_Energy(E_0, site_coordinates)
		
		# Update coordinates
		step = True
		for j in range(D):
			if step:
				site_coordinates[j] = site_coordinates[j] + 1
			if site_coordinates[j]>L:
				site_coordinates[j] = 1
				step = True
			else:
				step = False
		
	return energy_landscape

###########################################################################
# The following function finds all the neighbors of 
# a specified site. 
###########################################################################
def Find_Neighbors(site, L, D, periodic):
	
	site_coordinates = [int(x) for x in site.split("_")]
	neighbors = []
	
	# If the lattice is not periodic, neighbors beyond 1 or L
	# are not included. Otherwise each neighbor is defined as
	# the coordinates of each site a single step away from
	# the current site.
	
	for i in range(D):
		
		if site_coordinates[i]+1 <= L:
			neighbors.append("_".join([str(site_coordinates[j]+int(i == j)) for j in range(D)]))
		elif periodic:
			neighbors.append("_".join([str(site_coordinates[j]-(L-1)*int(i == j)) for j in range(D)]))
			
		if site_coordinates[i]-1 >= 1:
			neighbors.append("_".join([str(site_coordinates[j]-int(i == j)) for j in range(D)]))
		elif periodic:
			neighbors.append("_".join([str(site_coordinates[j]+(L-1)*int(i == j)) for j in range(D)]))
			
			
	return neighbors

###########################################################################
# To limit failed lattice development, this function
# checks various properties of the input parameters.
###########################################################################
def Check_Args(args):
	
	if args.lattice_size == 0:
		print("ERROR: Trying to generate a lattice of size 0.  Exiting.\n")
		exit()
		
	if args.dimension <= 0:
		print("ERROR: Dimension less than or equal to 0.  Exiting.\n")
		exit()
		
	if args.start_sites == None:
		args.start_sites = ",".join(["_".join([str(1) for d in range(args.dimension)])])
		
	if args.final_sites == None:
		args.final_sites = ",".join(["_".join([str(args.lattice_size) for d in range(args.dimension)])])
		
	if args.energy_file == None:
		args.energy_file = str(args.dimension)+"D_lattice"

def main():

	# Read user input lattice parameters
	
	parser = argparse.ArgumentParser(description="This script generates a multi-dimensional Euclidean lattice with a general energy landscape (default is a linear potential).")
	
	parser.add_argument("--lattice-size",		type=int,	default=10,	action="store",		help="Number of positions for each dimension  (default 10)")
	parser.add_argument("--dimension",		type=int, 	default=2,	action="store",		help="Dimension of Euclidean lattice (default 2)")
	parser.add_argument("--max-moment",		type=int, 	default=1,	action="store",		help="Maximum moment included for waiting time distributions at each site (default 1)")
	parser.add_argument("--start-sites",		type=str, 			action="store",		help="Location of each start site on the lattice written as x1_x2_...xD with xi between 1 and L, represented as a list of strings e.g. 1_1,1_2,1_3 (default 1_1_1_..._1). If more than one site is listed, each starting position will be given an equal probability.")	
	parser.add_argument("--final-sites",		type=str, 			action="store",		help="Location of each final site on the lattice written as x1_x2_...xD with xi between 1 and L, represented as a list of strings e.g. 8_10,9_10,10_10 (default L_L_L_..._L)")	
	parser.add_argument("--beta",			type=float, 	default=0.0,	action="store",		help="Effective temperature (default 0)")
	parser.add_argument("--E-scale",		type=float, 	default=1.0,	action="store",		help="Energy scale of landscape (default 1)")
	parser.add_argument("--rate-coefficient",	type=float, 	default=1.0,	action="store",		help="Rate coefficient (default 1)")
	parser.add_argument("--energy-file",		type=str, 			action="store",		help="Name of energy landscape file that will be used to compute rates between sites on the lattice. If the file does not exist, it will be generated as a linear potential proportional to the sum of the coordinates at each position. File extension is '.energy' (default [specified dimension]D_lattice.energy)")
	parser.add_argument("--periodic-boundary", 					action="store_true",	help="Set if the lattice has periodic boundary conditions")
	
	args = parser.parse_args()
	
	# Check validity of input and initialize parameters
	
	Check_Args(args)
	
	# Lattice properties:
	
	beta = args.beta
	E_0 = args.E_scale
	D = args.dimension
	L = args.lattice_size
	Gamma_0 = args.rate_coefficient
	periodic = args.periodic_boundary
	
	max_moment = args.max_moment
	file_name = args.energy_file
	
	start_sites = args.start_sites.split(",")
	final_sites = args.final_sites.split(",")
	
	
	# Check if .energy file already exists. If so, read-in energy dictionary,
	# if not, generate .energy file. A list of sites is also found.
	
	energy_landscape = {}

	try:
		energy_file = open(file_name + ".energy", "r")
		for line in energy_file.readlines():
			energy_landscape[line.split()[0]] = float(line.split()[1])
		
		sites = energy_landscape.keys()
		sites.sort()

	except IOError:
	
		energy_file = open(file_name + ".energy", "w")
		energy_landscape = Build_Energy_Landscape(L, D, E_0)
		
		sites = energy_landscape.keys()
		sites.sort()
		
		for site in sites:
			energy_file.write(site+"\t"+str(energy_landscape[site])+"\n")
		
	energy_file.close()	
	
	# Using the energy landscape, the rates from each site to its 
	# neighbors are found. The average time spent in each site
	# is then computed via the sum of the rates out of each site.
	# These values are then written to the .network file. 
	
	network_file = open(file_name + ".network", "w")
	
	for site in sites:
		
		# Find a list of number for each site
		
		neighbors = Find_Neighbors(site, L, D, periodic)
		transitions = []
		rate_total = 0.0
		
		# For each neighbor, calculate the rate to that 
		# neighbor via the energy change of the
		# transition and the inverse temperature.
		
		for neighbor in neighbors:
			
			Ei = energy_landscape[site]
			Ef = energy_landscape[neighbor]
		
			rate = Find_Rate(beta, Ei, Ef, Gamma_0)
			rate_total = rate_total + rate
			
			transitions.append(neighbor+","+str(rate))
		
		# Assuming exponentially-distributed waiting time
		# the specified number of moments are found.
		
		site_moments = ",".join([str(Exponential_Moment(n, 1./rate_total)) for n in range(1,max_moment + 1)])
		
		# The coordinates of each site is also included in the .network
		# file as a state dependent function for each dimension.
		
		site_coordinates = ",".join([x for x in site.split("_")])
		
		network_file.write(site + "\t" + ";".join(transitions))
		network_file.write("\t" + site_moments)
		network_file.write("\t" + site_coordinates)
		network_file.write("\n")
		
	network_file.close()

	# The boundary conditions file (.bc) is
	# written such that all starting sites have 
	# equal probability.

	bc_file = open(file_name + ".bc", "w")
	part_func = len(start_sites)
	
	bc_file.write("\t".join([start_site+","+str(1./part_func) for start_site in start_sites]))
	bc_file.write("\n")
	bc_file.write("\t".join(final_sites))
	bc_file.close()

if __name__ == '__main__':
	main()
