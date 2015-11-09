################################################################################
#    Plotting Script for Landscape and Properties of the RBM
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


from pylab import colorbar
import math
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
import argparse

pylab.rc('font', family='serif', size=14)


def Check_Args(args):
		
	if args.name == None:
		args.name = str(args.dimension)+"D_lattice"

def main():

	# Read user input lattice parameters
	
	parser = argparse.ArgumentParser(description="This script generates a 2 dimensional plot to visualize various properties of the RBM.")
	
	parser.add_argument("--name",		type=str, 			action="store",		help="Name of energy landscape file that will be visualized. The .network, .spatial, and .lengths files are assumed to start with the same string as the .energy file.")
	#parser.add_argument("--calc-boundary", 					action="store_true",	help="Set if the lattice has periodic boundary conditions")
	
	args = parser.parse_args()
	
	# Check validity of input
	
	Check_Args(args)
	
	time_per_site = {}
	visits_per_site = {}
	ave_wait_time = {}

	L = 10

	calc_action = False
	max_moment = 1

	fig = plt.figure(num=None, figsize=(16,4), facecolor='white')

	ave_path_color = "magenta"	
	cmap = plt.cm.jet
	cmap_energy = plt.cm.gist_gray

	line_width = 4

	letter_size = 30
	beta_size = 24
	num_size = 16
	power_size = 16
	title_size = 22
	x_y_size = 20

	file_name = args.name

	lengths = True
	pbcheck = False
	len_res = 10
	num_state_func = 2

	print('Searching for data file '+file_name+".network"+'... ')
	
	myFile = open(file_name+".network", "r")

	print('Extracting data... ')

	rates = {}
	

	for line in myFile.readlines():
		if not (line.startswith('#') or line.isspace()):	 	
			for nn in line.split()[1].split(";"):
				neighbor = nn.split(",")[0]
				rate = float(nn.split(",")[1])
				rates[line.split()[0]+"->"+neighbor] = rate
			ave_wait_time[line.split()[0]] = float(line.split()[2].split(",")[0])
		
	myFile.close()

	myFile = open(file_name+".energy", "r")

	energies = {}

	for line in myFile.readlines():
		 if not (line.startswith('#') or line.isspace()):
		 	energies[line.split()[0]+'-'+line.split()[1]] = float(line.split()[2])
		
	myFile.close()

	myFile = open(file_name+".spatial", "r")

	for line in myFile.readlines():
		if not (line.startswith('#') or line.isspace()):
			site_element = line.split()
			time_per_site[site_element[0]] = float(site_element[2])
			visits_per_site[site_element[0]] = float(site_element[1])

	myFile.close()



	try:
	
		myFile = open(file_name+".lengths", "r")
		lengths_file = True
	except IOError:
	
		lengths_file = False
	
	if lengths_file:	
		average_path_y = []
		average_path_x = []
		average_index = []
		res = 1
		res_count = res
		for line in myFile.readlines():
			if not (line.isspace() or line.startswith('#')):
				path_element = line.split()
	
				if res_count == res:
					res_count = 0
					average_path_x.append(float(path_element[(max_moment+1)*(1+int(calc_action))+1])-1.)
					average_path_y.append(float(path_element[(max_moment+1)*(1+int(calc_action))+2])-1.)
	
				res_count += 1

				if num_state_func>2:
					average_index.append(float(path_element[8]))
	
		myFile.close()

	print('Complete')

	time_per_site_values = [value for value in time_per_site.values()]
	visits_per_site_values = [value for value in visits_per_site.values()]
	ave_wait_time_values = [value for value in ave_wait_time.values()]
	rates_values = [value for value in rates.values()]
	energies_values = [value for value in energies.values()]

	max_time = max(time_per_site_values)
	min_time = min(time_per_site_values)
	max_visits = max(visits_per_site_values)
	min_visits = min(visits_per_site_values)
	max_awt = max(ave_wait_time_values)
	min_awt = min(ave_wait_time_values)
	max_rate = max(rates_values)
	min_rate = min(rates_values)
	max_energy = max(energies_values)
	min_energy = min(energies_values)

	mean_rate = sum(rates_values)/len(rates_values)
	mean_energy = sum(energies_values)/len(energies_values)
	mean_time = sum(time_per_site_values)/len(time_per_site_values)

	time_grid = [[(time_per_site[str(x)+"_"+str(y)]) for x in range(1, L+1)] for y in range(1, L+1)]
	log_time_grid = [[math.log(time_per_site[str(x)+"_"+str(y)])/math.log(10.) for x in range(1, L+1)  if not (time_per_site[str(x)+"_"+str(y)] == 0.0)] for y in range(1, L+1)]
	visits_grid = [[visits_per_site[str(x)+"_"+str(y)] for x in range(1, L+1)] for y in range(1, L+1)]
	log_visits_grid = [[math.log(visits_per_site[str(x)+"_"+str(y)])/math.log(10.) for x in range(1, L+1)] for y in range(1, L+1)]
	awt_grid = [[(ave_wait_time[str(x)+"_"+str(y)]) for x in range(1, L+1)] for y in range(1, L+1)]
	log_awt_grid = [[math.log(ave_wait_time[str(x)+"_"+str(y)])/math.log(10.) for x in range(1, L+1)] for y in range(1, L+1)]

	plt3 = plt.subplot(1,3,1)

	#plt.text(8.5, -30, 'Height of energy barriers:', fontsize=num_size)

	plt.xticks([0,1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9,10])
	plt.yticks([0,1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9,10])

	plt.xlim([-.5, L-.5])
	plt.ylim([-.5, L-.5])	

	max_log_awt = max([max(log_awt_grid[i]) for i in range(len(log_awt_grid))])
	min_log_awt = min([min(log_awt_grid[i]) for i in range(len(log_awt_grid))])
	
	log_scale = False

	plt.title('Mean waiting time '+r'$\theta^{(1)}(x,y)$'+'\n')

	if log_scale:
		cs_awt = plt.imshow(log_awt_grid, interpolation='none', cmap=cmap)
		tick_points = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000., 20000., 30000., 40000., 50000., 60000., 70000., 80000., 90000., 100000., 200000., 300000., 400000., 500000., 600000., 700000., 800000., 900000., 1000000.]
		log_tick_points = [math.log(tick, 10.) for tick in tick_points]
		cbar = colorbar(cs_awt, ticks=log_tick_points, fraction=0.045)
		plt.clim(min_log_awt, max_log_awt)
		cbar.ax.set_yticklabels([(r"$10^{%s}$" % str(int(round(tick,4))))*int(round(abs(tick),4) == float(math.ceil(abs(tick)))) for tick in log_tick_points if int(round(tick,4)) > min_log_awt and int(round(tick,4)) < max_log_awt], [tick for tick in log_tick_points], visible=True, fontsize=num_size)
	else:
		cs_awt = plt.imshow(awt_grid, interpolation='none', cmap=cmap)
		tick_points = np.linspace(min_awt, max_awt, 6, endpoint=True)
		plt.clim(min_awt, max_awt)
		cbar = colorbar(cs_awt, ticks=tick_points, fraction=0.045)
		cbar.formatter.set_powerlimits((0, 0))

	for i in range(0, L+1):
		rec_width = 0.2
		plt.gca().add_patch(Rectangle((i - .5 - 0.5*rec_width, -0.5), rec_width, L+1, color="black"))
		plt.gca().add_patch(Rectangle((-0.5, i - .5 - 0.5*rec_width), L+1, rec_width, color="black"))

	for i in range(1, L):

		for j in range(1, L+1):
			energy = energies[str(i)+"_"+str(j)+"-"+str(i+1)+"_"+str(j)]
			E_eff = (energy-min_energy)/(max_energy-min_energy)
			plt.gca().add_patch(Rectangle((i - 0.5 - 0.5*rec_width, j - 1.5 + rec_width), rec_width, 1-2*rec_width, color="white", alpha=E_eff))
		
			energy = energies[str(j)+"_"+str(i)+"-"+str(j)+"_"+str(i+1)]
			E_eff = (energy-min_energy)/(max_energy-min_energy)
			plt.gca().add_patch(Rectangle((j - 1.5 + rec_width, i - .5 - 0.5*rec_width), 1-2*rec_width, rec_width, color="white", alpha=E_eff))
	

	if lengths:
		plt.plot(average_path_x, average_path_y, lw=2., color=ave_path_color)
		plt.plot(average_path_x, average_path_y, 'mo', lw=line_width/2.)


	plt.xlabel(r'$x$', size=x_y_size)

	plt.ylabel(r'$y$', rotation=0, size=x_y_size, position=(0.5,.45))

	plt2 = plt.subplot(1,3,2, adjustable='box', aspect=1.0)


	plt.xticks([0,1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9,10])
	plt.yticks([0,1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9,10])

	plt.xlim([-.5, L-.5])
	plt.ylim([-.5, L-.5])	

	max_log_visits = max([max(log_visits_grid[i]) for i in range(len(log_visits_grid))])
	min_log_visits = min([min(log_visits_grid[i]) for i in range(len(log_visits_grid))])

	plt.title('Average visits '+r'$v(x,y)$'+'\n')

	if log_scale:
		cs_visits = plt.imshow(log_visits_grid, interpolation='none', cmap=cmap)
		tick_points = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000., 20000., 30000., 40000., 50000., 60000., 70000., 80000., 90000., 100000., 200000., 300000., 400000., 500000., 600000., 700000., 800000., 900000., 1000000.]
		log_tick_points = [math.log(tick, 10.) for tick in tick_points]
		cbar = colorbar(cs_visits, ticks=log_tick_points, fraction=0.045)
		plt.clim(min_log_visits, max_log_visits)
		cbar.ax.set_yticklabels([(r"$10^{%s}$" % str(int(round(tick,4))))*int(round(abs(tick),4) == float(math.ceil(abs(tick)))) for tick in log_tick_points if int(round(tick,4)) > min_log_visits and int(round(tick,4)) < max_log_visits], [tick for tick in log_tick_points], visible=True, fontsize=num_size)
	else:
	
		cs_visits = plt.imshow(visits_grid, interpolation='none', cmap=cmap)
		tick_points = np.linspace(min_visits, max_visits, 6, endpoint=True)
		plt.clim(min_visits, max_visits)
		cbar = colorbar(cs_visits, ticks=tick_points, fraction=0.045)
		cbar.formatter.set_powerlimits((0, 0))

	for i in range(0, L+1):
		rec_width = 0.2
		plt.gca().add_patch(Rectangle((i - .5 - 0.5*rec_width, -0.5), rec_width, L+1, color="black"))
		plt.gca().add_patch(Rectangle((-0.5, i - .5 - 0.5*rec_width), L+1, rec_width, color="black"))

	for i in range(1, L):

		for j in range(1, L+1):
			energy = energies[str(i)+"_"+str(j)+"-"+str(i+1)+"_"+str(j)]
			E_eff = (energy-min_energy)/(max_energy-min_energy)
			plt.gca().add_patch(Rectangle((i - 0.5 - 0.5*rec_width, j - 1.5 + rec_width), rec_width, 1-2*rec_width, color="white", alpha=E_eff))
		
			energy = energies[str(j)+"_"+str(i)+"-"+str(j)+"_"+str(i+1)]
			E_eff = (energy-min_energy)/(max_energy-min_energy)
			plt.gca().add_patch(Rectangle((j - 1.5 + rec_width, i - .5 - 0.5*rec_width), 1-2*rec_width, rec_width, color="white", alpha=E_eff))


	if lengths_file:
		plt.plot(average_path_x, average_path_y, lw=line_width/2., color=ave_path_color)
		plt.plot(average_path_x, average_path_y, 'mo', lw=line_width/2.)


	
	plt.xlabel(r'$x$', size=x_y_size)

	plt.ylabel(r'$y$', rotation=0, size=x_y_size, position=(0.5,.45))

	#position=fig.add_axes([0.475, 0.02, 0.168, 0.012])
	
	tick_points = [float(i) for i in range(10+1)]
	

	plt1 = plt.subplot(1,3,3, adjustable='box', aspect=1.0)

	plt.xticks([0,1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9,10])
	plt.yticks([0,1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9,10])

	plt.xlim([-.5, L-.5])
	plt.ylim([-.5, L-.5])	

	cs_time = plt.imshow(time_grid, interpolation='none', cmap=cmap)

	for i in range(0, L+1):
		rec_width = 0.2
		plt.gca().add_patch(Rectangle((i - .5 - 0.5*rec_width, -0.5), rec_width, L+1, color="black"))
		plt.gca().add_patch(Rectangle((-0.5, i - .5 - 0.5*rec_width), L+1, rec_width, color="black"))

	for i in range(1, L):

		for j in range(1, L+1):
			energy = energies[str(i)+"_"+str(j)+"-"+str(i+1)+"_"+str(j)]
			E_eff = (energy-min_energy)/(max_energy-min_energy)
			plt.gca().add_patch(Rectangle((i - 0.5 - 0.5*rec_width, j - 1.5 + rec_width), rec_width, 1-2*rec_width, color="white", alpha=E_eff))
	
			energy = energies[str(j)+"_"+str(i)+"-"+str(j)+"_"+str(i+1)]
			E_eff = (energy-min_energy)/(max_energy-min_energy)
			plt.gca().add_patch(Rectangle((j - 1.5 + rec_width, i - .5 - 0.5*rec_width), 1-2*rec_width, rec_width, color="white", alpha=E_eff))


	if lengths:
		plt.plot(average_path_x, average_path_y, lw=line_width/2., color=ave_path_color)
		plt.plot(average_path_x, average_path_y, 'mo', lw=line_width/2.)

	plt.clim(min_time, max_time)
	tick_points = np.linspace(min_time, max_time, 6, endpoint=True)
	cbar = colorbar(cs_time, ticks=tick_points, fraction=0.045)
	cbar.formatter.set_powerlimits((0, 0))
	#cbar.ax.set_yticklabels(["0.0", "0.3", "0.6", "0.9", "1.2", "1.5", "1.8"+r'$\times 10^{-2}$'], visible=True)

	plt.title('Fraction of time '+r'$\theta^{(1)}(x,y)v(x,y)/\bar{t}^{(1)}$'+'\n')
	
	plt.xlabel(r'$x$', size=x_y_size)

	plt.ylabel(r'$y$', rotation=0, size=x_y_size, position=(0.5,.45))
		
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.3)

	plt.show()	

if __name__ == '__main__':
	main()
