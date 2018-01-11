######!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to plot Pecube results in forward mode
By Xavier Robert
Grenoble, 2015.06.01

USAGE :
  1- Copy this file in the main Pecube/ folder where run.sh is stored
  2- Add the line "python Plot-Pecube-Forward.py at the end of the file run.sh
  3- Set which system you want to plot in the "# Define data to analysis" section (see below)
  4- Run ./run.sh
  2to4bis - Run in the terminal: $ python Plot-Pecube-Forward.py after the run of ./bin/Pecube
  5- The plots will be in the folder RUNXX/Graphs/ where RUNXX/ is the working folder
               set in topo_parameter.txt

INPUTS:
The inputs are in the script file, in the "# Define data to analysis" section. 
The different arguments are described.

xavier.robert@ujf-grenoble.fr

(c) licence CCby-nc : http://creativecommons.org/licenses/by-nc/3.0/ 2015

"""

###### To DO :  #######
#    - Test the MTL pdfs (I have not tested it, even if it should be OK)
#    - Add error bars on data graphs ? 
#      For that we should read the input data file given in topo_parameters.txt 
#                    (paragraphe N°12 of topo_parameters.txt)
#      
###### End To DO #######

from __future__ import  division
# This to be sure that the result of the division of integers is a real, not an integer

# Import modules
import sys
import os
import copy
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import shutil


#######################################
# Define data to analysis
# dataplot: List of data to plot ; By default, the altitude will be plotted
#           Do not forget the simple quotes !!!
#           ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL']
dataplot = ['AFT', 'ZFT']

# data: Name of the file Comparison.txt produced by Pecube
#       Usually you do not have to change it
data = 'Comparison.txt'

# graphpath: name of the folder where the plot will be written
#            Usually you do not have to change it
graphpath = 'Graphs'

# profiletype: type of profile = ['Latitude', 'Longitude', 'Altitude']
#              if [], no age profile is plotted
#profiletype = ['Latitude', 'Longitude', 'Altitude']
profiletype = ['Latitude', 'Altitude']

# end define the data and parameters
#######################################


# Dictionaries definitions:
# You may have to change it depending on your Pecube version, and on your settings...
# agecol: respective column number of the data in the file comparison.txt
#         needs to be changed depending on your settings in the last paragraph 
#         of the file topo_parameters.txt
agecol = {'alt' : 3,
          'AHe' : 5,
          'AFT' : 7,
          'ZHe' : 9,
          'ZFT' : 11,
          'KAr' : 13,
          'BAr' : 15,
          'MAr' : 17,
          'HbAr' : 19,
          'FTL' : 21
          }
# errcol: respective column number of the error on data in the data input file
errcol = {'AHe' : 5,
          'AFT' : 7,
          'ZHe' : 9,
          'ZFT' : 11,
          'KAr' : 13,
          'BAr' : 15,
          'MAr' : 17,
          'HbAr' : 19,
          'FTL' : 21
          }
# agename: legend of each data system         
agename = {'AHe' : 'AHe (Ma)',
           'AFT' : 'AFT (Ma)',
           'ZHe' : 'ZHe (Ma)',
           'ZFT' : 'ZFT (Ma)',
           'KAr' : 'KAr (Ma)',
           'BAr' : 'Biot. Ar (Ma)',
           'MAr' : 'Musc. Ar (Ma)',
           'HbAr' : 'Hb Ar (Ma)',
           'FTL' : 'FT length (µm)'
           }
# predname: legend of each predicted system         
predname = {'AHe' : 'Predicted AHe (Ma)',
            'AFT' : 'Predicted AFT (Ma)',
            'ZHe' : 'Predicted ZHe (Ma)',
            'ZFT' : 'Predicted ZFT (Ma)',
            'KAr' : 'Predicted KAr (Ma)',
            'BAr' : 'Predicted Biot. Ar (Ma)',
            'MAr' : 'Predicted Musc. Ar (Ma)',
            'HbAr' : 'Predicted Hb Ar (Ma)',
            'FTL' : 'Predicted FT length (µm)'
           }
# plotprof: respective column number of the data in the file comparison.txt
plotprof = {'Latitude' : 1, 
            'Longitude' : 0, 
            'Altitude' : 2
            }
# colores: Colors used for the different age system
colores = {'AHe' : 'y',
           'AFT' : 'r',
           'ZHe' : 'g',
           'ZFT' : 'b',
           'KAr' : 'k',
           'BAr' : 'c',
           'MAr' : 'm',
           'HbAr' : '0.75',
           'FTL' : 'y'
           }

############################################################################
######## Beginning of functions definition ###########
def plottot(datac, i, xlabel, ylabel, err = None):
	"""
	Function to plot the data
	"""
	if err == None:
		plt.plot(datac[:, i-1], datac[:, i], 'o', color = 'r')
	else:
		plt.errorbar(datac[:, i-1], datac[:, i], xerr = err, fmt = 'o', color = 'r')
		
	# Find axes limits to avoid the values < 0 which are NaN data
	if min(datac[:, i-1]) >= 0 and min(datac[:, i]) >= 0:
		# Check if there are NaN values
		# If not :
		minx = min(np.min(datac[:, i-1]), np.min(datac[:, i]))
		maxx = max(np.nanmax(datac[:, i-1]), np.nanmax(datac[:, i]))
	else:
		# If yes :
		# find min and max without NaN data
		datanan = copy.copy(datac)
		datanan[datanan[:, i-1] < 0.] = 'NaN'
		datanan[datanan[:, i] < 0.] = 'NaN'
		minx = min(np.nanmin(datanan[:, i-1]), np.nanmin(datanan[:, i]))
		maxx = max(np.nanmax(datanan[:, i-1]), np.nanmax(datanan[:, i]))
	
	# plot the line 1:1 for comparison
	plt.plot([minx,maxx], [minx,maxx], '-', color = 'b')
	plt.xlim(xmin = minx, xmax = maxx)
	plt.ylim(ymin = minx, ymax = maxx)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.grid(True)
	# Simplify the ticks (if not this is hard to read...)
	if i == 3:
		if (maxx-minx)/400 < 1:
			plt.xticks(np.arange(np.round(minx/10) * 10, 
			                    (np.round(maxx/10) +1) * 10, 
			                    np.round((maxx-minx)/40) * 10), 
			           rotation = 45)
			plt.yticks(np.arange(np.round(minx/10) * 10, 
			                     (np.round(maxx/10) +1) * 10, 
			                     np.round((maxx-minx)/40) * 10))
		elif (maxx-minx)/4000 < 1:
			plt.xticks(np.arange(np.round(minx/100) * 100, 
			           (np.round(maxx/100) +1) * 100, 
			            np.round((maxx-minx)/400) * 100), 
			            rotation = 45)
			plt.yticks(np.arange(np.round(minx/100) * 100, 
			                     (np.round(maxx/100) +1) * 100, 
			                     np.round((maxx-minx)/400) * 100))
		elif (maxx-minx)/40000 < 1:
			plt.xticks(np.arange(np.round(minx/1000) * 1000, 
			                     (np.round(maxx/1000) +1) * 1000, 
			                      np.round((maxx-minx)/4000) * 1000), 
			                      rotation = 45)
			plt.yticks(np.arange(np.round(minx/1000) * 1000, 
			                     (np.round(maxx/1000) +1) * 1000, 
			                      np.round((maxx-minx)/4000) * 1000))	
	else:
		if maxx-minx < 4:
			inter = 1
			if maxx-minx < 2:
				inter = 0.5
		else:
			inter = int((maxx-minx)/4)
		plt.xticks(np.arange(int(minx), int(maxx) + 1, inter))
		plt.yticks(np.arange(int(minx), int(maxx) + 1, inter))
	
	return
	# End of the function

############################################################################

def plotmtl(datac, iii, system, agecol, pgraph, pfolder):
	"""
	Function to plot the MTL data
	"""
	
	# open a new figure
	#plt.figure(iii+3)
	lrange = np.arange(0.5,17.5)
	mtl = datac[iii,agecol[system]:agecol[system]+17]
	plt.bar(lrange, mtl, 
			widtch = 1, 
			bottom = None, 
			label='MTL',
			color = 'r', 
			alpha=0.5)
	mtlp = datac[iii,agecol[system]+18:agecol[system]+18+17]
	plt.bar(lrange, mtlp, 
			widtch = 1, 
			bottom = None, 
			label='Predicted MTL',
			color = 'b', 
			alpha=0.5)
	plt.xlable(agecol[system])
	plt.ylabel('Frequency')
	plt.legend(loc='best', numpoints = 1)
	plt.set_title('Sample number ' + str(iii + 1))
	
	return
	# End of the function

############################################################################

def find_nbplot(nbplot):
	"""
	Function to find the number line/columns of subplot 
	      depending on the total number of plots
	      
	      Written by Xavier Robert 2015.06.11
	      
	INPUT:
	   nbplot = total number of plots   
	OUTPUT:
	    nbplotx = nb of lines in the subplot
	    nbploty = nb of columns in the subplots
	USAGE:
	  nbplotx, nbploty = find_nbplot(nbplot)
	"""
	
	if nbplot == 1:
		nbplotx = 1
		nbploty = 1
	elif nbplot == 2:
		nbplotx = 2
		nbploty = 1
	elif nbplot in [3, 4]:
		nbplotx = 2
		nbploty = 2
	elif nbplot in [5, 6]:
		nbplotx = 3
		nbploty = 2
	elif nbplot in [7, 9]:
		nbplotx = 3
		nbploty = 3		
	elif nbplot in [10, 12]:
		nbplotx = 4
		nbploty = 3
	else:
		# If there more plots to do, it has not be yet implemented
		# Raise an error and give the clues to solve it !
		try:
			raise Exception("foo")
		except:
			for frame in traceback.extract_tb(sys.exc_info()[2]):
				fname,lineno,fn,text = frame
				sys.exit('ERROR: The number of parameters is more than 12 (in fatc %i). \n'
			         'You need to change the code before line %d' % (nbplot, lineno))
			         
	return nbplotx, nbploty

######## END of the functions definition ###########
############################################################################

############################################################################
######### Main code ###############
if __name__ == "__main__":
	#### Read the data
	#print(' ')
	print('___________________________________________')
	print(' ')
	print('Plotting data for a Pecube forward model...')
	print('___________________________________________')
	
	# Read in which folder we are wrking from topo_parameters.txt
	f0r = open('input/topo_parameters.txt', "r")
	inputd = []
	for line in f0r:
		if line[0] != '$' and line[0] != '#' and line[0] != ' ' \
		   and len(line) != 1 and len(line) != 0 :
			inputd.append(line.strip())
	# close the file
	f0r.close()
	# get the working folder (first input in topo_parameters.txt)
	pfolder = inputd[0]
	# get the ages data file (3rd input from the end of the file topo_parameters.txt ; 
	#     -3 should be changed if the end of topo_parameters.txt is changed)
	datain = inputd[-3]
	
	print('   Working in the folder: ' + str(pfolder))
	pgraph = str(pfolder) + '/' + graphpath
	# Check if the Graphs/ folder exists, if not create it
	if os.path.exists(pgraph) == False:
		print('Output folder does not exist...')
		print('I am creating it...')
		print(' ')
		os.mkdir(pgraph)
	else:
		print('   Folder ' + pgraph +'/ already exists, I will write in it and erase previous pdf file')
	# Read file Comparison.txt
	pdata = pfolder + '/' + data
	# Check if file Comparison.txt exists, if not, raise an error
	if os.path.isfile(pdata) == False and os.access(pdata, os.R_OK) == False :
		#print('ERROR : File {FileNa} does not exist'.format(FileNa=str(pdata))) 
		sys.exit('ERROR : File {FileNa} does not exist'.format(FileNa=str(pdata)))
	print('    Reading data...')
	datac = np.loadtxt(pdata, skiprows = 1)
	if datain != 'Nil':
		if os.path.isfile('Data/'+ datain) == False and os.access('Data/'+ datain, os.R_OK) == False :
			sys.exit('ERROR : File {FileNa} does not exist'.format(FileNa=str('Data/'+ datain)))
		dataerr = np.loadtxt('Data/'+ datain, skiprows = 1)
		dataerr[dataerr[:,:] < 0.] = 'NaN'
	else:
		dataerr[0:datac.shape[0], 0:datac.shape[1]]	= 'NaN'
	
	if datain != 'Nil':
		print('    Plotting data')
		# Determine the number of plot
		nbplot = int(len(dataplot)) + 1
		nbplotx, nbploty = find_nbplot(nbplot)
		# Open the new figure
		plt.figure(1)
		# Adjust the subplots
		plt.subplots_adjust(wspace = (nbplot + 2)/10., hspace = (nbplot + 2)/10.)
		
		#### Plot the Altitude comparison
		plt.subplot(nbploty, nbplotx, 1, aspect='equal')
		plottot(datac, agecol['alt'], 'Elevation (m)', 'Predicted elevation (m)', err = None)
		
		for system in dataplot:
			if system != 'FTL':
				print('    Plotting ' + system + ' system')
				plt.subplot(nbploty, nbplotx, dataplot.index(system) + 2, aspect='equal')
				# use the dictionary defined at th beginning of the code
				plottot(datac, agecol[system], agename[system], predname[system], err = dataerr[:, errcol[system] - 1])												
		# Add main title
		plt.suptitle('Run ' + str(pfolder), size=16)
		plt.savefig(pgraph + '/Plot_forward' + pfolder +'.pdf')
		plt.close(1)
		
		if 'FLT' in dataplot: 
			# Normalized FTL distribution as an array of length 17. 
			# ftld(k) is the proportion of tracks with length between k-0.5 and k+0.5 µm
			# the sum of ftlf(k) for k = 1 to 17 must be equal to 1
			# needs to be done for each given data > 0
			plt.figure(2)
			# Find the number of data to plot
			i = 0
			for iii in range (0, nbplot-1):
				if datac[iii, agecol[system]] >= 0:
					i += 1
			nbplotx, nbploty = find_nbplot(nbplot)
	
			# and only them define the subplots
			plt.subplots_adjust(wspace = (nbplot + 2)/10., hspace = (nbplot + 2)/10.)
			print('    Plotting FTL distributions')
			for iii in range (0, nbplot-1):
				plt.subplot(nbploty, nbplotx, iii + 1, aspect='equal')
				# stepping over data
				if datac[iii, agecol[system]] >= 0:
					print('       Doing data number ' + str(iii) + ' on ' + str(nbplot))
					# if there are FTL datan call plotmtl
					plotmtl(datac, iii, system, agecol, pgraph, pfolder)
				else:
					print('       Data number %i on %i without MTL data' % (iii + 1, nbplot))
			plt.suptitle('MTL distributions')
			plt.savefig(pgraph + '/MTL_' + pfolder +'.pdf')
			plt.close(2)
		print(' ')
	else:
		print('No input data, no data/prediction graphs')
	
	# Plot data and prediction along simple transects (Altitude, Longitude or Latitude
	for plotplot in profiletype:
		# The next few line are used to plot ages along a transect to compare data and model predictions. 
		# It could be changed, depending on your model settings
		# plot ages in function of latitude
		print('    Plotting age transects %s...' % plotplot)
		plt.figure()
		if plotplot != 'Altitude':
			for system in dataplot:
				if system != 'FTL':
					datanan = copy.copy(datac)
					if datain != 'Nil': 
						plt.errorbar(datac[:, plotprof[plotplot]], datac[:, agecol[system]-1],
						             yerr =  dataerr[:, errcol[system] - 1],
						             fmt = 'o', 
						             color = colores[system], 
						             label = agename[system])
						# remove predictions where there is no data
						datanan[datanan[:, agecol[system]-1] < 0.] = 'NaN'
						for i in range (0,datanan.shape[0]):
							if datanan[i, agecol[system] - 1] == 'NaN':
								datanan[i, agecol[system]] == 'NaN'
								
					plt.plot(datac[:, plotprof[plotplot]], datanan[:, agecol[system]], 
					               's', 
					               color = colores[system], 
					               label = predname[system], 
					               alpha = 0.5)
			
		else:
			for system in dataplot:
				if system != 'FTL':
					datanan = copy.copy(datac)
					if datain != 'Nil':
						plt.errorbar(datac[:, agecol[system] - 1], datac[:, plotprof[plotplot]],
						             xerr =  dataerr[:, errcol[system] - 1],
						             fmt = 'o', 
						             color = colores[system], 
						             label = agename[system])
						# remove predictions where there is no data
						datanan[datanan[:, agecol[system]-1] < 0.] = 'NaN'
						for i in range (0,datanan.shape[0]):
							if datanan[i, agecol[system] - 1] == 'nan':
								datanan[i, agecol[system]] == 'NaN'
								
					plt.errorbar(datac[:, agecol[system]], datanan[:, plotprof[plotplot] + 1],
					         yerr = np.abs(datac[:, agecol['alt'] - 1] - datac[:, agecol['alt']]),
					         fmt = 's', 
					         color = colores[system], 
					         label = predname[system], 
					         alpha = 0.5)

		plt.grid(True)
		plt.legend(loc='best', numpoints = 1)
		# Find the age range. For that, choose 0 for the minimum, 
		#      and the max age whatever is the system used for the max range
		maxage = max(max(datac[:, agecol['AHe']]), 
		             max(datac[:, agecol['AFT']]), 
		             max(datac[:, agecol['ZFT']]),
		             max(datac[:, agecol['ZHe']]), 
		             max(datac[:, agecol['KAr']]), 
		             max(datac[:, agecol['MAr']]), 
		             max(datac[:, agecol['BAr']]), 
		             max(datac[:, agecol['HbAr']]))
		if plotplot != 'Altitude':
			# If the range of ages is not goof on the graph, it can be set manually here
			# plt.xlim(0, your maxage)
			plt.ylim(0, maxage)
			plt.xlabel('%s ($^o$N)' % plotplot)
			plt.ylabel('Age (Ma)')
		else:
			# If the range of ages is not goof on the graph, it can be set manually here
			# plt.xlim(0, your maxage)
			plt.xlim(0, maxage)
			plt.ylim(0, max(datac[:, plotprof[plotplot]]))
			plt.ylabel('%s (m)' % plotplot)
			plt.xlabel('Age (Ma)')
		plt.title('%s transect' % plotplot)
		plt.savefig(pgraph + '/Plot_forward_transect' + pfolder + '_' + plotplot +'.pdf')
		plt.close()
	
	print('    Saving files...')
	# save input files in the working folder
	print('       --> copying input/topo_parameters.txt to ' + pgraph)
	shutil.copy2('input/topo_parameters.txt', pgraph)
	print('       --> copying input/fault_parameters.txt to ' + pgraph)
	shutil.copy2('input/fault_parameters.txt', pgraph)
	
	print('__________________________________________')
	print(' ')	
	print(' ')
	
#### END ####