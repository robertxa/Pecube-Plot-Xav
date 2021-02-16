######!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to plot Pecube results in NA mode
By Xavier Robert
Grenoble, 2015.06.11

USAGE :
  1- Copy this file in the main Pecube/ folder
  2- After an inversion, run python Plot-Pecube-NA.py
  3- The plots will be in the folder NA/Graphs/prefix 
  4- The script copy Pecube's output files in the folder NA/Graphs/prefix 

INPUTS:
The inputs are in the script file, in the "# Define data to analysis" section. 
The different arguments are described.

xavier.robert@ujf-grenoble.fr

(c) licence CCby-nc : http://creativecommons.org/licenses/by-nc/3.0/ 2015

"""

from __future__ import  division
# This to be sure that the result of the division of integers is a real, not an integer

# Import modules
import sys
import os
import traceback
import copy
import numpy as np
import scipy as sp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import shutil


#######################################
# Define data to analysis
# param: list of the name of the different parameters
#        given in the same order than the columns in NA_results.txt
#        or described in na.sum
param = [u'Thermal diffusivity (km2/Myr)', 
         u'Basal Temperature (°C)', 
         u'Heat production (°C/My)', 
         u'V_MBT (km/My)', 
         u'V_MFT (km/My)']
# dataplot: List of couple of parameters to plot
#           [(1, 2), (2, 3), (3, 4), ...]
#           If [], do not plot couple of parameters
#dataplot = [(4,1), (4,2), (4,3), (4,5), (1,2), (2,3), (1,3)]
dataplot = [(4,1), (4,2), (4,3), (4,5)]
# suffix: suffix to use to save the graph and backup files
suffix = u'Xav'
# nasum: name or path (inside the folder NA/) of the na.sum file
#nasum = 'na.sum'
# ipecube: False if we are working outside the classic Pecube directory
#          True if we are working inside the classic Pecube directory
ipecube = False
# NApath: Path where is stored the NA_results.txt file
#         Should be u'' if ipecube = True and if Na_results.txt is in NA/
NApath = u'Test_Fin3'
# naresults: name or path (inside the folder NA/) of the NA_results.txt file
naresults = u'NA_results.txt'
# graphpath: name of the folder where the plot will be written
#            Usually you do not have to change it
graphpath = u'Graphs'
### Options to plot grahs
# sizex, sizey: size of the plot in cm;
#               set them to None if you want automatic scaling
sizex = None
sizey = None
#sizex = 16
#sizey = 13.5
# symbolsize: size of the symbols
symbolsize = 15
#sizefont: size of the fonts for the graphs. 
#          Set it to None if you want to use default size
sizefont = None
#sizefont = 16
# end define the data and parameters
#######################################


######## Beginning of functions definition ###########
def plotindiv(datac, nbplotx, nbploty, nbplot, i, param, symbolsize, sizefont = None):
	"""
	Function to plot the data individually in function of misfit

	      Written by Xavier Robert 2015.06.11
	
	INPUT:
	  datac = array to plot
	  nbplotx = number of subplots in the x direction 
	  nbploty = number of subplots in the y direction
	  i = number of the subplot
	  param = list of the name of the parameters that will be plotted
	  sizefont = size of the font
	  symbolsize = size of the symbols
	OUTPUT:
	  graph
	USAGE:
	  plotindiv(datac, nbplotx, nbploty, nbplot, i, param)
	"""
	
	ax = plt.subplot(nbploty, nbplotx, i+1)
	plt.plot(datac[:, i+1], datac[:, 0], 
	         'o', 
	         color = 'r',
	         markersize = symbolsize)
	if sizefont != None:
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
		            ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(sizefont)
	plt.xlabel(param)
	plt.ylabel('Misfit')
	plt.grid(True)
	
	return
	# End of the function


def plotcouple(datac, dataplot, couple, nbplotx, nbploty, param):
	"""
	Function to plot the couple data (with a least 2 parameters)
	
	      Written by Xavier Robert 2015.06.11
	
	INPUT:
	  datac = array to plot
	  dataplot = list of couple to plot (ex : [(1,2),(3,4)])
	  couple = couple of parameters to plot (ex : (1,2))
	  nbplotx = number of subplots in the x direction 
	  nbploty = number of subplots in the y direction
	  param = list of the name of the parameters that will be plotted
	OUTPUT:
	  graph
	USAGE:
	  plotcouple(datac, dataplot, couple, nbplotx, nbploty, param)
	"""
	
	# Find the subplot number : dataplot.index(couple)
	plt.subplot(nbploty, nbplotx, dataplot.index(couple) + 1)
	# plot all the data
	#cm = plt.cm.cool
	cm = plt.cm.get_cmap('RdYlBu')
	sc = plt.scatter(datac[:, couple[0]],datac[:, couple[1]], 
	                 c = datac[:, 0],
	                 s = 50,
	                 cmap=cm,)
	
	# plot the lowest misfit with a red star
	plt.plot(datac[np.argmin(datac[:,0]), couple[0]], 
	         datac[np.argmin(datac[:,0]), couple[1]], 
	         marker='*', 
	         c='r',
	         markersize=20,
	         markeredgewidth=1,
	         label=u'lowest misfit')
	plt.colorbar(sc, label = u'Misfit')
	plt.xlabel(param[couple[0] - 1])
	plt.ylabel(param[couple[1] - 1])
	
	# find the min and max values of axes
	# Initiate
	minx = np.zeros(2)
	maxx = np.zeros(2)  
	# Find the range and precision of the parameters
	# i.e. if it is between 0.001 and 0.01 or between 10 and 1000 for instance
	for i in range (0,2):
		rangex = max(datac[:, couple[i]]) - min(datac[:, couple[i]])
		dec = 0
		if rangex < 0.0001 and rangex > 0.00001:
			prec = 0.00001
			dec = 5		
		elif rangex < 0.001 and rangex > 0.0001:
			prec = 0.0001
			dec = 4		
		elif rangex < 0.01 and rangex > 0.001:
			prec = 0.001
			dec = 3	
		elif rangex < 0.1 and rangex > 0.01:
			prec = 0.01
			dec = 2		
		elif rangex < 1 and rangex > 0.1:
			prec = 0.1
			dec = 1		
		elif rangex < 10 and rangex > 1:
			prec = 1
		elif rangex < 100 and rangex > 10:
			prec = 10
		elif rangex < 1000 and rangex > 100:
			prec = 100
		elif rangex < 10000 and rangex > 1000:
			prec = 1000
		elif rangex < 100000 and rangex > 10000:
			prec = 10000	
		else:
			# If the range of the parameter is out of defined bounds, 
			#   raise an Error to force the user to program is own bounds
			try:
				raise Exception("foo")
			except:
				for frame in traceback.extract_tb(sys.exc_info()[2]):
					fname,lineno,fn,text = frame
					raise ParamError('ERROR: The range of the parameter %s is out of bounds. \n'
			    	                'You need to change the code before line %d \n' 
			    	                'Add the range you need following the existing syntax'
			    	                % (param[couple[0] - 1 + i], lineno))
		# Calcul the min 
		minx[i] = np.around(min(datac[:, couple[i]]), decimals=dec)
		# Check if the computed min is lower than the real min
		if minx[i] > min(datac[:, couple[i]]):
			# if not, substract one unit (calculated some lines higher) to the computed min
			minx[i] = minx[i] - prec
		# Do the same with the maximum
		maxx[i] = np.around(max(datac[:, couple[i]]), decimals=dec)
		if maxx[i] < max(datac[:, couple[i]]):
			maxx[i] = maxx[i] + prec
	
	# plot contours
	plot_coutour(datac[:, couple[0]], datac[:, couple[1]], datac[:, 0], minx, maxx, levels = [0.7])

	plt.xlim(minx[0], maxx[0])
	plt.ylim(minx[1], maxx[1])	
	
	return
	# End of the function


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
				raise ParamError('ERROR: The number of parameters is more than 12 (in fatc %i). \n'
			         'You need to change the code before line %d' % (nbplot, lineno))
			         
	return nbplotx, nbploty
	# End of the function


def cm2inch(value):
	"""
	function to transform the size from cm to inches
	
	INPUT:
	  - Value = real in cm
	OUTPUT:
	  - value in inches
	"""
	
	inch = 2.54
	return value/inch
	
	# to work with *tupl
	#if isinstance(tupl[0], tuple):
	#	return tuple(i/inch for i in tupl[0])
	#else:
	#	return tuple(i/inch for i in tupl)


def check_files(NApath = '', naresults = u'NA_results.txt', graphpath = u'Graphs', suffix = '', ipecube = True):
	"""
	Function to check where are the files
	
	INPUTS:
	  - NApath = 
	  - naresults =
	  - graphpath =
	  - suffix =
	  - ipecube =
	OUTPUTS:
	  - pdata =
	  - graphpath =
	USAGE:
	  
	"""
	
	if NApath != '' and NApath[-1] != u'/':
		NApath = NApath + u'/'
		
	if ipecube:	
		if os.path.exists(u'NA/') == False:
			raise NameError('   ERROR : Folder NA/ does not exist')
		else:
			print('   Folder NA/ already exists, I will use it')
			#if os.path.isfile(u'NA/' + NApath + naresults) == True:
			pdata = u'NA/' + NApath + naresults
			#elif os.path.isfile(u'NA/' + naresults) == True:
			#	pdata = u'NA/' + naresults
					
		if os.path.exists(u'NA/' + graphpath) == False:
			print('   Output folder does not exist...')
			print('   I am creating it...')
			print(' ')
			os.mkdir(u'NA/' + graphpath)
			graphpath = u'NA/' + graphpath
		else:
			print('    Folder NA/' + graphpath + '/ already exists, I will use it and potentially erase previous pdf file')
			graphpath = u'NA/' + graphpath
		# Check if the Graphs/suffix folder exists
		if suffix != '':
			if os.path.exists(graphpath + '/' + suffix) == False:
				print('   Output folder does not exist...')
				print('   I am creating it...')
				print(' ')
				os.mkdir(graphpath + u'/' + suffix)
			else:
				print('   Folder ' + graphpath + '/' + suffix + '/ already exists, '
				      'I will write in it and erase previous pdf file')	
	else:
		pdata = NApath + naresults
		if os.path.exists(NApath + graphpath) == False:
			print('   Output folder does not exist...')
			print('   I am creating it...')
			print(' ')
			os.mkdir(NApath + graphpath)
			graphpath = NApath + graphpath
		else:
			print('    Folder ' + NApath + graphpath + '/ already exists, I will use it and potentially erase previous pdf file')
			graphpath = NApath + graphpath
		# Check if the Graphs/suffix folder exists
		#if suffix != '':
		#	if os.path.exists(graphpath + '/' + suffix) == False:
		#		print('   Output folder does not exist...')
		#		print('   I am creating it...')
		#		print(' ')
		#		os.mkdir(graphpath + u'/' + suffix)
		#	else:
		#		print('   Folder ' + graphpath + '/' + suffix + '/ already exists, '
		#		      'I will write in it and erase previous pdf file')
	
	# Check if file NA_results.txt exists, if not, raise an error
	if os.path.isfile(pdata) == False and os.access(pdata, os.R_OK) == False :
		raise NameError('ERROR : File {FileNa} does not exist'.format(FileNa=str(pdata)))
	
	return pdata, graphpath

def plot_coutour(x, y, z, minx, maxx, levels = [1]):
	"""
	
	"""
	# define grid.
	xi = np.linspace(minx[0], maxx[0], 100)
	yi = np.linspace(minx[1], maxx[1], 100)
	# grid the data.
	zi = mlab.griddata(x, y, z, xi, yi)
	# contour the gridded data, plotting dots at the randomly spaced data points.
	#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
	CS = plt.contour(xi, yi, zi, linewidths=2, colors='k', levels = levels)
	#plt.clabel(CS, inline=1, fontsize=10)
	#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
	#plt.colorbar() # draw colorbar
	
	return

######## END of the functions definition ###########


######### Main code ###############
if __name__ == "__main__":
	#### Read the data
	#print(' ')
	print('___________________________________________')
	print(' ')
	print('Plotting data for a Pecube inverse model...')
	print('___________________________________________')
	
	
	pdata, graphpath = check_files(NApath, naresults, graphpath, suffix, ipecube)
	
	if suffix != '': suffix = u'_' + suffix
	
	print('    Reading data...')
	# load data
	datac = np.loadtxt(pdata)
	
	print('    Plotting individual data')
	# Determine the number of plot (max 8)
	nbplot = datac.shape[1] - 1   # -1 to remove the first column
	print('         %i parameter(s) to plot' % nbplot)
	nbplotx, nbploty = find_nbplot(nbplot)
	
	# Size of the plot
	if sizex != None and sizey != None:
		plt.figure(1, figsize = (cm2inch(sizex), cm2inch(sizey)))
	else:
		plt.figure(1)
	
	if len(param) != nbplot:
		print('         number of axis legend in "param": %s' % nbplot)
		raise AxisError('Error: number of axis legend should be identical to the number of subplot/parameters')
	
	plt.subplots_adjust(wspace = (nbplot + 2)/10., hspace = (nbplot + 2)/10.)
	for i in range (0,nbplot):
		# call plotindiv
		plotindiv(datac, nbplotx, nbploty, nbplot, i, param[i], symbolsize, sizefont)
	
	plt.savefig(graphpath + u'/NAgraph_indiv' + suffix + u'.pdf')
	
	# Plotting scatter plots
	if len(dataplot) != 0:
		print('    Plotting data by couple')
		nbplot = len(dataplot)
		print('         %i couple(s) to plot' % nbplot)	
		nbplotx, nbploty = find_nbplot(nbplot)
		
		plt.figure(2)
		plt.subplots_adjust(wspace = (nbplot + 2)/10., hspace = (nbplot + 2)/10.)
		for couple in dataplot:
			# call plotcouple
			plotcouple(datac, dataplot, couple, nbplotx, nbploty, param)
		plt.savefig(graphpath + u'/NAgraph_couple' + suffix + u'.pdf')
	else:
		print('    No couple to plot')


	# If we are working in the classic Pecube directory, save Pecube's input files in the graph folder	
	if naresults == 'NA_results.txt' and ipecube == True:
		if suffix != '': suffix = suffix[1:]
		print('    Saving files...')
		print('       --> copying input/topo_parameters.txt to ' + graphpath + u'/' + suffix)
		shutil.copy2(u'input/topo_parameters.txt', graphpath + u'/' + suffix)
		print('       --> copying input/fault_parameters.txt to ' + graphpath + u'/' + suffix)
		shutil.copy2(u'input/fault_parameters.txt', graphpath + u'/' + suffix)
		print('       --> copying NA/'+ nasum + ' in NA/'+ graphpath + '/' + suffix)
		shutil.copy2(u'NA/' + nasum, u'NA/'+ graphpath + u'/' + suffix)
		print('       --> copying NA/'+ naresults + ' in NA/'+ graphpath + '/' + suffix)
		shutil.copy2(u'NA/' + naresults, u'NA/'+ graphpath + u'/' + suffix)
		print('       --> copying NA/na.nad in NA/'+ graphpath + '/' + suffix)
		shutil.copy2(u'NA/na.nad', u'NA/'+ graphpath + u'/' + suffix)
	else:
		print('    No file to save')
	print('__________________________________________')
	print(' ')	
	print(' ')
	
#### END ####