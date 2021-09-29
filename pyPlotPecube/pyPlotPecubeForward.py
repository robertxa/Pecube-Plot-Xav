######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2021 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Script to plot Pecube 4 results in forward mode
By Xavier Robert
Grenoble, 2021.07.16

xavier.robert@ird.fr

"""

###### To DO :  #######
#   - Write the code to plot the MTL pdfs
#	- Write the code to plot the time-temperature comparisons
#	- 
#   
###### End To DO #######

from __future__ import  division
# This to be sure that the result of the division of integers is a real, not an integer
# Normally not needed for Python 3

# Import modules
import os
import numpy as np
import matplotlib.pyplot as plt


#################################################
def dict_pecube():
	"""
	Definition of dictionnaries

	You may have to change it depending on your Pecube version, and on your settings...
	"""

	# agecol: respective column number of the data in the file comparison.txt
	#         needs to be changed depending on your settings in the last paragraph 
	#         of the file topo_parameters.txt
	agecol = {'alt' : 'HEIGHT',
    	      'AHE' : 'AHE',
        	  'AFT' : 'AFT',
			  'ZHE' : 'ZHE',
        	  'ZFT' : 'ZFT',
	          'KAR' : 'KAR',
    	      'BAR' : 'BAR',
        	  'MAR' : 'MAR',
	          'HBAR' : 'HAR',
    	      'FTL' : 'FT'
        	  }

	# errname: respective column number of the error on data in the data input file 
	errname = {'AHE' : 'DAHE',
			   'AFT' : 'DAFT',
			   'ZHE' : 'DZHE',
			   'ZFT' : 'DZFT',
			   'KAR' : 'DKAR',
        	   'BAR' : 'DBAR',
	           'MAR' : 'DMAR',
    	       'HBAR' : 'DHAR'
			   }

	# agename: legend of each data system         
	agename = {'AHE' : 'AHe (Ma)',
    	       'AFT' : 'AFT (Ma)',
        	   'ZHE' : 'ZHe (Ma)',
	           'ZFT' : 'ZFT (Ma)',
    	       'KAR' : 'KAr (Ma)',
        	   'BAR' : 'Biot. Ar (Ma)',
	           'MAR' : 'Musc. Ar (Ma)',
    	       'HBAR' : 'Hb Ar (Ma)',
        	   'FTL' : 'FT length (µm)',
			   'alt' : 'elevation (m)'
	           }

	# predname: legend of each predicted system         
	predname = {'AHE' : 'Predicted AHe (Ma)',
    	        'AFT' : 'Predicted AFT (Ma)',
        	    'ZHE' : 'Predicted ZHe (Ma)',
	            'ZFT' : 'Predicted ZFT (Ma)',
    	        'KAR' : 'Predicted KAr (Ma)',
        	    'BAR' : 'Predicted Biot. Ar (Ma)',
	            'MAR' : 'Predicted Musc. Ar (Ma)',
    	        'HBAR' : 'Predicted Hb Ar (Ma)',
        	    'FTL' : 'Predicted FT length (µm)',
				'alt' : 'Predicted elevation (m)'
	           }

	# colores: Colors used for the different age system
	colores = {'AHE' : 'y',
    	       'AFT' : 'r',
        	   'ZHE' : 'g',
	           'ZFT' : 'b',
    	       'KAR' : 'k',
        	   'BAR' : 'c',
	           'MAR' : 'm',
    	       'HBAR' : '0.75',
        	   'FTL' : 'y',
			   'alt' : 'y'
	           }
	
	profdict = {'Latitude'  : 'LAT', 
	            'Longitude' : 'LON', 
				'Altitude'  : 'HEIGHT', 
				'Projected' : 'PROJ'}

	return agecol, errname, agename, predname, colores, profdict


#################################################
def calc_length(XA,YA,XB,YB):       
	"""
	function to calcule the distance on a sphere
	between two points A and B given by their long/lat coordinates,   

	Read the page http://gis.stackexchange.com/questions/44064/how-to-calculate-distances-in-a-point-sequence

	"""
	# Earth radius in km
	R = 6371
	
	calc_lengthd = R * np.arccos(np.sin(YA * np.pi / 180) * np.sin(YB * np.pi / 180) + \
	                          np.cos(YA * np.pi / 180) * np.cos(YB * np.pi / 180) * \
	                          np.cos(-XA * np.pi / 180 + XB * np.pi / 180))
	return calc_lengthd


#################################################       
def bearing(XA,YA,XB,YB):       
	"""
	function to calcule the bearing on a sphere
	between two points A and B given by their long/lat coordinates,   
    """   
	
	cosdeltaAB = np.sin(YA * np.pi / 180) * np.sin(YB * np.pi / 180) + \
	                    np.cos(YA * np.pi / 180) * np.cos(YB * np.pi / 180) * \
	                    np.cos(-XA * np.pi / 180 + XB * np.pi / 180)
	bearingd = np.arccos((np.sin(YB * np.pi / 180) - cosdeltaAB * np.sin(YA * np.pi / 180)) \
	             / (np.sqrt(1 - cosdeltaAB * cosdeltaAB) * np.cos(YA * np.pi / 180)))

	return bearingd


#################################################
def project(A, B, datac):
	"""
	Function to project lat/long data along a line defined by the lat/long coordinates
	   of the two points defining the line ends.

	Args:
		A (2*1 np array of floats): beginning of the line on which to project
		B (2*1 np array of floats): end of the line on which to project
		datac (np array): data to project
	
	(c) licence CCby-nc-sa : http://creativecommons.org/licenses/by-nc-sa/4.0/ 2021

	"""

	# Earth radius in km
	R = 6371

	# Begin stepping on data
	longM = datac['LON']
	latM = datac['LAT']

	# If lat long         
	# compute on great circles
	AM = calc_length(A[0], A[1], longM, latM)
	# Calcul azimuth AB
	AzAB = bearing(A[0], A[1], B[0], B[1])
	# Calcul azimuth MM'
	AzMM = AzAB + np.pi / 2 # Verifier sens du signe
	# Calcul azimuth AM
	AzAM = bearing(A[0], A[1], longM, latM)
	# calcul Bheta = angle AM-MM from bearings of AM and MM'
	Bheta = (np.pi - AzMM + AzAM)            
	# calcul de la longueur
	coordproj = R * np.arcsin(np.sin(AM / R) * np.sin(Bheta))

	return coordproj


#################################################
def plotTransect(inputdata, profiletype, dataplot,
                 datac, inputc = None, coordproj = None, coordprojinputc = None,
				 agerange = [0, 100],
				 agecol = None, errname = None, colores = None,
				 predname = None, agename = None, profdict = None,
				 graphpath = 'Graphs', graphtitle = None):
	"""

	Args:
		inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars

		profiletype (list, optional): type of profile = ['Latitude', 'Longitude', 'Altitude', 'Projected']
						              If [], no age profile is plotted.

		dataplot (list): List of data to plot ; By default, the altitude will be plotted
						 ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
						 Rem: For the moment, MTL and TTp not implemented.

		datac (string): results of Pecube forward modeling

		inputc (array, optional): Array of data to plot.
		                          Defaults to None.

		coordproj (array, optional): Predictions projected along A-B.
		                             Defaults to None.
		
		coordprojinputc (array, optional): Observations projected along A-B.
		                                   Defaults to None.

		agerange (2*1 array of floats, optional): range of the ages to plot on the profiles
												  [min, max]. Defaults to None.

		agecol (dict, optional): respective column number of the data in the file comparison.txt.
		                         Defaults to None.
		
		errname (dict, optional): respective column number of the error on data in the data input file. 
		                          Defaults to None.
		
		colores (dict, optional): Colors used for the different age system. 
		                          Defaults to None.
				
		predname (dict, optional): legend of each predicted system. 
		                           Defaults to None.

		agename (dict, optional): legend of each data system. 
		                          Defaults to None.

		profdict (dict, optional): Columns of profile types. 
		                           Defaults to None.

		graphpath (str, optional): name of the folder where the plot will be written
						            Usually you do not have to change it. 
									Defaults to 'Graphs'.

		graphtitle (str, optional): title to write on the graph.
		                            Defaults to None.

	"""

	# Dictionnary to build axis' legends
	profname = {'Latitude'  : 'Latitude (°)', 
	            'Longitude' : 'Longitude (°)', 
				'Altitude'  : 'Elevation (m)', 
				'Projected' : 'Distance along profile (km)'}

	fig1 = plt.figure()
	for profile in profiletype:
		print(u'\tPlotting %s transect' %(profile))
		for item in dataplot:
			if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
				item = item.upper()
				if inputdata:
					# Plot input data with Error bars
					if profile != 'Projected' and profile != 'Altitude':
						plt.errorbar(inputc[profdict[profile]], inputc[agecol[item]], yerr = inputc[errname[item]], 
			    			         fmt = 'o', label = agename[item], color = colores[item])
					elif profile == 'Altitude':
						plt.errorbar(inputc[agecol[item]], inputc['HEIGHT'], xerr = inputc[errname[item]], 
			    		         fmt = 'o', label = agename[item], color = colores[item])
					else:
						plt.errorbar(coordprojinputc, inputc[agecol[item]], yerr = inputc[errname[item]], 
			    		         fmt = 'o', label = agename[item], color = colores[item])
				else:
					# Plot input data with NO error bars
					if profile != 'Projected' and profile != 'Altitude':
						plt.plot(datac[profdict[profile]], datac[agecol[item]+'OBS'], 
			    			     marker = 'o', linestyle = 'None', 
								 label = agename[item], color = colores[item])
					elif profile == 'Altitude':
						plt.plot(datac[agecol[item]+'OBS'], datac['HEIGHTOBS'], 
			    			     marker = 'o', linestyle = 'None', 
								 label = agename[item], color = colores[item])
					else:
						plt.plot(coordproj, datac[agecol[item]+'OBS'], 
			    			     marker = 'o', linestyle = 'None', 
								 label = agename[item], color = colores[item])
				# Plot predictions
				if profile != 'Projected' and profile != 'Altitude':
					plt.plot(datac[profdict[profile]], datac[agecol[item]+'PRED'], 
		    			     marker = 's', linestyle = 'None', 
							 label = predname[item], color = colores[item], alpha = 0.3)
				elif profile == 'Altitude':
					plt.plot(datac[agecol[item]+'PRED'], datac['HEIGHTPRED'], 
		    	    		 marker = 's', linestyle = 'None',
							 label = predname[item], color = colores[item], alpha = 0.3)	
				else:
					plt.plot(coordproj, datac[agecol[item]+'PRED'], 
		    			     marker = 's', linestyle = 'None', 
							 label = predname[item], color = colores[item], alpha = 0.3)

		# Remove -9999 values from the graph (= no data values)	
		if profile != 'Altitude':
			plt.ylim(bottom = agerange[0], top = agerange[1])
		else:
			plt.ylim(bottom = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']), 0),
		 		        top = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED'])))
			plt.xlim(left = agerange[0], right = agerange[1])

		# set legend
		plt.legend(loc = 'best')
		#plt.xlabel(u'Longitude (°)')
		if profile != 'Altitude':
			plt.xlabel(profname[profile])
			plt.ylabel(u'Age (Ma)')
		else:
			plt.ylabel(profname[profile])
			plt.xlabel(u'Age (Ma)')
		plt.title(graphtitle)

		plt.savefig(graphpath + '/profile/' + graphtitle +'_' + profile + '.pdf')
		fig1.clear()

	return


#################################################
def plotComparisons(inputdata, dataplot,
                    datac, inputc = None,
					agecol = None, errname = None, colores = None, agename = None,
					graphpath = 'Graphs', graphtitle = None):
	"""
	Function to plot the comparison between data and predictions

	Args:
		inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars

		dataplot (list): List of data to plot ; By default, the altitude will be plotted
						 ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
						 Rem: For the moment, MTL and TTp not implemented.
						 Defaults to ['AHe','AFT'].

		datac (string): results of Pecube forward modeling

		inputc (array, optional): Array of data to plot.
		                          Defaults to None.

		agecol (dict, optional): respective column number of the data in the file comparison.txt.
		                         Defaults to None.
		
		errname (dict, optional): respective column number of the error on data in the data input file. 
		                          Defaults to None.
		
		colores (dict, optional): Colors used for the different age system. 
		                          Defaults to None.
		
		agename (dict, optional): legend of each data system. 
		                          Defaults to None.

		graphpath (str, optional): name of the folder where the plot will be written
						           Usually you do not have to change it. 
								   Defaults to 'Graphs'.

		graphtitle (str, optional): title to write on the graph. 
		                            Defaults to None.

	"""

	for item in dataplot:
		if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
			if item == 'Altitude':
				inputdata = None
				print(u'\tPlotting %s comparison' %(str(item)))
				item = 'alt'
			else:
				item = item.upper()
				print(u'\tPlotting %s age comparison' %(str(item)))
			fig2 = plt.figure()
			plt.gca().set_aspect('equal')
			# Plot 1:1 line
			plt.plot([0,8000], [0,8000], marker = 'None', linestyle = '-', 
					 label = '1:1 line', color = 'lightgrey', alpha = 1)
			# Plot the age comparison
			# And add error bars on data
			if inputdata:
				plt.errorbar(datac[agecol[item]+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)], 
							 datac[agecol[item]+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)],
							 xerr = inputc[errname[item]][np.logical_not(np.isnan(inputc[errname[item]]))],
			    	    	 fmt = 'o', label = 'Sample', color = colores[item])
			else:
				plt.plot(datac[agecol[item]+'OBS'], datac[agecol[item]+'PRED'], 
					 		 marker = 'o', linestyle = 'None',
							 label = 'Sample', 
							 color = colores[item])

			# Write the 1:1 text over the line
			# We may need to revise the way to find the x/y
			plt.text((max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2,
					 (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2 - 
					      (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/20,
				 	'1:1', rotation = 45, color = 'lightgrey', alpha = 1, ha = 'center', va = 'center')

			plt.ylim(bottom = min(min(abs(datac[agecol[item]+'OBS'])), 
						 		  min(abs(datac[agecol[item]+'PRED']))), 
	    	    	 top = max(max(datac[agecol[item]+'OBS']), 
				 			   max(datac[agecol[item]+'PRED'])))
			#plt.ylim(bottom = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']), 0), 
	        # top = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED'])))
			plt.xlim(left = min(min(abs(datac[agecol[item]+'OBS'])), 
							min(abs(datac[agecol[item]+'PRED']))), 
				 right = max(max(datac[agecol[item]+'OBS']), 
				 			 max(datac[agecol[item]+'PRED'])))
			#plt.xlim(left = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']), 0), 
			# right = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED'])))

			#plt.xlabel(u'%s Observed (Ma)' %(str(item)))
			#plt.ylabel(u'%s Predicted (Ma)' %(str(item)))
			plt.xlabel(u'Observed %s' %(str(agename[item])))
			plt.ylabel(u'Predicted %s' %(str(agename[item])))
			plt.legend(loc = 'best')
			plt.title(graphtitle)

			plt.savefig(graphpath + '/Compare/' + graphtitle +'_Compare_' +str(item) + '.pdf')
			fig2.clear()

	return

#################################################
def plotTTp(inputdataPTt, outputdataPTt, graphpath = 'Graphs', graphtitle = None):
	"""
	Function to plot Time-Temperature paths

	Args:
		inputdataPTt (str): path and name of the input TTp file
		outputdataPTt (str): path and name of the output TTpfile.
		                     Generally, this is 'output/CompareTT.csv'
		graphpath (str, optional): Path to the graphs folder.
		                           Defaults to 'Graphs'.
		graphtitle (str, optional): Title of the graph.
		                            Defaults to None.
	"""

	# Read the inputdataPTt (.csv) file
	inputPTt = np.genfromtxt(fname = inputdataPTt, delimiter = ',', names = True, case_sensitive = 'upper')
	# Read the Predicion PT-t (.csv) file
	# PROBLEM TO READ THE OUTPUT FILE COMPARE-TT.CSV if not worked in a text editor ?
	# Clean it:
	with open(outputdataPTt, 'r') as f1:
		with open(outputdataPTt + 'touched', 'w') as f2:
			lines = f1.readlines()
			changes = False
			for k in range (0, len(lines)-1):
				if lines[k] != '':
					# Remove the ',' at the eand of the line if there is one
					if lines[k][-2] == ',' :
						if lines[k][-1] != ',' :
							lines[k] = lines[k][0:-2]+lines[k][-1]
							changes = True
					# merge line k and line k+1 if lat/long not followed by data
					if len(lines[k].split(',')) == 2 and lines[k+1].split(',')[0] == '' and lines[k+1].split(',')[1] == '':
						lines[k] = lines[k].rstrip('\n') + lines[k+1][1:]
						lines[k+1] = ''
						changes = True
			if changes:
				# Write the lines a new file, but not the empty lines
				for line in lines:
					if line.strip('\n') != '':
						f2.write(line)
				# update the name of the output PTt file with the corrected file
				outputdataPTt = outputdataPTt + 'touched'
				print('\t\033[91mWarning:\033[00m Predicted input file modified to %s to be plotted' %(outputdataPTt))

	# Read the cleaned Predicion PT-t (.csv) file
	outputPTt = np.genfromtxt(fname = outputdataPTt, delimiter = ',', names = True, case_sensitive = 'upper')

	# Find number of samples in the input file
	nPTtsamples = max(inputPTt[np.isnan(inputPTt['LAT']) == False].shape[0],
	    				inputPTt[np.isnan(inputPTt['LON']) == False].shape[0],
						inputPTt[np.isnan(inputPTt['HEIGHT']) == False].shape[0],
						inputPTt[np.isnan(inputPTt['SAMPLE']) == False].shape[0])
	# Find indexes of sample names/beginning
	indexPTt = np.where(np.isnan(inputPTt['LAT']) == False)[0]
	
	if nPTtsamples > 0:
		fig3 = plt.figure()
		# For each sample with PT-t, 
		for k in range (0, nPTtsamples):
			#	Extract the PTt path from the file
			if k == (nPTtsamples-1):
				ppt = inputPTt[['TIMEH', 'TEMPH', 'DTEMPH']][indexPTt[k]:]
				pptOut = outputPTt[['TIME', 'TEMP', 'TEMPPRED']][indexPTt[k]:]
			else:
				ppt = inputPTt[['TIMEH', 'TEMPH', 'DTEMPH']][indexPTt[k]:indexPTt[k + 1]]
				pptOut = outputPTt[['TIME', 'TEMP', 'TEMPPRED']][indexPTt[k]:indexPTt[k + 1]]
			# Plot the enveloppe deduced from the error bars on T
			plt.fill_between(x = ppt['TIMEH'], 
    		                 y1 = ppt['TEMPH'] + ppt['DTEMPH'],
    		                 y2 = ppt['TEMPH'] - ppt['DTEMPH'],
    		                 alpha=0.20, 
    		                 color='lightblue', 
    		                 interpolate=True,
    		                 label = 'Acceptable paths')
			# plot the max enveloppe
			plt.plot(ppt['TIMEH'],
			         ppt['TEMPH'] + ppt['DTEMPH'],
					 marker = 'None', linestyle = '-', 
				 	 color = 'lightblue', alpha = 1,
					 label = 'Min-Max')
			# Plot the min enveloppe
			plt.plot(ppt['TIMEH'],
			         ppt['TEMPH'] - ppt['DTEMPH'],
					 marker = 'None', linestyle = '-', 
				 	 color = 'lightblue', alpha = 1)
			# 	Plot the mean PT-t path
			plt.plot(ppt['TIMEH'],
			         ppt['TEMPH'],
					 marker = 'None', linestyle = '-', 
				 	 color = 'blue', alpha = 1,
					 label = 'Mean')
			#	Plot the prediction
			plt.plot(pptOut['TIME'],
			         pptOut['TEMPPRED'],
					 marker = 'None', linestyle = '-', 
				 	 color = 'green', alpha = 1,
					 label = 'Predictions')

			#	Save the graph
			plt.xlabel(u'Time before present (Ma)')
			plt.ylabel(u'Temperature (Ma)')
			plt.legend(loc = 'best')
			plt.title(graphtitle)
			# Invert x- and y-axis
			plt.axis([max(ppt['TIMEH']), min(ppt['TIMEH']),
					  max(ppt['TEMPH'] + ppt['DTEMPH']), min(ppt['TEMPH'] - ppt['DTEMPH'])])
			plt.savefig(graphpath + '/TTpaths/' + graphtitle +'_ttpath_sample' + str(k+1) + '.pdf')
			fig3.clear()
	else:
		print('\033[91mWarning:\033[00m No PT-t paths to plot...')	

	return


#################################################
def plotMTL():

	print('\t\033[91mWarning:\033[00m MTL not implemented for now...')	
	
	# TO DO

	# Read input file

	# Read Output file

	# Find number of sample with MTL
	# For each sample with MTL,
	#	Build input histogram
	#	Build output histogram
	#	Save the graph

	return

#################################################
################# Main code #####################
def PlotPecubeForward(datafnme, inputdata, inputdataPTt = None, outputdataPTt = None,
					  dataplot = ['AHe','AFT'],
					  graphpath = 'Graphs', graphtitle = None, agerange = None,
					  profiletype = [], A = None, B = None,
					  size_x = 15, size_y = 15,
					  agename = None, predname = None, colores = None):
	"""[summary]

	Args:
		datafnme (string): results of Pecube forward modeling

		inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars

		inputdataPTt (str, optional): Path and file name of the input PTt data file 
									  Need to be given if 'TTp' in dataplot
        							  Default = None 

    	outputdataPTt (str, optional): Path and file name of the output PTt prediction/comparison file
									   Usually, this is 'output/CompareTT.csv'
									   Need to be given if 'TTp' in dataplot
        							   Default = None

		dataplot (list, optional): List of data to plot ; By default, the altitude will be plotted
				        		    Do not forget the simple quotes !!!
						   			['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
				 					Rem: For the moment, MTL and TTp not implemented.
									Defaults to ['AHe','AFT'].

		graphpath (str, optional): name of the folder where the plot will be written
						            Usually you do not have to change it. Defaults to 'Graphs'.

		graphtitle (str, optional): title to write on the graph. Defaults to None.

		agerange (2*1 array of floats, optional): range of the ages to plot on the profiles
												  [min, max]. Defaults to None.

		profiletype (list, optional): type of profile = ['Latitude', 'Longitude', 'Altitude', 'Projected']
						              If [], no age profile is plotted. Defaults to [].

		size_x      : Font size for x-axis. Defaults to 15.
        
		size_y      : Font size for y-axis. Defaults to 15.

		A, B (floats, optional): If need of a projected transect, define the line along which we will project
						 		  With the coordinate of the point A and B defining espectivelly
						         the begining and the end of the transect, in lat-long/WGS84. Defaults to None.
		
		agename (dict, optional): legend of each data system         
								  If None, this is set to
									{'AHE' : 'AHe (Ma)',
									'AFT' : 'AFT (Ma)',
									'ZHE' : 'ZHe (Ma)',
									'ZFT' : 'ZFT (Ma)',
									'KAR' : 'KAr (Ma)',
									'BAR' : 'Biot. Ar (Ma)',
									'MAR' : 'Musc. Ar (Ma)',
									'HBAR' : 'Hb Ar (Ma)',
									'FTL' : 'FT length (µm)'}
								  Default = None. 

		predname (dict, optional): legend of each predicted system.
								   If None, this is set to
									{'AHE' : 'Predicted AHe (Ma)',
									'AFT' : 'Predicted AFT (Ma)',
									'ZHE' : 'Predicted ZHe (Ma)',
									'ZFT' : 'Predicted ZFT (Ma)',
									'KAR' : 'Predicted KAr (Ma)',
									'BAR' : 'Predicted Biot. Ar (Ma)',
									'MAR' : 'Predicted Musc. Ar (Ma)',
									'HBAR' : 'Predicted Hb Ar (Ma)',
									'FTL' : 'Predicted FT length (µm)'}         
								   Default = None. 
								   
		colores (dict, optional): Colors used for the different age system
								  If None, this is set to
									{'AHE' : 'y',
									'AFT' : 'r',
									'ZHE' : 'g',
									'ZFT' : 'b',
									'KAR' : 'k',
									'BAR' : 'c',
									'MAR' : 'm',
									'HBAR' : '0.75',
									'FTL' : 'y',
									'alt' : 'y'}
								  Default = None.

	Raises:
		NameError: Problem with input files !
	"""

	print('###########################################################################\n')
	print(u'\tPlotting Pecube V4+ forward modelling results...\n')
	print('\t\t\xa9Xavier Robert - IRD-ISTerre\n')
	print('###########################################################################\n')

	# Check if the input files exist,
	if not os.path.isfile(datafnme):
		raise NameError(u'\033[91mERROR:\033[00m F** input file %s does not exist' % datafnme)
	if not os.path.isfile(inputdata):
		print(u'\n \033[91mWarning:\033[00m No %s file, I am skipping the plot of the error bars...\n' % inputdata)
		inputdata = None
	# Check if the Graphs/ folder exists, if not create it
	if os.path.exists(graphpath) == False:
		print(u'Output folder %s/ does not exist...' %(str(graphpath)))
		print(u'I am creating it...')
		os.mkdir(graphpath)
	else:
		print(u'\033[91mWarning:\033[00m Folder %s/ already exists, I will write in it and erase previous pdf files' %(str(graphpath)))
	
	# Read data files
	datac = np.genfromtxt(fname = datafnme, delimiter = ',', names = True, case_sensitive = 'upper')
	if inputdata: inputc = np.genfromtxt(fname = inputdata, delimiter = ',', names = True, case_sensitive = 'upper')

	# Read dictionnaries
	if not agename:
		if not predname:
			if not colores:
				agecol, errname, agename, predname, colores, profdict =  dict_pecube()
			else:
				agecol, errname, agename, predname, junk, profdict =  dict_pecube()
		else:
			agecol, errname, agename, junk1, junk, profdict =  dict_pecube()
	else:
		agecol, errname, junk2, junk1, junk, profdict =  dict_pecube()

	# Plot profiles
	if profiletype:
		if os.path.exists(graphpath + '/Profile') == False:
			os.mkdir(graphpath + '/Profile')
		coordproj = None
		coordprojinputc = None
		if 'Projected' in profiletype and A and B:
			# plot Projected transect
			# Compute projected coordinates
			coordproj = project(A, B, datac)
			if inputdata: coordprojinputc = project(A, B, inputc)
		else:
			print(u'\t\033[91mWarning:\033[00m No points given to project along a profile, I am skipping it...')
			profiletype.remove('Projected')

		plotTransect(inputdata, profiletype, dataplot,
                 datac, inputc = inputc, coordproj = coordproj, coordprojinputc = coordprojinputc,
				 agerange = agerange,
				 agecol = agecol, errname = errname, colores = colores,
				 predname = predname, agename = agename, profdict = profdict,
				 graphpath = graphpath, graphtitle = graphtitle)

	########
	# Plot Altitude comparison
	if os.path.exists(graphpath + '/Compare') == False:
			os.mkdir(graphpath + '/Compare')
	# Call plot function for altitude
	#print(u'\tPlotting altitude comparison')
	plotComparisons(inputdata = inputdata, dataplot = ['Altitude'],
                    datac = datac, inputc = inputc, 
					agecol = agecol, errname = errname, colores = colores,
					agename = agename,
					graphpath = graphpath, graphtitle = graphtitle)
	# Do also Age comparisons
	# Call plot function for dataplot items
	plotComparisons(inputdata = inputdata, dataplot = dataplot,
                    datac = datac, inputc = inputc, 
					agecol = agecol, errname = errname, colores = colores,
					agename = agename,
					graphpath = graphpath, graphtitle = graphtitle)

	# Check if input and output PTt path files are present
	if (inputdataPTt and outputdataPTt) and (not os.path.isfile(inputdataPTt) or not os.path.isfile(outputdataPTt)):
		print(u'\n \033[91mWarning:\033[00m No %s and/or %s file, I am skipping the PTt plot...\n' % (inputdataPTt, outputdataPTt))
		inputdataPTt = None
		outputdataPTt = None
	if inputdataPTt and outputdataPTt and \
	   ('TTp' in dataplot or 'Ttp' in dataplot or 'ttp' in dataplot or 'TTP' in dataplot or
	    'tTP' in dataplot or 'ttP' in dataplot or 'TtP' in dataplot or 'tTp' in dataplot):
		print('\tPlotting PTt paths')
		if os.path.exists(graphpath + '/TTpaths') == False:
			os.mkdir(graphpath + '/TTpaths')
		# Call the plot PTt function
		plotTTp(inputdataPTt = inputdataPTt, 
		        outputdataPTt =outputdataPTt,
				graphpath = graphpath, graphtitle = graphtitle)	
	else:
		print('\033[91mWarning:\033[00m No PT-t paths to plot...')

	# TO DO !!!!
	if ('MTL' in dataplot or 'MTl' in dataplot or 'mtl' in dataplot or 'mtL' in dataplot or
	    'mTL' in dataplot or 'MtL' in dataplot or 'mtl' in dataplot or 'mTl' in dataplot):
		print('\tPlotting MTL')
		if os.path.exists(graphpath + '/MTL') == False:
			os.mkdir(graphpath + '/MTL')
		plotMTL()
	# END - TO DO !!!!

	print('\n###########################################################################\n')
	print('\t\t\t\tEND\n')
	print('###########################################################################\n')
	### End


#################################################
#################################################
if __name__ == "__main__":
	# Define data to plot
	# dataplot: List of data to plot ; By default, the altitude will be plotted
	#           Do not forget the simple quotes !!!
	#           ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
	# 			Rem: For the moment, MTL and TTp not implemented
	dataplot = ['AHe','AFT']

	# graphpath: name of the folder where the plot will be written
	#            Usually you do not have to change it
	graphpath = '../Tests/Forward/Graphs'

	# Data to plot
	#	datafnme:  results of Pecube forward modeling
	#	inputdata: input data declared in Pecube.in; This is used to plot the errorbars
	datafnme = '../Tests/Forward/Data/CompareAGE.csv'
	inputdata = '../Tests/Forward/Data/Trujillo.csv'
	inputdataPTt = '../Tests/Forward/Data/TrujilloPTt.csv'
	outputdataPTt = '../Tests/Forward/Data/TimeTemperaturePaths.csv'

	# graphtitle: title to write on the graph
	graphtitle = 'Trujillo transect'

	#agerange = range of the ages to plot on the profiles
	#			[min, max]
	agerange = [0, 100]

	# profiletype: type of profile = ['Latitude', 'Longitude', 'Altitude', 'Projected']
	#              if [], no age profile is plotted
	#              For the moment, NO there is no projected profile; This is one thing to add ?
	#profiletype = ['Longitude']
	#profiletype = ['Projected']
	profiletype = ['Latitude', 'Longitude', 'Altitude', 'Projected']
	#profiletype = ['Latitude', 'Altitude']

	# A, B =  If need of a projected transect, define the line along which we will project
	# 		  With the coordinate of the point A and B defining espectivelly
	#         the begining and the end of the transect, in lat-long/WGS84
	A = [-79.1, -8.21]
	B = [-78.41, -7.83]

	# Font size of x and y axis
	size_x = 15
	size_y = 15

	# end define the data and parameters
	#######################################

	PlotPecubeForward(dataplot = dataplot,
					  graphpath = graphpath,
					  datafnme = datafnme,
					  inputdata = inputdata,
					  inputdataPTt = inputdataPTt,
					  outputdataPTt = outputdataPTt,
					  graphtitle = graphtitle,
					  agerange = agerange,
					  profiletype = profiletype,
					  size_x = size_x, size_y = size_y,
					  A = A, B = B)
