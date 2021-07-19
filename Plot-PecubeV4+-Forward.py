######!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to plot Pecube 4 results in forward mode
By Xavier Robert
Grenoble, 2021.07.16

USAGE :
  1- Copy this file in the main Pecube/RUNXX folder
  2- Set (edit) which system you want to plot in the "# Define data to analysis" section (see below)
  3 - Run in the terminal: $ python Plot-Pecube-Forward.py after the run of ./bin/Pecube
  4- The plots will be in the folder RUNXX/Graphs/ where RUNXX/ is the working folder

INPUTS:
The inputs are in the script file, in the "# Define data to analysis" section. 
The different arguments are described.

xavier.robert@ird.fr

(c) licence CCby-nc-sa : http://creativecommons.org/licenses/by-nc-sa/4.0/ 2021

"""
###### History :  #######
#    - ??? : First release for Pecube V2 and V3
#    - 2021/07/16: Rewrite the script for Pecube V4 outputs
###### End History #######

###### To DO :  #######
#   - Test the MTL pdfs (I have not tested it, even if it should be OK)
#   - Add the plot of a projected profile between point A and point B ?
#	- Revise x & ylim
#	- 
#      
###### End To DO #######

from __future__ import  division
# This to be sure that the result of the division of integers is a real, not an integer

# Import modules
import os
import numpy as np
import matplotlib.pyplot as plt

#######################################
# Define data to plot
# dataplot: List of data to plot ; By default, the altitude will be plotted
#           Do not forget the simple quotes !!!
#           ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL']
dataplot = ['AHe','AFT']

# graphpath: name of the folder where the plot will be written
#            Usually you do not have to change it
graphpath = 'Graphs'

# Data to plot
#	datafnme:  results of Pecube forward modeling
#	inputdata: input data declared in Pecube.in
datafnme = 'Data/CompareAGE.csv'
inputdata = 'Data/Trujillo.csv'

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

# end define the data and parameters
#######################################

# Dictionaries definitions:
# You may have to change it depending on your Pecube version, and on your settings...
# agecol: respective column number of the data in the file comparison.txt
#         needs to be changed depending on your settings in the last paragraph 
#         of the file topo_parameters.txt
agecol = {'alt' : 'HEIGHT',
          'AHe' : 'AHE',
          'AFT' : 'AFT',
          'ZHe' : 'ZHE',
          'ZFT' : 'ZFT',
          'KAr' : 'KAR',
          'BAr' : 'BAR',
          'MAr' : 'MAR',
          'HbAr' : 'HAR',
          'FTL' : 'FT'
          }
# errname: respective column number of the error on data in the data input file
#		   Need to be updated for other 
errname = {'AHe' : 'DAHE',
		   'AFT' : 'DAFT',
		   'ZHe' : 'DZHE',
		   'ZFT' : 'DZFT',
		   'KAr' : 'DKAR',
           'BAr' : 'DBAR',
           'MAr' : 'DMAR',
           'HbAr' : 'DHAR'
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

# Earth radius in km
R = 6371

#################################################
#     Functions     #
def calc_length(XA,YA,XB,YB):       
	"""
	function to calcule the distance on a sphere
	between two points A and B given by their long/lat coordinates,   

	Read the page http://gis.stackexchange.com/questions/44064/how-to-calculate-distances-in-a-point-sequence

	"""

	R = 6371  # Earth diameter/2 (km)
	
	calc_lengthd = R * np.arccos(np.sin(YA * np.pi / 180) * np.sin(YB * np.pi / 180) + \
	                          np.cos(YA * np.pi / 180) * np.cos(YB * np.pi / 180) * \
	                          np.cos(-XA * np.pi / 180 + XB * np.pi / 180))
	return calc_lengthd
# -----------
       
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

#####################################

def project(A, B, datac):
	"""[summary]

	Args:
		A ([type]): [description]
		B ([type]): [description]
		datac ([type]): [description]
	"""

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

######### Main code ###############
if __name__ == "__main__":

	print('Plotting Pecube forward modelling results...')
	print('\t\tWritten by Xavier Robert, 07/2021\n')
	# Check if the Graphs/ folder exists, if not create it
	if os.path.exists(graphpath) == False:
		print('Output folder %s/ does not exist...' %(str(graphpath)))
		print('I am creating it...')
		os.mkdir(graphpath)
	else:
		print('Folder %s/ already exists, I will write in it and erase previous pdf files' %(str(graphpath)))
	
	########
	# Read data files
	datac = np.genfromtxt(fname = datafnme, delimiter = ',', names = True)
	if inputdata: inputc = np.genfromtxt(fname = inputdata, delimiter = ',', names = True)

	########
	if 'Longitude' in profiletype:
		print('\tPlotting longitudinal transect')
		# plot longitude transect
		# Loop on the data to plot
		fig1 = plt.figure()
		for item in dataplot:
			if inputdata:
				plt.errorbar(datac['LON'], datac[agecol[item]+'OBS'], yerr = inputc[errname[item]], 
			    	         fmt = 'o', label = agename[item], color = colores[item])
			else:
				plt.plot(datac['LON'], datac[agecol[item]+'OBS'], 
			    	     marker = 'o', linestyle = 'None', 
						 label = agename[item], color = colores[item])
		
			plt.plot(datac['LON'], datac[agecol[item]+'PRED'], 
		    	     marker = 's', linestyle = 'None', 
					 label = predname[item], color = colores[item], alpha = 0.3)

		# Remove -9999 values from the graph (= no data values)	
		plt.ylim(bottom = agerange[0], top = agerange[1])
	
		# set legend
		plt.legend()
		plt.xlabel('Longitude (°)')
		plt.ylabel('Age (Ma)')
		plt.title(graphtitle)

		plt.savefig(graphpath + '/' + graphtitle +'_long.pdf')
		fig1.clear()

	if 'Latitude' in profiletype:
		print('\tPlotting latitudinal transect')
		# plot latitude transect
		# Loop on the data to plot
		fig1 = plt.figure()
		for item in dataplot:
			if inputdata:
				plt.errorbar(datac['LAT'], datac[agecol[item]+'OBS'], yerr = inputc[errname[item]], 
			    	         fmt = 'o', label = agename[item], color = colores[item])
			else:
				plt.plot(datac['LAT'], datac[agecol[item]+'OBS'], 
			    	     marker = 'o', linestyle = 'None', 
						 label = agename[item], color = colores[item])
		
			plt.plot(datac['LAT'], datac[agecol[item]+'PRED'], 
		    	     marker = 's', linestyle = 'None', 
					 label = predname[item], color = colores[item], alpha = 0.3)

		# Remove -9999 values from the graph (= no data values)	
		plt.ylim(bottom = agerange[0], top = agerange[1])
	
		# set legend
		plt.legend()
		plt.xlabel('Latitude (°)')
		plt.ylabel('Age (Ma)')
		plt.title(graphtitle)

		plt.savefig(graphpath + '/' + graphtitle +'_lat.pdf')
		fig1.clear()

	if 'Altitude' in profiletype:
		print('\tPlotting Age-Elevation transect')
		# plot Altitude transect
		# Loop on the data to plot
		fig1 = plt.figure()
		for item in dataplot:
			if inputdata:
				plt.errorbar(datac[agecol[item]+'OBS'], datac['HEIGHTOBS'], xerr = inputc[errname[item]], 
			    	         fmt = 'o', label = agename[item], color = colores[item])
			else:
				plt.plot(datac[agecol[item]+'OBS'], datac['HEIGHTOBS'], 
			    	     marker = 'o', linestyle = 'None', 
						 label = agename[item], color = colores[item])
		
			#plt.plot(datac[agecol[item]+'PRED'], datac['HEIGHTOBS'], 
			plt.plot(datac[agecol[item]+'PRED'], datac['HEIGHTPRED'], 
		    	     marker = 's', linestyle = 'None', 
					 label = predname[item], color = colores[item], alpha = 0.3)

		# Remove -9999 values from the graph (= no data values)	
		plt.ylim(bottom = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']), 0),
		         top = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED'])))
		plt.xlim(left = agerange[0], right = agerange[1])
	
		# set legend
		plt.legend()
		plt.ylabel('Elevation (m)')
		plt.xlabel('Age (Ma)')
		plt.title(graphtitle)

		plt.savefig(graphpath + '/' + graphtitle +'_Alt.pdf')
		fig1.clear()

	if 'Projected' in profiletype and A and B:
		print('\tPlotting projected transect')
		# plot Projected transect
		# Compute projected coordinates
		coordproj = project(A, B, datac)
		if inputdata: coordprojinputc = project(A, B, inputc)

		# Loop on the data to plot
		fig1 = plt.figure()
		for item in dataplot:
			if inputdata:
				plt.errorbar(coordprojinputc, datac[agecol[item]+'OBS'], yerr = inputc[errname[item]], 
			    	         fmt = 'o', label = agename[item], color = colores[item])
			else:
				plt.plot(coordproj, datac[agecol[item]+'OBS'], 
			    	     marker = 'o', linestyle = 'None', 
						 label = agename[item], color = colores[item])
		
			plt.plot(coordproj, datac[agecol[item]+'PRED'], 
		    	     marker = 's', linestyle = 'None', 
					 label = predname[item], color = colores[item], alpha = 0.3)

		# Remove -9999 values from the graph (= no data values)	
		plt.ylim(bottom = agerange[0], top = agerange[1])
	
		# set legend
		plt.legend()
		plt.xlabel('Distance along transect (km)')
		plt.ylabel('Age (Ma)')
		plt.title(graphtitle)

		plt.savefig(graphpath + '/' + graphtitle +'_Proj.pdf')
		fig1.clear()

	########
	# Plot Altitude comparison
	print('\tPlotting altitude comparison')
	fig2 = plt.figure()
	plt.gca().set_aspect('equal')
	# Plot 1:1 line
	plt.plot([0,8000], [0,8000], marker = 'None', linestyle = '-', 
			 label = '1:1 line', color = 'lightgrey', alpha = 1)
	# Plot the altitude comparison
	plt.plot(datac[agecol['alt']+'OBS'], datac[agecol['alt']+'PRED'], 
			 marker = 'o', linestyle = 'None', 
			 label = agename[dataplot[0]], color = colores[dataplot[0]])
	# Write the 1:1 text over the line
	plt.text((max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED']))-0)/2,
			 (max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED']))-0)/2 - 
			      (max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED']))-0)/20,
			 '1:1', rotation = 45, color = 'lightgrey', alpha = 1, ha = 'center', va = 'center')

	plt.ylim(bottom = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']), 0), 
	         top = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED'])))
	plt.xlim(left = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']), 0), 
			 right = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED'])))

	plt.xlabel('Observed elevation (m)')
	plt.ylabel('Predicted elevation (m)')
	plt.title(graphtitle)

	plt.savefig(graphpath + '/' + graphtitle +'_Compare_Alt.pdf')
	fig2.clear()

	### End
