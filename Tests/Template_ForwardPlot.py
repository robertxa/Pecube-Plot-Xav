######!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyPlotPecube import pyPlotPecubeForward as Fplot

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

# graphtitle: title to write on the graph
graphtitle = 'Trujillo transect'

#agerange = range of the ages to plot on the profiles
#			[min, max]
agerange = [0, 100]

# profiletype: type of profile = ['Latitude', 'Longitude', 'Altitude', 'Projected']
#              if [], no age profile is plotted
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

Fplot.PlotPecubeForward(dataplot = dataplot,
				  graphpath = graphpath,
				  datafnme = datafnme,
				  inputdata = inputdata,
				  graphtitle = graphtitle,
				  agerange = agerange,
				  profiletype = profiletype,
				  A = A, B = B)