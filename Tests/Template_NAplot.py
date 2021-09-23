######!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyPlotPecube import pyPlotPecubeNA as NAplot

    
############# User defined variables #################

#Define as many variable as you have, with their unit
# Check the order in na.sum (open with text editor) or in NA_Results
# If you use the later, first column always is the misfit
# The others paramaters are from topo_parameters.txt then from fault_parameters.txt
# Exemple : param=['Offset (km)','Basal Temperature (°C)','Slip rate (km∕Ma)']

# /!\ IF YOU WANT TO USE A SLASH USE THIS ONE --> '∕' <-- , IT'S A UNICODE DIVISION SYMBOL
# WINDOWS AND OSX DON'T ALLOW THE USE OF THE REGULAR SLASH 

param = ['Ignimbrite filling time (Ma)', 
         'Initiation for ignimbrite carving (Ma)',
         'x1 fault (20 km depth) (km)',
         'x2 fault (4 km high) (km)',
         'Block exhumation 39-14 (km/Ma)',
         'Block exhumation 14-0 (km/Ma)',
         'Initiation for faulting (Ma)',
         'Fault velocity (km/Ma)']

# Set the couple of variables to plot against each other
# Exemple : Offset vs Slip rate 
# dataplot = [(1,3)] -->  plot=(1,3)
# If you plot 2D pdfs (contours), please, CHECK that the couple of parameter to plot
#     are the same and in the same order than in the nab.in file.
#     This python script check it and will insult you if this is not compatible !!!
dataplot = [(1,2), (3,4), (5,6), (7,8)]

# Give the name of the file with the inversion results
# Usually is NA_Results and nab.out
#inv_results = 'NA/NA_results.txt'
inv_results = 'NA/NA_results.csv'
data_nab    = 'NA/NAB/nab.out'
graph_path  = 'NA/Graphs/'    # Do not forget the '/' at the end.

# Choose if you want the PDFs (Probability Density Function) 
# plots 1D or 2D contour with the misfit plots (True of False)
PDF_1D = True
PDF_2D = False

# Set the space between ticks for x and y axes for each parameters
# (same order than the list param)
#   If the tick format does not fit your variables, 
#   you may need to modify the dictionnary tick_order
#   in the function multiplot
tick_space =[2, 2, 25, 5, 0.05, 0.1, 5, 0.5]

# Set the size of the font for the x and y axes label
size_x = 15
size_y = 15
size_m = 15

# Set the size of the markers of the misfit
size_mis = 50

# Parameters to find peaks
peak_thres = 0.05
peak_min_dist = 30

# If you want to print the pdfs in a text file, just modify the 2 next lines
pdf1d_results = None #'NA/NAB/PDF_DATA.txt'
pdf2d_results = None #'NA/NAB/PDF_2D_DATA.txt'
    
    
NAplot.NAplot_Pecube(param = param,
              dataplot = dataplot,
              inv_results = inv_results,
              data_nab = data_nab,
              graph_path  = graph_path,
              PDF_1D = PDF_1D, PDF_2D = PDF_2D,
              pdf1d_results = pdf1d_results, pdf2d_results = pdf2d_results,
              tick_space = tick_space,
              size_x = size_x, size_y = size_y, size_m = size_m,
              size_mis = size_mis,
              peak_thres = peak_thres, peak_min_dist = peak_min_dist
              )