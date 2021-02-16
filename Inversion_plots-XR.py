######!/usr/bin/env python
# -*- coding: utf8 -*-
# coding: utf8

# Copyright (c) 2021 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later


"""
	!---------------------------------------------------------!
	!                                                         !
	!                   Nab 2 pdfs plots                      !
	!                                                         !
	!        Code to plot the NAB output into pdf plots       !
	!                                                         !
	!              Written by Xavier Robert                   !
	!                                                         !
	!---------------------------------------------------------!

	 This code is to plot scatter and pdf plots from the nab.out file
	 
     History:
        - version 0.0.1 (2021/02/11): initial preliminary release
        - 

	 TODOS:
        - 

"""
# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division
#from __future__ import unicode_literals
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
# To catch warnings
import warnings
warnings.simplefilter("error", UserWarning)


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
         'Last fault-point distance (km)',
         'First fault segment depth (km)',
         'Global block e-rate (km/Ma)',
         'Global block e-rate (km/Ma)',
         'Initiation for faulting (Ma)',
         'Fault velocity (km/Ma)']

# Set the couple of variables to plot against each other
# Exemple : Offset vs Slip rate 
# dataplot = [(1,3)] -->  plot=(1,3)
# If you plot 2D pdfs (contours), please, CHECK that the couple of parameter to plot
#     are the same and in the same order than in the nab.in file.
#     This python script DOES NOT check that!
dataplot = [(1,2), (3,4), (5,6), (7,8)]

# Give the name of the file with the inversion results
# Usually is NA_Results and nab.out
# /!\ THE FILE HAS TO BE LOCATED INSIDE THE "NA" FOLDER /!\
#inv_results = 'NA/NA_results.txt'
inv_results = 'NA/NA_results.csv'
data_nab = 'NA/NAB/nab.out'

# Choose if you want the PDFs (Probability Density Function) 
# plots 1D or 2D contour with the misfit plots
PDF_1D = 'yes'
pdf1d_results = 'NA/NAB/PDF_DATA.txt'

# Set the space between ticks for x and y axes for each parameters
# (same order than the list param)
#   If the tick format does not fit your variables, 
#   you may need to modify the dictionnary tick_order
#   in the function multiplot
tick_space =[2, 2, 25, 5, 0.02, 0.1, 5, 0.5]

# Set the size of the font for the x and y axes label
size_x = 15
size_y = 15
size_m = 15

#Set the size of the markers of the scatter plot and of the misfit
size_plo = 50
size_mis = 50

#Give the size of the number of bins for the 2d pdfs
size_2d = 100
PDF_2D = 'yes'
pdf2d_results = 'NA/NAB/PDF_2D_DATA.txt'

##########################################################################################
def read_NABout(data_nab, pdf1d_results, pdf2d_results, dataplot):
    """
    Read the nad.out file and write the PDF_DATA.txt and PDF_2D_DATA.txt
    
    INPUTS:
    	data_nab       : path and name of the nab.out file
        pdf1d_results  : path and name of the file with the PDF 1D data
        pdf2d_results  : path and name of the file with the PDF 2D data
        dataplot       : list of couple of variables we want to plot against each other
    
    OUTPUTS:
    	For the moment, no outputs
    		
    USAGE:
        read_NABout(data_nab, pdf1d_results, pdf2d_results, dataplot)
    
    Author: Xavier Robert, Grenoble 2021/02/10
    
    """
    # define test variable
    test = True

    # Open the nab.out file
    f1r = open(data_nab, 'r')

    # Read line to line
    lines = f1r.readlines()
    for line in lines:
        # Get number of parameters
        if u'Number of dimensions' in line: 
            Nparam = int(line[-13:len(line)])
            if test: print('Nparam: ', Nparam)
        if u'Number of bins per axis for 1D marginals' in line: 
            Nbins1D = int(line[-13:len(line)])
            if test: print('Nbins1D: ', Nbins1D)

        # Get number of parameters 2D
        if u'Number of 2D marginal pdfs to be calculated' in line: 
            Nparam2D = int(line[-13:len(line)])
            if test: print('Nparam2D: ', Nparam2D)
        if u'Number of bins per axis for 2D marginals' in line: 
            Nbins2D = int(line[-13:len(line)])
            if test: print('Nbins2D: ', Nbins2D)
    
    if os.path.isfile(pdf1d_results):
        print('\t\tWARNING: File %s already exists \n \t\tI do not erase it\n' %(pdf1d_results))
    else:
        # intitiate variable
        data1D = np.zeros((Nbins1D, 2*Nparam))
        for i in range(0, Nparam, 1):
            # Find indexes in lines of the
            #indexparam = lines.index('  Marginal for parameter :           ' + str(i+1) + '\n')
            indexparam = [item.replace(u' ', u'') for item in lines].index('Marginalforparameter:' + str(i+1) +'\n')
            #print(indexparam, indexparam2)
            # Copy the lines until the marginal i+1 in the variable ?
            for k in range (indexparam + 2, indexparam + 2 + Nbins1D):
                data1D[k - indexparam -2, 2*i:2*i+2] = lines[k].split()[0:2]
    
        np.savetxt(pdf1d_results, data1D)
        print('\tFile %s written\n' %(pdf1d_results))

    if os.path.isfile(pdf2d_results):
        print('\t\tWARNING: File %s already exists \n \t\tI do not erase it\n' %(pdf2d_results))
    else:
        # intitiate variable
        #data2D = np.zeros((Nbins2D, 4*Nparam2D))
        data2D = np.zeros((Nbins2D, Nbins2D*len(dataplot)))
        for i in range (0, len(dataplot)):
            # Find indexes in lines of the
            try:
                #indexparam = lines.index('  2D Marginal for parameters :           ' + str(dataplot[i][0]) + '           ' + str(dataplot[i][1]) + '\n')
                indexparam = [item.replace(u' ', u'') for item in lines].index('2DMarginalforparameters:' + str(dataplot[i][0]) + str(dataplot[i][1]) + '\n')
            except ValueError:
                raise ValueError('2D pdf couple asked do not agree with 2D pdf computed (%s, %s) ' %(str(dataplot[i][0]), str(dataplot[i][1])))
            # Copy the lines until the marginal i+1 in the variable ?
            for k in range (indexparam + 4, indexparam + 4 + Nbins2D):
                data2D[k - indexparam - 4, (Nbins2D * i):(Nbins2D * i + Nbins2D)] = lines[k].split()

        np.savetxt(pdf2d_results, data2D)
        print('\tFile %s written \n' %(pdf2d_results))

    # Close the nab.out file
    f1r.close()

    # find parameter ranges (Uncomment the block if needed)
    #param_ranges = np.zeros((Nparam, 3))
    #indexparam = lines.index('  Parameter ranges\n')
    #for i in range (0, Nparam):
    #    param_ranges[i, 0:3] = lines[indexparam+2+i].split()[0:3]
    #
    #return param_ranges

    return

##########################################################################################
def multiplot(param, nb_var, plot,
              tick_space, 
              pdf1d_results, pdf2d_results,
              size_x = 15, size_y = 15, size_m = 15, size_plo = 50, size_mis = 50, 
              size_2d = 100, PDF_1D = 'yes', PDF_2D = 'yes',
              i_param = None):
    """[summary]

    Args:
        param ([type]): [description]
        tick_space ([type]): [description]
        inv_results ([type]): [description]
        plot ([type]): [description]
        i_param ([type], optional): [description]. Defaults to None.
    """

    # Define dictionnary to manage tick format
    tick_order = {1000  : u'%4.0f',
                  100   : u'%4.0f',
                  10    : u'%4.0f',
                  1     : u'%4.1f',
	              0.1   : u'%4.1f',
                  0.01  : u'%4.2f',
                  0.001 : u'%4.3f'}

    # New names for the variables to make the code easy to read
    misfit = nb_var[0]
    value_x = nb_var[plot[0]]
    value_y = nb_var[plot[1]]
    
    # Set up x_ and y_tick_space
    x_tick_space = tick_space[plot[0]-1]
    y_tick_space = tick_space[plot[1]-1]

    # Plot the results as a scatter plot
    fig, ax = plt.subplots(figsize = (7, 6))

    # set the misfit value as a variable for the color
    scatter_mis = ax.scatter(value_x, 
                             value_y, 
                             s = size_mis,            
                             linewidth = 0.1,   
                             edgecolor = 'k',   
                             c = np.log10(misfit),  # data use for the color bar, here as a log10
                             cmap = 'Spectral')
                           
    # Find the lowest misfit run and return the line/run number
    if np.isnan(min(misfit)):
        raise ValueError('ERROR: one of the misifit is NaN; \n\tthis is probably because you had an issue with the multithread inversion process...\n\tRe-run the inversion with 1 core only.')
        # This error is known with clusters using OAR protocols.
        # It seems that there are errors when nodes discuss with themselves, 
        # Some misfit became NaN, and the misfits printed in the results files 
        # does not trully correspond to the parameters values associated in those files.
        # BE CAREFULL !
    else:
        bfit = np.where(misfit == min(misfit))[0][0]
        print("Best model is number " + str(bfit + 1))
                
    # Plot the best (lowest misfit) run as a yellow star
    plt.plot(value_x[bfit],
             value_y[bfit],
             '*',
             markeredgewidth = 1,
             markeredgecolor = 'k',
             markersize = 25,
             c = '#ffff00')

    # Set the misfit legend as a color bar
    cbar = plt.colorbar(scatter_mis, pad = 0.025)
    cbar.ax.set_ylabel('Misfit',
                       rotation = 270,
                       size = size_m,
                       labelpad = 20)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
    
    ################
    # Axes settings
    ################

    # Set the starting value of the axesnumber of ticks on the axes 
    plt.xticks(np.arange(value_x.min(),
                         value_x.max()),
                         fontsize = 10)
        
    plt.yticks(np.arange(value_y.min(),
                         value_y.max()),
                         fontsize = 10)
    
    # find the number of digits to print for the ticks
    for x_order in (1000, 100, 10, 1, 0.1, 0.01, 0.001):
        if int(x_tick_space/x_order) in range (1,10):
            ax.xaxis.set_major_formatter(FormatStrFormatter(tick_order[x_order]))
        if int(y_tick_space/x_order) in range (1,10):
            ax.yaxis.set_major_formatter(FormatStrFormatter(tick_order[x_order]))
    ax.xaxis.set_major_locator(mpltick.MultipleLocator(x_tick_space))
    ax.yaxis.set_major_locator(mpltick.MultipleLocator(y_tick_space))

    # set the name of the axes
    plt.xlabel(param[plot[0]-1], fontsize=size_x)
    plt.ylabel(param[plot[1]-1], fontsize=size_y)
    
    #################
    # PDFs Histogram plot code
    #################
    
    # Load the na_bayes data
    if PDF_2D == 'yes' or PDF_1D == 'yes': 
        pdf_data = np.loadtxt(pdf1d_results, unpack = True)

    #Start 1D PDF or/and 2D PDF    
    if PDF_2D == 'yes':
        print('Build 2d graphs')  

        # Load the 2D pdf data and assing variables to it
        val_z2d = np.zeros((size_2d, size_2d)) # empty array to stock the pdf results
        data_2pdf = np.loadtxt(pdf2d_results, unpack = True)        
        val_x2d = pdf_data[(plot[0]-1)*2]
        val_y2d = pdf_data[(plot[1]-1)*2]
        if i_param or i_param == 0:
            val_z2d[0:100, 0:100] = data_2pdf[(i_param+1) * np.shape(data_2pdf)[1] - np.shape(data_2pdf)[1]:(i_param+2) * np.shape(data_2pdf)[1] - np.shape(data_2pdf)[1], :]
        else:
            val_z2d = data_2pdf

        val_xx2d, val_yy2d = np.meshgrid(val_x2d, val_y2d)

        try:
            CS = plt.contour(val_xx2d, val_yy2d, val_z2d.T, 
                             levels = [0.4], cmap = 'CMRmap')
            plt.clabel(CS, inline = 2, fontsize = 12)
        except UserWarning:
            print('WARNING: No contour level')
    
    # Start 1D PDF or/and 2D PDF    
    if PDF_1D == 'yes':
        print('build 1d graphs')
        
        # Position the 2 new figures     
        divider = make_axes_locatable(ax)
        axHistx = divider.append_axes("top", 0.7, pad = 0.3, sharex = ax)
        axHisty = divider.append_axes("right", 0.7, pad = 0.3)
    
        # Disable the labels for the axes they share with the scatter plot    
        plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
                 visible = False)
    
        # Loop to select the right columns for the right parameters
        if plot[0] == 1:
            axHistx.plot(pdf_data[0], pdf_data[1])
        else:
            axHistx.plot(pdf_data[plot[0] + (plot[0] - 2)], 
                         pdf_data[plot[0] + (plot[0] - 1)])
        if plot[1] == 1:
            axHisty.plot(pdf_data[1], pdf_data[0])
        else:
            
            axHisty.plot(pdf_data[plot[1] + (plot[1] - 1)], 
                         pdf_data[plot[1] + (plot[1] - 2)])
          
            axHistx.set_ylim(ymin = 0)
            axHisty.set_xlim(xmin = 0)
    #param_1 = param[plot[0]-1]
    #param_2 = param[plot[1]-1]
    
    #nab_best=np.loadtxt('NA/nab_best.txt', unpack = True)
    
    # Saving the plots as a pdf file
    if i_param or i_param == 0:
        plt.savefig('NA/Graphs/PDF_' + str(i_param+1) + '.pdf')
        print('Plotting : ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (NA/Graphs/PDF_' + str(i_param+1) + '.pdf)')
    else: 
        plt.savefig('NA/Graphs/PDF.pdf')
        print('Plotting : ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (NA/Graphs/PDF.pdf)')
    #plt.show()
    
    return




##########################################################################################
if __name__ == u'__main__':
    """
    main code
    """

    # Check that nab.out exists
    if not os.path.isfile(data_nab): raise NameError('ERROR : File {FileNa} does not exist'.format(FileNa=str(data_nab)))
    # Check that NA_results exists
    if not os.path.isfile(inv_results): raise NameError('ERROR : File {FileNa} does not exist'.format(FileNa=str(inv_results)))
    if not os.path.exists('NA/Graphs'): os.mkdir('NA/Graphs')

    # Build Nab files for plotting
    read_NABout(data_nab, pdf1d_results, pdf2d_results, dataplot)

    # Loading the data, calling as many variables as given in param + 1 for the misfit
    if inv_results[-3:] == 'txt':
        nb_var = np.zeros(len(param) + 1)
        nb_var = np.loadtxt(inv_results, unpack = True)
    elif inv_results[-3:] == 'csv':
        nb_var = np.genfromtxt(inv_results, delimiter = ',', skip_header = 1, unpack = True)
        np.savetxt('test.txt', nb_var)
    else:
        raise ImportError('NA_results format not supported...\n\n')

    # plot data
    if len(dataplot) > 1:
        for i in range(len(dataplot)):
            print()
            multiplot(param, nb_var, dataplot[i], 
                      tick_space, pdf1d_results, pdf2d_results, 
                      size_x, size_y, size_m, size_plo, size_mis, 
                      size_2d, PDF_1D, PDF_2D, i)
    else:
        print()
        multiplot(param, nb_var, dataplot[0],
                  tick_space, pdf1d_results, pdf2d_results, 
                   size_x, size_y, size_m, size_plo, size_mis, 
                   size_2d, PDF_1D, PDF_2D)
