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
        - version 1.0.1 (2021/02/18): initial preliminary release
        - 

	 TODOS:
        - 
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

#Set the size of the markers of the scatter plot and of the misfit
size_plo = 50
size_mis = 50

# If you want to print the pdfs in a text file, just modify the 2 next lines
pdf1d_results = None #'NA/NAB/PDF_DATA.txt'
pdf2d_results = None #'NA/NAB/PDF_2D_DATA.txt'

##########################################################################################
def read_NABout(data_nab, dataplot, pdf1d_results = None, pdf2d_results = None):
    """
    Read the nad.out file and write the PDF_DATA.txt and PDF_2D_DATA.txt
    
    INPUTS:
    	data_nab       : path and name of the nab.out file
        dataplot       : list of couple of variables we want to plot against each other
        pdf1d_results  : path and name of the file with the PDF 1D data if file saved ;
                         None by default
        pdf2d_results  : path and name of the file with the PDF 2D data if file saved ;
                         None by default
    
    OUTPUTS:
    	data1D : 1D marginals
        data2D : 2D marginals
    		
    USAGE:
        data1D, data2D = read_NABout(data_nab, dataplot, pdf1d_results, pdf2d_results)
    
    Author: Xavier Robert, Grenoble 2021/02/10
    
    """

    # Open the nab.out file and store the content in lines
    with open(data_nab, 'r') as f1r:
        lines = f1r.readlines()

    test_couple = []
    for line in lines:
        # Get number of parameters
        if u'Number of dimensions' in line: 
            Nparam = int(line[-13:len(line)])
        if u'Number of bins per axis for 1D marginals' in line: 
            Nbins1D = int(line[-13:len(line)])

        # Get number of parameters 2D
        if u'Number of 2D marginal pdfs to be calculated' in line: 
            Nparam2D = int(line[-13:len(line)])
            # Check if the number of plots asked is compatible with the number of 2D marginals computed
            if len(dataplot) > Nparam2D:
                raise ValueError('ERROR: Number of plots (%s) > Number of 2D marginals computed (%s)' %(len(dataplot), Nparam2D))
        # Record the 2D marginals couple present in the nab.out file
        if u'2D Marginal for parameters :' in line:
            test_couple.append((int(line.split()[-2]), int(line.split()[-1])))
        # Find the number of bins for 2D marginals
        if u'Number of bins per axis for 2D marginals' in line: 
            Nbins2D = int(line[-13:len(line)])

    # check if the couples asked to plot are OK with the marginals computed
    for elem in dataplot:
        if elem not in test_couple:
            raise ValueError('ERROR: F*** You did nto read all the comments in the beginning !!\n2D plot %s not in 2D marginals computed...\n\t\tPlease, rerun NAD with correct couples\n' %(str(elem)))
    
    # intitiate variable
    data1D = np.zeros((Nbins1D, 2*Nparam))
    for i in range(0, Nparam, 1):
        # Find index of the beginning of the last iteration for the parameter i+1
        zzz = 0
        indexparam = 0
        for zzz, z_str in enumerate([item.replace(u' ', u'') for item in lines]):
            if z_str == 'Marginalforparameter:' + str(i+1) +'\n':
                indexparam = zzz
        # Copy the lines until the marginal i+1 in the variable data1D
        for k in range (indexparam + 2, indexparam + 2 + Nbins1D):
            data1D[k - indexparam -2, 2*i:2*i+2] = lines[k].split()[0:2]
    
    # If asked, save the data1D variable in a text file
    if pdf1d_results:
        if os.path.isfile(pdf1d_results):
            print('\t\tWARNING: File %s already exists \n \t\tI do not erase it\n' %(pdf1d_results))
        else:
            np.savetxt(pdf1d_results, data1D)
            print('\tFile %s written\n' %(pdf1d_results))

    # intitiate variable
    data2D = np.zeros((Nbins2D, Nbins2D*len(dataplot)))
    for i in range (0, len(dataplot)):
        # Find index of the beginning of the last iteration for the couple of parameters in dataplot
        try:
            # We kept this line because it permits to raise an error we do not find the couple do plot.
            # Normaly, this is already done before the call to thsi function...
            #indexparam = [item.replace(u' ', u'') for item in lines].index('2DMarginalforparameters:' + str(dataplot[i][0]) + str(dataplot[i][1]) + '\n')
            zzz = 0
            indexparam = 0
            # Find index of the beginning of the last iteration for the couple of parameters
            for zzz, z_str in enumerate([item.replace(u' ', u'') for item in lines]):
                if z_str == '2DMarginalforparameters:' + str(dataplot[i][0]) + str(dataplot[i][1]) + '\n':
                    indexparam = zzz
        except ValueError:
            raise ValueError('2D pdf couple asked do not agree with 2D pdf computed (%s, %s) ' %(str(dataplot[i][0]), str(dataplot[i][1])))
        # Copy the lines until the marginal i+1 in the variable ?
        for k in range (indexparam + 4, indexparam + 4 + Nbins2D):
            data2D[k - indexparam - 4, (Nbins2D * i):(Nbins2D * i + Nbins2D)] = lines[k].split()
    
    # If asked, save the data2D variable in a text file
    if pdf2d_results:
        if os.path.isfile(pdf2d_results):
            print('\t\tWARNING: File %s already exists \n \t\tI do not erase it\n' %(pdf2d_results))
        else:
            np.savetxt(pdf2d_results, data2D)
            print('\tFile %s written \n' %(pdf2d_results))

    # Find parameter ranges (Uncomment the block if needed)
    #param_ranges = np.zeros((Nparam, 3))
    #indexparam = lines.index('  Parameter ranges\n')
    #for i in range (0, Nparam):
    #    param_ranges[i, 0:3] = lines[indexparam+2+i].split()[0:3]
    #
    #return data1D, data2D, param_ranges

    return data1D, data2D



##########################################################################################
def multiplot(param, nb_var, plot,
              tick_space, 
              data1D = None, data2D = None,
              size_x = 15, size_y = 15, size_m = 15, size_plo = 50, size_mis = 50, 
              PDF_1D = True, PDF_2D = True,
              i_param = None):
    """
    Function to plot the scatter plots with the 1D and 2D pdfs if asked.

    INPUTS:
        param       : name of the parameters
        nb_var      : NA_Results data
        plot        : parameter to plot
        tick_space  : Space between the ticks
        data1D      : 1D marginal. Defaults to None.
        data2D      : 2D marginal. Defaults to None.
        size_x      : Font size for x-axis. Defaults to 15.
        size_y      : Font size for y-axis. Defaults to 15.
        size_m      : Font size for ???. Defaults to 15.
        size_plo    : Size of the plot. Defaults to 50.
        size_mis    : Size of the plot. Defaults to 50.
        PDF_1D      : Boolean to tell if we plot 1D-pdf (True) or not (False).
                      Defaults to True. Optional.
        PDF_2D      : Boolean to tell if we plot 2D-pdf (True) or not (False). 
                      Defaults to True. Optional.
        i_param     : Index of couple of parameters to plot. 
                      Defaults to None. Optional.

    OUTPUTS:
        Just pdf graphs !
    """

    # Define dictionnary to manage tick format
    # You may add you own settings if needed
    tick_order = {1000  : u'%4.0f',
                  100   : u'%4.0f',
                  10    : u'%4.0f',
                  1     : u'%4.1f',
	              0.1   : u'%4.1f',
                  0.01  : u'%4.2f',
                  0.001 : u'%4.3f'}

    # Plot the results as a scatter plot
    fig, ax = plt.subplots(figsize = (7, 6))

    # set the misfit value as a variable for the color
    scatter_mis = ax.scatter(nb_var[plot[0]], 
                             nb_var[plot[1]], 
                             s = size_mis,            
                             linewidth = 0.1,   
                             edgecolor = 'k',   
                             c = np.log10(nb_var[0]),  # data use for the color bar, here as a log10
                             cmap = 'Spectral')
                           
    # Find the lowest misfit run and return the line/run number
    if np.isnan(min(nb_var[0])):
        raise ValueError('ERROR: one of the misifit is NaN; \n\tthis is probably because you had an issue with the multithread inversion process...\n\tRe-run the inversion with 1 core only.')
        # This error is known with clusters using OAR protocols.
        # It seems that there are errors when nodes discuss with themselves, 
        # Some misfit became NaN, and the misfits printed in the results files 
        # does not trully correspond to the parameters values associated in those files.
        # BE CAREFULL !
    else:
        bfit = np.where(nb_var[0] == min(nb_var[0]))[0][0]
        print("Best model is number %s with misfit of %s" %(str(bfit + 1), str(min(nb_var[0]))))
                
    # Plot the best (lowest misfit) run as a yellow star
    plt.plot(nb_var[plot[0]][bfit],
             nb_var[plot[1]][bfit],
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
    # Find the right misfit digits to print
    for x_order in (1000, 100, 10, 1, 0.1, 0.01, 0.001):
        if int((np.log(max(nb_var[0]))-np.log(min(nb_var[0])))/x_order) in range (1,10):
            cbar.ax.yaxis.set_major_formatter(FormatStrFormatter(tick_order[x_order/10]))
    #cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
    
    ################
    # Axes settings
    ################

    # Set the starting value of the axesnumber of ticks on the axes 
    plt.xticks(np.arange(nb_var[plot[0]].min(),
                         nb_var[plot[0]].max()),
                         fontsize = 10)
        
    plt.yticks(np.arange(nb_var[plot[1]].min(),
                         nb_var[plot[1]].max()),
                         fontsize = 10)
    
    # find the number of digits to print for the ticks
    for x_order in (1000, 100, 10, 1, 0.1, 0.01, 0.001):
        if int(tick_space[plot[0]-1]/x_order) in range (1,10):
            ax.xaxis.set_major_formatter(FormatStrFormatter(tick_order[x_order]))
        if int(tick_space[plot[1]-1]/x_order) in range (1,10):
            ax.yaxis.set_major_formatter(FormatStrFormatter(tick_order[x_order]))
    ax.xaxis.set_major_locator(mpltick.MultipleLocator(tick_space[plot[0]-1]))
    ax.yaxis.set_major_locator(mpltick.MultipleLocator(tick_space[plot[1]-1]))

    # set the name of the axes
    plt.xlabel(param[plot[0]-1], fontsize=size_x)
    plt.ylabel(param[plot[1]-1], fontsize=size_y)
    
    ############################
    # PDFs Histogram plot code #
    ############################

    # Start 1D PDF or/and 2D PDF    
    if PDF_2D:
        print('Build 2D pdfs countours')  

        # Load the 2D pdf data and assing variables to it
        val_z2d = np.zeros((np.shape(data2D.T)[1], np.shape(data2D.T)[1]))  # empty array to stock the pdf results
        val_x2d = data1D.T[(plot[0]-1)*2]
        val_y2d = data1D.T[(plot[1]-1)*2]
        if i_param or i_param == 0:
            val_z2d[0:100, 0:100] = data2D.T[(i_param+1) * np.shape(data2D.T)[1] - np.shape(data2D.T)[1]:(i_param+2) * np.shape(data2D.T)[1] - np.shape(data2D.T)[1], :]
        else:
            val_z2d = data2D.T

        val_xx2d, val_yy2d = np.meshgrid(val_x2d, val_y2d)

        try:
            CS = plt.contour(val_xx2d, val_yy2d, val_z2d.T, 
                             #levels = [0.4], cmap = 'CMRmap')
                             #levels = 3, cmap = 'CMRmap')
                             #levels = 3, cmap = 'gray_r')
                             levels = 3, cmap = 'hot')
            plt.clabel(CS, inline = 2, fontsize = 12)
        except UserWarning:
            print('WARNING: No contour levels plotted')
    
    # Start 1D PDF or/and 2D PDF    
    if PDF_1D:
        print('Build 1D pdfs')
        
        # Position the 2 new figures     
        divider = make_axes_locatable(ax)
        axHistx = divider.append_axes("top", 0.7, pad = 0.3, sharex = ax)
        axHisty = divider.append_axes("right", 0.7, pad = 0.3)
    
        # Disable the labels for the axes they share with the scatter plot    
        plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
                 visible = False)
    
        # Loop to select the right columns for the right parameters
        if plot[0] == 1:
            axHistx.plot(data1D.T[0], data1D.T[1])
        else:
            axHistx.plot(data1D.T[plot[0] + (plot[0] - 2)], 
                         data1D.T[plot[0] + (plot[0] - 1)])
        if plot[1] == 1:
            axHisty.plot(data1D.T[1], data1D.T[0])
        else:
            axHisty.plot(data1D.T[plot[1] + (plot[1] - 1)], 
                         data1D.T[plot[1] + (plot[1] - 2)])
          
            axHistx.set_ylim(ymin = 0)
            axHisty.set_xlim(xmin = 0)
    
    # Saving the plots as a pdf file
    if i_param or i_param == 0:
        plt.savefig('NA/Graphs/PDF_' + str(i_param+1) + '.pdf')
        print('Plotting : ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (NA/Graphs/PDF_' + str(i_param+1) + '.pdf)\n')
    else: 
        plt.savefig('NA/Graphs/PDF.pdf')
        print('Plotting : ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (NA/Graphs/PDF.pdf)\n')
    
    # If you want to plot on screen each plot
    #plt.show()
    
    return




##########################################################################################
if __name__ == u'__main__':
    """
    main code
    """

    print('###########################################################################\n')
    print('\t\tPlot results from NA inversions\n')
    print('###########################################################################\n')
    # Check that NA_results exists
    if not os.path.isfile(inv_results):
        raise NameError('ERROR: File %s does not exist' %(str(inv_results)))
    # Check if the number of parameters is compatible with the number ticks_space
    if len(param) != len(tick_space):
        raise ValueError('ERROR: the number of parameters is different than the number of ticks_space asked !\n Check it !')
    
    # Loading the data, calling as many variables as given in param + 1 for the misfit
    if inv_results[-3:] == 'txt':
        nb_var = np.loadtxt(inv_results, unpack = True)
    elif inv_results[-3:] == 'csv':
        nb_var = np.genfromtxt(inv_results, delimiter = ',', skip_header = 1, unpack = True)
    else:
        raise ImportError('NA_results format not supported...\n\n')
    
    # Check if misfits are not with NaN or infinite values...
    # Raise error if this is the case, and recheck your inversion parameters
    if np.any(np.isnan(nb_var)) or not np.all(np.isfinite(nb_var)):
        print("\nInfinity: ", np.inf in nb_var)
        raise ValueError('NA_results.csv contains NaN or infinite numbers...\n\n')

    # Sort nb_var in function of misfit to get the best misfits over the lower ones 
    #          --> Better scatter plot visualisation
    nb_var = nb_var.T[nb_var.T[:,0].argsort()[::-1]].T
    #nb_var.T[:] = nb_var.T[::-1]

    # Check that nab.out exists
    if not os.path.isfile(data_nab): 
        print ('WARNING: File %s does not exist\n\n' %(str(data_nab)))
        PDF_1D = False
        PDF_2D = False
        
    if not PDF_1D: pdf1d_results = None
    if not PDF_2D: pdf2d_results = None

    # Build Nab files for plotting
    if PDF_1D or PDF_2D:
        data1D, data2D = read_NABout(data_nab, dataplot, pdf1d_results, pdf2d_results)

    # Plot data
    # Check if there is a folder to store the outputs
    if not os.path.exists('NA/Graphs'): os.mkdir('NA/Graphs')
    if len(dataplot) > 1:
        for i in range(len(dataplot)):
            # Call the plot function for each parameter value to plot
            multiplot(param, nb_var, dataplot[i], 
                      tick_space, 
                      data1D, data2D,
                      size_x, size_y, size_m, size_plo, size_mis, 
                      PDF_1D, PDF_2D, i)
    else:
        multiplot(param, nb_var, dataplot[0],
                  tick_space, 
                  data1D, data2D,
                  size_x, size_y, size_m, size_plo, size_mis, 
                  PDF_1D, PDF_2D)

    print('###########################################################################\n')
    print('\t\tEnd...\n')
    print('###########################################################################\n')
