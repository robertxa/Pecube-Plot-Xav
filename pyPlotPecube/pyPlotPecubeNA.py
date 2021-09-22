######!/usr/bin/env python
# -*- coding: utf8 -*-
# coding: utf8

# Copyright (c) 2021 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later

################################################################
# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division
#from __future__ import unicode_literals
import sys, os
import numpy as np
from scipy.optimize import curve_fit    # To find gaussian parameters
import peakutils    # To find number of peaks
from peakutils.plot import plot as pplot
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
# To catch warnings
import warnings
warnings.simplefilter("error", UserWarning)

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
              graph_path = 'NA/Graphs/',
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
        graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                     Do not forget the '/' at the end.
                                     Defaults to 'NA/Graphs/'.
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

    if i_param or i_param == 0:
        print('\tPlotting : ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (' + graph_path + 'PDF_' + str(i_param+1) + '.pdf)')
    else: 
        print('\tPlotting : ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (' + graph_path + 'PDF.pdf)')

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
        raise ValueError('\t\tERROR: one of the misifit is NaN; \n\tthis is probably because you had an issue with the multithread inversion process...\n\tRe-run the inversion with 1 core only.')
        # This error is known with clusters using OAR protocols.
        # It seems that there are errors when nodes discuss with themselves, 
        # Some misfit became NaN, and the misfits printed in the results files 
        # does not trully correspond to the parameters values associated in those files.
        # BE CAREFULL !
    else:
        bfit = np.where(nb_var[0] == min(nb_var[0]))[0][0]
        print("\t\tBest model is number %s with misfit of %s" %(str(bfit + 1), str(min(nb_var[0]))))
                
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
        print('\t\tBuild 2D pdfs countours')  

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
            print('\t\tWARNING: No contour levels plotted')
    
    # Start 1D PDF or/and 2D PDF    
    if PDF_1D:
        print('\t\tBuild 1D pdfs')
        
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
        plt.savefig(graph_path + 'PDF_' + str(i_param+1) + '.pdf')
    else: 
        plt.savefig(graph_path + 'Graphs/PDF.pdf')
    print(' ')
    
    return


##########################################################################################
def gaussian(x, ampg, meang, sigma1g):
    """
    Function to calculate the Gaussian with constants ampg, meang, and sigmag

    Args:
        x ([float])      : data to describe
        ampg ([float])   : amplitude of the Gaussian
        meang ([float])  : mean of the Gaussian
        sigma1g ([float]): 1-sigma1 of the Gaussian

    Returns:
        1 Gaussian pdf
    """
    
    return ampg * np.exp(-np.power(x - meang, 2)/(2 * np.power(sigma1g, 2)))
    

##########################################################################################
def gaussians(x, *gaussians_param):#amp1, cen1, sigma1, amp2, cen2, sigma2):
    """
    Compute the sum of several Gaussian functions

    INPUTS:
        x ([float])                               : data to desccribe
        gaussians_param (list of floats triplets) : list of tripplets that describes each gaussian pdf
                                                    with the amplitude, the mean and the 1 sigma for each pdf

    RETURNS:
        gaussians_results (np.array of floats): results of the sum of several Gaussian functions.
                                                Th size is the same than the input vector x.
    """
    
    # Clear the variable
    gaussians_results = np.zeros(x.shape)
    # Do a loop on the number of peaks --> Determines the number of gaussians to stack
    for h in range(0, int(len(gaussians_param)/3)):
        gaussians_results = gaussians_results + gaussian(x, gaussians_param[3*h], 
                                                            gaussians_param[3*h+1], 
                                                            gaussians_param[3*h+2])

    return gaussians_results


##########################################################################################
def plotgfit(data1D, i, pars, ipeak = None):
    """
    Plot the gaussian fitted with the peak and the 1 sigma error

    Args:
        data1D ([type]): [description]
        i ([type]): [description]
        pars ([type]): [description]
        ipeak ([type], optional): [description]. Defaults to None.
    """
    
    if ipeak == 1 : labelname = '1-peak Gaussian fit'
    else: labelname = None
    plt.plot(data1D.T[2*i], 
            gaussian(data1D.T[2*i], pars[0], pars[1], pars[2]),
            "-b",
            label = labelname)

    # Colorize the gaussian between the 1 sigmas
    if ipeak == 1 : labelname = 'Acceptable values\n' + '(1\u03C3)'
    else: labelname = None
    plt.fill_between(data1D.T[2*i], 
                     gaussian(data1D.T[2*i], *pars),
                     0,
                     where = ((data1D.T[2*i] >= pars[1] - pars[2]) & (data1D.T[2*i] <= pars[1] + pars[2])),
                     alpha=0.30, 
                     color='green', 
                     interpolate=True,
                     label = labelname)
            
    # Plot the mean
    plt.vlines(x = pars[1], 
               ymin = 0, 
               ymax = max(gaussian(data1D.T[2*i], *pars)),
               colors = 'red',
               label = 'Mean value %s (%0.2f %s)\n%0.2f +/- %0.2f' %(ipeak, pars[0]*100, chr(37), pars[1], pars[2]))

    # End of plot
    return


##########################################################################################
def statsparam(data1D, param, graph_path  = 'NA/Graphs/', size_x = 15, size_y = 15):
    """
    Function to find the number of gaussian functions, to compute their mean and sigma,
    and to plot the gaussians fit

    INPUTS:
        data1D ([np.array])   : Array of 1D-pdfs for each parameter
        param (list of string): List the name of the parameters 
                                in the same order than define in Pecube inversion
        graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                     Do not forget the '/' at the end.
                                     Defaults to 'NA/Graphs/'.
        size_x (int, optional): Size of the font for the x axes label. 
                                Defaults to 15.
        size_y (int, optional): Size of the font for the y axes label. 
                                Defaults to 15.
    """
    
    # Open the text file where to put the results
    fstats_w = open(graph_path + 'Nab-stats.txt', 'w')
    # Write the Header
    fstats_w.write('Param \t Mean \t Mean_stdev \t Std_Err \t Std_Err_Stdev \n')

    # Faire la boucle sur i in range (0, len(param)):
    for i in range (0, len(param)):
        print('\t%s' %(param[i]))
        # initiate fig
        plt.clf()
        # print the 1D-pdf
        plt.plot(data1D.T[2*i], data1D.T[2*i+1], "mo", label = 'NA results')

        # find number of peaks/gaussian
        indexes = peakutils.indexes(data1D.T[2*i+1], 
                                    thres = 0.05, # Normalized threshold. Only the peaks with amplitude higher than the threshold will be detected.
                                    min_dist = 30)
                                    #min_dist = (max(data1D.T[2*i])-min(data1D.T[2*i]))/2) # Minimum distance between each detected peak.
                                    #min_dist = max(data1D.T[2*i])) # Minimum distance between each detected peak.
        #pplot(data1D.T[2*i], data1D.T[2*i+1], indexes)
        print('\t\tNumber of peaks : %s' %(indexes.shape))

        # Fit the Gaussian data
        # Compute the mean and 1 sigma erro used for the first guess of the fit
        mean = sum(data1D.T[2*i] * data1D.T[2*i+1]) / sum(data1D.T[2*i+1])
        sigma = np.sqrt(sum(data1D.T[2*i+1] * (data1D.T[2*i] - mean) ** 2) / sum(data1D.T[2*i+1]))
        # Do the fitting
        pars, cov = curve_fit(f = gaussians,
                            xdata = data1D.T[2*i], ydata = data1D.T[2*i+1], 
                            p0 = np.array([np.stack([data1D.T[2*i+1][indexes[k]], data1D.T[2*i][indexes[k]], sigma]) 
                                          for k in range(0, indexes.shape[0])]).reshape(-1),
                            bounds=(-np.inf, np.inf))
        # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
        stdevs = np.sqrt(np.diag(cov))
        # Calculate the residuals
        res = data1D.T[2*i+1] - gaussians(data1D.T[2*i], *pars)

        for k in range(0,indexes.shape[0]):
            # Compute the result for each peak
            pars_1 = pars[(3*k) : (3*k+3)]
            # print on screen results
            print('\t\tMean & sigma (Gaussian fit %s, %0.2f %s) : %0.2f +/- %0.2f' %(str(k+1), pars_1[0]*100, chr(37), pars_1[1], pars_1[2]))
            # Save the results in text file for each peak
            line = (str(param[i]) + 'peak ' + str(k+1) + '\t' + str(pars[1]) + '\t' + str(stdevs[1]) +
                                   '\t' + str(pars[2]) + '\t' + str(stdevs[2]) +
                                   '\n')
            fstats_w.write(line)
            # Plot the gaussian fitted
            plotgfit(data1D = data1D, i = i, pars = pars_1, ipeak = k+1)

        # Set axis names
        plt.xlabel(param[i], fontsize=size_x)
        plt.ylabel('Probability', fontsize=size_y)
        plt.ylim = 0
        plt.legend(loc = 'best')

        # Saving the plots as a pdf file
        plt.savefig(graph_path + 'PDF-1D_param' + str(i+1) + '.pdf')
        print('\t\tPlotting PDF: '+ graph_path + 'PDF-1D_param' + str(i+1) + '.pdf\n')

    # close the output text file
    fstats_w.close()

    return


##########################################################################################
##########################################################################################
def NAplot_Pecube(param, dataplot, tick_space,
                  inv_results = 'NA/NA_results.csv', data_nab = 'NA/NAB/nab.out',
                  graph_path  = 'NA/Graphs/',
                  PDF_1D = True, PDF_2D = False, pdf1d_results = None, pdf2d_results = None,
                  size_x = 15, size_y = 15, size_m = 15, size_plo = 50, size_mis = 50):
    
    """
    Main code of NA Plot for Pecube inversions with NA

    INPUTS:
        param (list of str): Define as many variable as you have, with their unit
                             Check the order in na.sum (open with text editor) or in NA_Results
                             If you use the later, first column always is the misfit
                             The others paramaters are from topo_parameters.txt then from fault_parameters.txt
                             Exemple : param=['Offset (km)','Basal Temperature (°C)','Slip rate (km∕Ma)']
                             /!\ IF YOU WANT TO USE A SLASH USE THIS ONE --> '∕' <-- , IT'S A UNICODE DIVISION SYMBOL
                            WINDOWS AND OSX DON'T ALLOW THE USE OF THE REGULAR SLASH 
        dataplot (list of couple of integers): Set the couple of variables to plot against each other
                                               Exemple : Offset vs Slip rates
                                               dataplot = [(1,3)] -->  plot=(1,3)
                                               If you plot 2D pdfs (contours), please, CHECK that the couple of parameters to plot
                                               are the same and in the same order than in the nab.in file.
                                               This python script checks it and will insult you if this is not compatible !!!
        tick_space (array of floats): Set the space between ticks for x and y axes for each parameters
                                      (same order than the list param).
                                      If the tick format does not fit your variables, 
                                      you may need to modify the dictionnary tick_order in the function multiplot
        inv_results (str, optional): Name of the NA file with the inversion results; Usually NA_Results.
                                     Defaults to 'NA/NA_results.csv'.
        data_nab (str, optional): Name of the NAB file with the inversion results. Usually nab.out.
                                  Defaults to 'NA/NAB/nab.out'.
        graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                     Do not forget the '/' at the end.
                                     Defaults to 'NA/Graphs/'.
        PDF_1D (bool, optional): Choose if you want the 1-PDFs (Probability Density Function) 
                                 Defaults to True.
        PDF_2D (bool, optional): Choose if you want the 2-PDFs (Probability Density Function) 
                                 Defaults to False.
        pdf1d_results (string, optional): Print the 1-pdfs in a text file.
                                          Defaults to None.
        pdf2d_results (string, optional): Print the 2-pdfs in a text file.
                                          Defaults to None.
        size_x (int, optional): Size of the font for the x axes label. 
                                Defaults to 15.
        size_y (int, optional): Size of the font for the y axes label. 
                                Defaults to 15.
        size_m (int, optional): Size of the font for the markers label. 
                                Defaults to 15.
        size_plo (int, optional): Size of the markers of the scatter plot. 
                                    Defaults to 50.
        size_mis (int, optional): Size of the markers of the misfits. 
                                    Defaults to 50.
    
    Raises:
        NameError: Input file names does not exists
        ValueError: Problems with number of values or with NaN
        ImportError: File format not supported, please, change the format
    """

    # If you want to print the pdfs in a text file, just modify the 2 next lines
    pdf1d_results = None #'NA/NAB/PDF_DATA.txt'
    pdf2d_results = None #'NA/NAB/PDF_2D_DATA.txt'

    print('###########################################################################\n')
    print('\t\tPlot results from NA inversions\n')
    print('\t\t  \xa9Xavier Robert - IRD-ISTerre\n')
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
        
    if not PDF_1D: 
        pdf1d_results = None
        data1D = None
    if not PDF_2D: 
        pdf2d_results = None
        data2D = None

    print('Building NAB files for plotting...')
    # Build Nab files for plotting
    if PDF_1D or PDF_2D:
        data1D, data2D = read_NABout(data_nab, dataplot, pdf1d_results, pdf2d_results)

    # Plot data
    # Check if there is a folder to store the outputs
    print('\nPlotting results...')
    if not os.path.exists(graph_path): os.mkdir(graph_path)
    if len(dataplot) > 1:
        for i in range(len(dataplot)):
            # Call the plot function for each parameter value to plot
            multiplot(param, nb_var, dataplot[i], 
                      tick_space, 
                      data1D, data2D,
                      size_x, size_y, size_m, size_plo, size_mis, 
                      PDF_1D, PDF_2D, graph_path, i)
    else:
        multiplot(param, nb_var, dataplot[0],
                  tick_space, 
                  data1D, data2D,
                  size_x, size_y, size_m, size_plo, size_mis, 
                  PDF_1D, PDF_2D, graph_path)
    
    # Compute the stats, print them
    print('Computing statistics...')
    if not PDF_1D:
        print('WARNING: No 1D-pdf computed --> No stats computed !')
    else:
        statsparam(data1D, param, graph_path, size_x, size_y)

    print('###########################################################################\n')
    print('\t\tEnd...\n')
    print('###########################################################################\n')

    return

##########################################################################################
##########################################################################################
if __name__ == u'__main__':
    
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
    inv_results = '../Tests/NA/NA_results.csv'
    data_nab    = '../Tests/NA/NAB/nab.out'
    graph_path  = '../Tests/NA/Graphs/'    # Do not forget the '/' at the end.

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

    # Set the size of the markers of the scatter plot and of the misfit
    size_plo = 50
    size_mis = 50

    # If you want to print the pdfs in a text file, just modify the 2 next lines
    pdf1d_results = None #'NA/NAB/PDF_DATA.txt'
    pdf2d_results = None #'NA/NAB/PDF_2D_DATA.txt'
    
    
    NAplot_Pecube(param = param,
                  dataplot = dataplot,
                  inv_results = inv_results,
                  data_nab = data_nab,
                  graph_path  = graph_path,
                  PDF_1D = PDF_1D, PDF_2D = PDF_2D,
                  pdf1d_results = pdf1d_results, pdf2d_results = pdf2d_results,
                  tick_space = tick_space,
                  size_x = size_x, size_y = size_y, size_m = size_m,
                  size_plo = size_plo, size_mis = size_mis
                  )
