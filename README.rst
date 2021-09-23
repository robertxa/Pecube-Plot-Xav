Python scripts to plot PECUBE modeling results
==============================================

Xavier Robert, Grenoble,                         

This tool have been set up to help the plot of PECUBE (Braun, 2003; Braun et al., 2012) modeling results.

    - The Plot-PecubeV4+-Forward.py permits to plot predictions versus observations for Pecube V4+ forward modeling.

    - The Inversion_plots-XR.py permits to plot the results of NA Inversions (Sambridge, 1999 a & b) with scatter plots and 1D / 2D marginals. To plot the marginals, the results of NA Inversions should be analyzed with NAB.

Usage
-----
	
After copying the scripts in your folder ant editing the main sections to fit with your project, each script can be ran in a terminal with the command : 

.. code-block:: bash

	~$ python pyPlotPecubeForward.py

or:

.. code-block:: bash

    ~$ python pyPlotPecubeNA.py

Please, edit the Python files in your text-editor and read the headers of the main section in each files to know how to parametrize them. 

The second solution to install the module with pip:

.. code-block:: bash

    ~$ pip install pyPlotPecube

And then, inside a Python environement, import the function you need:

.. code-block:: python

    >>> from pyPlotPecube import pyPlotPecubeForward as Fplot

    >>> from pyPlotPecube import pyPlotPecubeNA as NAplot

and finally run the function in your Python environement :

.. code-block:: python

    >>> Fplot.PlotPecubeForward(dataplot = dataplot, graphpath = graphpath, datafnme = datafnme, inputdata = inputdata, graphtitle = graphtitle, agerange = agerange, profiletype = profiletype, size_x = size_x, size_y = size_y, A = A, B = B)

    >>> NAplot.NAplot_Pecube(param = param, dataplot = dataplot, inv_results = inv_results, data_nab = data_nab, graph_path  = graph_path, PDF_1D = PDF_1D, PDF_2D = PDF_2D, pdf1d_results = pdf1d_results, pdf2d_results = pdf2d_results, tick_space = tick_space, size_x = size_x, size_y = size_y, size_m = size_m, size_mis = size_mis)

Parameters
----------

For Forward plots, arguments are :

    1. ``datafnme`` (string): results of Pecube forward modeling.

        No default value.

	2. ``inputdata`` (string): input data declared in Pecube.in. This is used to plot the errorbars.

        No default value.

	3. ``dataplot`` (list, optional): List of data to plot: ``['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']``; by default, the altitude will be plotted; Do not forget the simple quotes !!! 
        
        Note: For the moment, MTL and TTp not implemented.
        
        Defaults = ``['AHe','AFT']``.

	4. ``graphpath`` (str, optional): name of the folder where the plot will be written. Usually you do not have to change it.
        
        Defaults = ``'Graphs'``.

	5. ``graphtitle`` (str, optional): title to write on the graph. 
        
        Defaults = ``None``.

	6. ``agerange`` (2*1 array of floats, optional): range of the ages to plot on the profiles ``[min, max]``.
        
        Defaults = ``None``.

	7. ``profiletype`` (list, optional): type of profile, could be one or more of ``['Latitude', 'Longitude', 'Altitude', 'Projected']``. If ``[]``, no age profile is plotted.
        
        Defaults = ``[]``.

    7. ``size_x`` (int, optional): Font size for x-axis. 
    
        Defaults = ``15``.
        
	8. size_y (int, optional): Font size for y-axis. 
        
        Defaults = ``15``.

	9. ``A``, ``B`` (floats, optional): If need of a projected transect, define the line along which we will project with the coordinate of the point A and B defining, respectivelly, the begining and the end of the transect, in lat-long/WGS84.
        
        Defaults = ``None``.
    
    10. ``agename`` (dict, optional): legend of each data system         
		
        If ``None``, this is set to
        
            {'AHe' : 'AHe (Ma)',
		
            'AFT' : 'AFT (Ma)',
		
            'ZHe' : 'ZHe (Ma)',
			
            'ZFT' : 'ZFT (Ma)',
			
            'KAr' : 'KAr (Ma)',
			
            'BAr' : 'Biot. Ar (Ma)',
			
            'MAr' : 'Musc. Ar (Ma)',
			
            'HbAr' : 'Hb Ar (Ma)',
			
            'FTL' : 'FT length (µm)'}
			
        Default = ``None``. 
	
    11. ``predname`` (dict, optional): legend of each predicted system.
        
        If ``None``, this is set to:
            
            {'AHe' : 'Predicted AHe (Ma)',
            
            'AFT' : 'Predicted AFT (Ma)',
            
            'ZHe' : 'Predicted ZHe (Ma)',
            
            'ZFT' : 'Predicted ZFT (Ma)',
            
            'KAr' : 'Predicted KAr (Ma)',
            
            'BAr' : 'Predicted Biot. Ar (Ma)',
            
            'MAr' : 'Predicted Musc. Ar (Ma)',
            
            'HbAr' : 'Predicted Hb Ar (Ma)',
            
            'FTL' : 'Predicted FT length (µm)'}         
        
        Default = ``None``. 
	
    12. ``colores`` (dict, optional): Colors used for the different age system
    
        If ``None``, this is set to:
            
            {'AHe' : 'y',
            
            'AFT' : 'r',
            
            'ZHe' : 'g',
            
            'ZFT' : 'b',
            
            'KAr' : 'k',
            
            'BAr' : 'c',
            
            'MAr' : 'm',
            
            'HbAr' : '0.75',
            
            'FTL' : 'y'}
        
        Default = ``None``.


For Inverse plots, arguments are :

    1. ``param`` (list of str): Define as many variable as you have, with their unit. Check the order in ``na.sum`` (open it with a text editor) or in ``NA_Results``. If you use the later, first column is always the misfit. /!\ If you want to use a ``slash`,  use this one --> ``'∕'`` <-- . The regular slash it is a unicode division symbol windows and OSX do not allow the use of it.
        
        Exemple : ``param = ['Offset (km)','Basal Temperature (°C)','Slip rate (km∕Ma)']``
                             
    
    2. ``dataplot`` (list of couple of integers): Set the couple of variables to plot against each other. If you plot 2D pdfs (contours), please, CHECK that the couple of parameters to plot are the same and in the same order than in the nab.in file. This python script checks it and will insult you if this is not compatible !!!
        
        Exemple : Offset vs Slip rates, ``dataplot = [(1,3)]`` -->  plot=(1,3); 
        
        No default value.
    
    3. ``tick_space`` (array of floats): Set the space between ticks for x and y axes for each parameters (same order than the list param). If the tick format does not fit your variables, you may need to modify the dictionnary ``tick_order`` in the function multiplot.

        No default value.
    
    4. ``inv_results`` (str, optional): Name of the NA file with the inversion results, usually ``NA_Results``.
        
        Defaults = ``'NA/NA_results.csv'``.
    
    5. ``data_nab`` (str, optional): Name of the NAB file with the inversion results, usually ``nab.out``. 
        
        Defaults = ``'NA/NAB/nab.out'``.
    
    6. ``graph_path`` (str, optional) : Path where to save graphs and results Usually NA/Graphs. 
        
        Do not forget the ``'/'`` at the end. 
        
        Defaults = ``'NA/Graphs/'``.
    
    7. ``PDF_1D`` (bool, optional): Choose if you want the 1-PDFs (Probability Density Function); 
        
        Defaults to ``True``.
    
    8. ``PDF_2D`` (bool, optional): Choose if you want the 2-PDFs (Probability Density Function); 
        
        Defaults = ``False``.
    
    9. ``pdf1d_results`` (string, optional): Print the 1-pdfs in a text file. 
        
        Defaults = ``None``.
    
    10. ``pdf2d_results`` (string, optional): Print the 2-pdfs in a text file. 

        Defaults = ``None``.
    
    11. ``size_x`` (int, optional): Size of the font for the x axes label. 

        Defaults = ``15``.
    
    12. ``size_y`` (int, optional): Size of the font for the y axes label. 
        
        Defaults = ``15``.
    
    12. ``size_m`` (int, optional): Size of the font for the markers label. 
    
        Defaults = ``15``.
    
    13. ``size_mis`` (int, optional): Size of the markers of the misfits. 
        
        Defaults = ``50``.
    
    14. ``peak_thres`` (float, optional): Threshold to find peaks; between 0. and 1. See peakutils documentation.
                                      
        Default = 0.05.

    15. peak_min_dist (interger, optional): Minimum distance between the peaks. See peakutils documentation.
        
        Default = 30.

Examples
--------

Two example's sripts are in the ``Tests/`` folder. Just run them from their location :

.. code-block:: bash
    
    ~$ python Template_ForwardPlot.py
    
    ~$ python Template_NAPlot.py

For instance, the Template_ForwardPlot.py script gives the plot:

.. image:: https://github.com/robertxa/Pecube-Plot-Xav/tree/master/Tests/Graphs/Forward.png
   :scale: 50 %
   :align: center
   :alt: Data and predictions along a projected transect.


and the Template_NAPlot.py script permits to build the plots:

.. image:: https://github.com/robertxa/Pecube-Plot-Xav/tree/master/Tests/Graphs/NA-1.png
   :scale: 50 %
   :align: center
   :alt: NA inversion results for 2 parameters; on the scatter plot, each point corresponds to a model, and the color corresponds to the value of the misfit for that model; The curves on the sides shows the 1D-pdf of each parameter.

.. image:: https://github.com/robertxa/Pecube-Plot-Xav/tree/master/Tests/Graphs/NA-2.png
   :scale: 50 %
   :align: center
   :alt: curve fitting of the 1D-pdf of one parameter after NA inversion.

How to cite
-----------

Please, if you use this module, cite :
``Robert X., pyPlotPecube, a python module to plot PECUBE forward and inverse modeling results (2021), DOI:10.5281/zenodo.5521061``

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5521061.svg
    :target: https://doi.org/10.5281/zenodo.5521061

Contact
-------

If needed, do not hesitate to add a new branch or to contact the author. 
Please, use `https://www.isterre.fr/identite_id135055.html# <https://www.isterre.fr/identite_id135055.html#>`_

Licence
-------

Copyright (c) 2021 Xavier Robert <xavier.robert@ird.fr>

This package is licenced with <SPDX-License-Identifier: GPL-3.0-or-later>

