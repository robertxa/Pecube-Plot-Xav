######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2021 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later


from __future__ import  division
# This to be sure that the result of the division of 2 integers is a real, not an integer
from __future__ import absolute_import
from __future__ import print_function

__version__ = '1.2.0'

# Import modules
import sys, os, time
#import copy
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

# Import all the functions
from .pyPlotPecubeNA import *
from .pyPlotPecubeForward import *
