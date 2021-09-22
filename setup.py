######!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
#from importlib.metadata import version

# Import of the lib pyRRIM
import pyPlotPecube

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='pyPlotPecube',
	version=pyPlotPecube.__version__,
	description='package that provides tools to plot forward and NA Pecube Results',
	#long_descritpion=open('README.rst','r').read(),
	url='https://github.com/robertxa/Pecube-Plot-Xav',
	download_url='https://github.com/robertxa/Pecube-Plot-Xav/archive/master.zip',
	author='Xavier Robert',
	author_email='xavier.robert@ird.fr',
	license='GPL-V3.0',
	packages=find_packages(),
	#include_package_data=True,	# What is the use of it ?
	install_requires=[
	      'scipy',
	      'peakutils',
	      'numpy',
	      'matplotlib'
	],
	classifiers=[
		"Operating System :: OS Independent"#,
		#"Topic :: Scientific/Engineering :: GIS"
	],
	include_package_data=True,
	zip_safe=False)
      