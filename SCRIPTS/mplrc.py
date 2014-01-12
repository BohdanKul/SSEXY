"""
mplrc.py

Adrian Del Maestro
07.27.2012

matplotlib rc params and axes rectangles to generate figures of appropriate
size for different types of publication.

see: http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
"""
from math import sqrt
fig_width_pt = 246.0                    # Get this from LaTeX using \showthe\columnwidth
fig_width_pt = 510.0
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5.0)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
aps = {'params': {'axes.labelsize': 10,
                  'text.fontsize': 15,
                  'legend.fontsize': 15,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'font.family': 'serif',
                  'font.serif': 'Computer Modern Roman',
                  'test.usetex': True,
                  'figure.figsize': fig_size,
                  'xtick.major.size': 4,
                  'xtick.minor.size': 2,
                  'xtick.major.pad': 4,
                  'xtick.minor.pad': 4,
                  'ytick.major.size': 4,
                  'ytick.minor.size': 2,
                  'ytick.major.pad': 4,
                  'ytick.minor.pad': 4,
                  'axes': [0.13,0.2,0.95-0.13,0.95-0.2]}}

fig_width_pt = 510.0
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5.0)+1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width/golden_mean*1.2      # height in inches
fig_size =  [fig_width,fig_height]
PRB = {'params': {'axes.labelsize': 10,
                  'text.fontsize': 15,
                  'legend.fontsize': 15,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'font.family': 'serif',
                  'font.serif': 'Computer Modern Roman',
                  'test.usetex': True,
                  'figure.figsize': fig_size,
                  'xtick.major.size': 4,
                  'xtick.minor.size': 2,
                  'xtick.major.pad': 4,
                  'xtick.minor.pad': 4,
                  'ytick.major.size': 4,
                  'ytick.minor.size': 2,
                  'ytick.major.pad': 4,
                  'ytick.minor.pad': 4,
                  'axes': [0.13,0.2,0.95-0.13,0.95-0.2]}}
