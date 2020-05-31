""" This script creates the file e_and_semi_axes.txt, which contains,
for each eccentricity e, the semi-major and semi-minor
axes for an ellipse with eccentricity e, and circumference 1
"""
from scipy import special #library for elliptical integrals
import numpy as np #library for math manipulations and functions
import math #another library for math manipulations and functions
import pandas as pd # for dataframe manipulations and creation, mostly to organize everything

semi_axis_dict = {}

for e in np.arange(0,1,0.01):
    a = 1/(4*special.ellipe(e**2)) # e**2 is our modulus
    b = a*math.sqrt((1-(e**2)))
    semi_axis_dict['%.2f'% e] = [a,b]

pd.DataFrame(semi_axis_dict).to_csv("e_and_semi_axes.txt", sep = '\t') #make the semi axis dict into a file
