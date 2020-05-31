"""This script creates files containing the collision points for
period q in [2,maxq) for eccentricities in [0,1)

"""

from scipy import special #library for elliptical integrals
import numpy as np #library for math manipulations and functions
import math #another library for math manipulations and functions
import pandas as pd # for dataframe manipulations and creation, mostly to organize everything

maxq=500

semi_axes = pd.read_csv("e_and_semi_axes.txt", sep = '\t') # file containing the semi-axes associated with each eccentricity
semi_axes.drop("Unnamed: 0", axis = 1, inplace = True)

def rotation_no(l,e):
    a = semi_axes[e][0]
    b = semi_axes[e][1]
    k_l_sq = (a**2 - b**2)/(a**2 - l**2)
    return special.ellipkinc(np.arcsin(l/b),k_l_sq)/(2*special.ellipk(k_l_sq))

def find_lambda(w):
    """ This function finds an approximate lambda associated with a rotation number w
    """
    l_eccen_dict = {}
    for e in semi_axes:
        a = semi_axes[e][0]
        b = semi_axes[e][1]
        start = 0
        end = b
        l = (start + end)/2
        while True:
            w_0 = rotation_no(l,e)
            if abs(w-w_0) < 0.0000000001:
                l_eccen_dict[e] = l
                break
            elif w > w_0:
                start = l
                l = (start + end)/2
            else:
                end = l
                l = (start + end)/2
    return(l_eccen_dict)


period_lambda_dict = {} # make the dictionary and find the lambdas for orbits of different periods
for q in np.arange(3,maxq): # these are the periods, q =1 and q=2 are treated separately in subsequent script
   # w = 1/q
    period_lambda_dict[str(q)] = find_lambda(1/q)
    print(q,"\r");
print("done finding λ")
eccen_col_dict = {}

def find_collision_pts(e,q):
    """ Finds the collision points for a period q, given an eccentricity e
    """
    collisions_dict = {}
    a = semi_axes[e][0]
    b = semi_axes[e][1]
    l = period_lambda_dict[q][e]
    k_l_sq = ((a**2)-(b**2))/((a**2)-(l**2))
    for j in range(int(q)):
        d_l_q = (4*(special.ellipk(k_l_sq)))/int(q)
        t_j = (special.ellipk(k_l_sq))+j*d_l_q
        collisions_dict[str(j).zfill(2)] = special.ellipj(t_j,k_l_sq)[3]
    return (collisions_dict)

for e in semi_axes:
    eccen_row_dict ={}
    # add by hand the bouncing ball orbit
    eccen_row_dict["02"] = {"00" : math.pi/2,
                            "01" : 3*math.pi/2}
    for q in np.arange(3,maxq):
        eccen_row_dict[str(q).zfill(2)] = find_collision_pts(e,str(q)) #the zfill() makes sure there are only 2 significant figures
        print (q,"\r");
    df = pd.DataFrame(eccen_row_dict)
    df.to_csv("./all_periods_{}e_col_amplitudes.txt".format(e),sep='\t') # creating the file containing the collision pts for every period given an eccentricity

print("Done finding collision points")
