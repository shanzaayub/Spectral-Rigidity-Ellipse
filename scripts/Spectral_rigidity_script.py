""" This script generates the norm terms for each eccentricity using the files created in the previous steps,
based on the terms generated, we can infer dynamical spectral rigidity
"""

from scipy import special
from scipy import integrate
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import time

maxq = 500 # initialize max period number

semi_axes = pd.read_csv("e_and_semi_axes.txt", sep = '\t') # the file that contains the semi-axes associated with each eccentricity
semi_axes.drop("Unnamed: 0", axis = 1, inplace = True)

def collision_period(q):
    amplitudes = collision_amplitude(q)
    a = semi_axes[str_e][0]
    b = semi_axes[str_e][1]
    result = []
    for j in range(0,q):
        result.append([a*math.sin(amplitudes[j]),-b*math.cos(amplitudes[j])])
    return (result)

def collision_amplitude(q):
    """
    Returns an array of the q collision points relative to period q
    """

    result = []
    points = collision_pts[str(q).zfill(2)].dropna()
    for j in range(0,q):
        result.append(float(points[j]))
    return (result)

################### The following functions find the angle associated with each collision point #################################


def find_vector1(u,v): # u is the earlier point, v is the next point/the point we are looking at
    v1 = np.subtract(v,u)
    mag_v1 = math.sqrt((v1[0]**2)+(v1[1]**2))
    v1_unit = v1 / mag_v1
    return (v1_unit)

def find_tangent_vector(x,y,str_e):
    a = semi_axes[str_e][0]
    b = semi_axes[str_e][1]
    y_prime = (b*x)/a
    x_prime = -(a*y)/b
    mag = np.sqrt(x_prime**2+y_prime**2)
    y_tan = y_prime/mag
    x_tan = x_prime/mag
    return (x_tan, y_tan)

def sinphi_lst(q, str_e):
    '''Returns a list of sin(phi) for all point in a given periodic orbit'''
    sinphi_list = []
    for pt in collision_period(q):
        v1_unit = tuple(find_vector1(pt, collision_period(q)[(collision_period(q).index(pt)-1)%q]))
        v2_unit = find_tangent_vector(pt[0],pt[1],str_e)
        cosphi = np.dot(v1_unit,v2_unit)
        sinphi_list.append(np.sqrt(1-cosphi**2))
    return(sinphi_list)

################## The following functions find the lazutkin coordinate associated with each collision point ###########################

def lazutkin_coordinate_analytic(amplitude,e):
    return 0.25*(special.ellipkinc(amplitude,e**2)/special.ellipk(e**2)-1)

def mu_analytic(amplitude,e):
    return 2*special.ellipk(e**2)*np.sqrt((1-e**2)/(1-e**2*math.sin(amplitude)**2))


############# The following functions are needed for our inference of dynamical spectral rigidity ###########################

def T_of_q_j(q,j,str_e,e):
    '''T_j_q is a matrix, by varying q and j, you obtain an entry for that matrx.
    Think of the q as the rows, while the j are the columns. Fixing q and j gives you one entry.
    Fixing q and varying j will give you a row, fixing j and varying q will give you a column.
    Make sure to keep this in mind when you use this function, you will have to make changes to the code accordingly
    before you run this function.'''
    k_sum = []
    if q==1:
        mu_k = mu_analytic(np.pi/2,e)
        return 1./mu_k
    col_pt = collision_period(q)
    col_amp = collision_amplitude(q)
    sinphi_list_q = sinphi_lst(q,str_e)

    for k in range(q):
        sinphi = sinphi_list_q[k]

        laz_k = lazutkin_coordinate_analytic(col_amp[k],e)

        mu_k = mu_analytic(col_amp[k],e)

        sum_exp = sinphi*(np.cos(2*np.pi*j*laz_k)/mu_k)
        k_sum.append(sum_exp)
    return(sum(k_sum))

def lambda_marvizi_melrose(j,str_e,e):
    """
    This computes an approximate limit for q→∞ of T_qj
    it computes the magic_q'th element of the vector

    """
    if (j%2==1): # Ellipse has additional symmetry: every odd term is 0
        return 0
    magic_q=min(j+10,999)
    q=j;
    mmc=(q**2*T_of_q_j(q,j,str_e,e));
    continuing=True
    accord=0.000001
    while continuing:
        q=q+1;
        mmc_new=(q**2*T_of_q_j(q,j,str_e,e));
        continuing=(np.abs(mmc-mmc_new)<accord)
    print ("mmc:",j," ",q," ",magic_q,"   -   ",mmc_new)
    return mmc_new
# We probably need to justify the choice of magic_q


def norm_partial_q(gamma,q,str_e,e):
    terms=[]
    for j in np.arange(1,q * arbitrary_accuracy):
        delta = 1 if q==j else 0
        terms.append(j**(-gamma)*abs(T_of_q_j(q,j,str_e,e)-delta-lambda_MM[j-1]/(q**2)))
    return (q**gamma)*sum(terms)

########## Some considerations in the interest of computational time ###########################

magic_j=750 #arbitrary choice, consider this to be the point when the terms are close enough to their "actual" values
gamma=3.5 #for faster decay, can change this to get smaller values but keep in mind the decay may be slower


norm_terms_circle=[] # these are the terms we get for the case of the circle (eccentricity 0) to compare to our values
for q in np.arange(1,magic_j):
    cq = np.sinc(1/q) if q>1 else 1/np.pi
    norm_terms_circle.append(1+cq*(special.zeta(gamma)-2))
arbitrary_accuracy=100 # ideally this should give an error of at most 1/accuracy²

######################### Putting it all together for the computation ###################################

sampled_e=[]
norm_e=[]
norm_dict = {}
for e in np.arange(0,1,0.1): # the eccentricities to consider, change this to look at other ones
    t0 = time.time()
    str_e = '%.2f'%e
    sampled_e.append(e)
    collision_pts = pd.read_csv("all_periods_{}e_col_amplitudes.txt".format(str_e), sep = '\t') # the file that contains collision points for every period given an eccentricity
    collision_pts.drop("Unnamed: 0", axis = 1, inplace = True)
    lambda_MM=[]
    print("\n\n---\neccentricity",str_e)
    print("Collision points loaded")

    # caching the Marvizi-Melrose coefficients up to magic_j
    for j in np.arange(1,magic_j):
        l = lambda_marvizi_melrose(j,str_e,e)
        if (j%2 == 0) and (abs(l) < 0.0000001): #stopping the computation after the coefficients for even j are close enough to 0
            break
        lambda_MM.append(l)
    for j in np.arange(len(lambda_MM)+1,maxq*arbitrary_accuracy): #at this point the terms were close enough to 0, hence no need to compute
        lambda_MM.append(0)
    ##
    print("Marvizi-Melrose coefficents cached")

    norm_terms=[]
    for q in np.arange(1,magic_j):
        print("Generating norm term for q =",q)
        norm_term = norm_partial_q(gamma,q,str_e,e)
        norm_terms.append(norm_term)
        print("Norm term:",norm_terms[q-1], "(reference for the circle:", norm_terms_circle[q-1], ")")
        accord=abs((norm_terms[q-1]-norm_terms_circle[q-1])/norm_terms_circle[q-1])

        if (q > 3) and (accord < .5): #difference between the actual value and the approximate value is less than 5%
            break
        elif (q > 20) and (norm_term < 0.5): #the terms are less than 1
            break
        elif (q > 30) and (norm_term > 1): #high period but still not less than 1, meaning the information generated is not reliable
            print("This eccentricity {} was giving unreliable information".format(str_e))
            break
    temp_lst = []
    for i in norm_terms:
        temp_lst.append(abs(i))
    norm_e.append(max(temp_lst))
    norm_dict[str_e] = max(temp_lst)
    t1 = time.time()
    print("done finding norms for",str_e,":",norm_e[-1], "\ntime taken:", t1-t0)

series = pd.Series(norm_dict)
df = pd.DataFrame({'Eccentricity':series.index, 'Norm term':series.values})
step = float(df['Eccentricity'][1]) - float(df['Eccentricity'][0])
df.to_csv("Norms_till_{}".format(str_e)+"_{}step_e".format(str(step))+"_{}gamma.txt".format(str(gamma)), sep = '\t') # creating the file that contains the max norm term reached for eccentricity considered
print("All done!")
