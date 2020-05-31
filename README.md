# Spectral-rigidity-ellipse

We use the scripts contained herein to provide numerical evidence for dynamical spectral rigidity among ellipses of
various eccentricities. The code contained in these scripts follows the technique implemented in a paper by J. De Simoi, V. Kaloshin and Q. Wei (a link to the paper is here: [Dynamical_Spectral_Rigidity](https://annals.math.princeton.edu/2017/186-1/p07)).  
There are 3 scripts:
1. `e_and_semi_axes_file.py`
2. `Col_pts_find.py`
3. `Spectral_rigidity_script.py`

### e_and_semi_axes_file.py
_Objective: Finds the semi-axes associated with a given eccentricity._  
Fixes the circumference of an ellipse to be 1, and the eccentricity, `e` in the interval \[0,1) using a step-size of 0.01. The script finds the semi-major axis, denoted `a`, using the Complete Elliptic Integral of the Second kind, and then finds the semi-minor axis, denoted `b`, using the formula for the eccentrcity of an ellipse.  
Generates one file:
* e_and_semi_axes.txt 

### Col_pts_find.py
_Objective: Finds the collision points of orbits of periods 1 to 500, for each eccentricity._  
First finds a sequence of values \lambda_q corresponding to periodic orbits of rotation number 1/q (where q is in 1 to 500). This is done by numerically inverting the formula for rotation number, \omega_\lambda of the orbit associated to the caustic C_\lambda, using the bisecting method. Then, finds the collision points.  
Generates one file:
* all_periods_`<eccentricity>`e_col_amplitudes.txt  

### Spectral_rigidity_script.py
_Objective: Finds the norm terms for each eccentricity to infer dynamical spectral rigidity for the associated ellipse._  
Uses the files generated from the earlier scripts, and applies the technique implemented in the paper by J. De Simoi, V. Kaloshin and Q. Wei to determine spectral rigidity.  
Generates one file:  
* Norms_till_`<max_eccentricity_considered>`\_`<step_size>`step_e\_`<gamma_used>`gamma.txt  

*Note: This file will be named differently based on the max eccentricity, step-size between the eccentricities considered and the value of \gamma used when the script is run*.
