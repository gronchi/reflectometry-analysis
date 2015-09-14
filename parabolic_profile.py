# -*- coding: utf-8 -*-
"""
    Module for approximated parabolic density profile
    ------ --- ------------ --------- ------- -------
    
    Dependencies : Numpy.py
    
    Developed by : Alex V. Andriati - USP.
    
    Main functions of interesting to optimize curve fitting by numeric
    integration of parabolic density profile to Group Delay.
    
    Numeric integration with trapezium method. More Datils see ProfileInt
    function.
"""

import numpy as np, scipy.integrate as integ

# Useful constants in IS system
# ------ --------- -- -- ------

a = 0.18; Rwall = 0.22;
epslon_o = 8.854187817E-12; e = 1.6021766E-19; 
c = 299792458.0; me = 9.109389E-31;

# Functions
# ---------

def Profile (x, alpha) :
    """ return (1 - (x/a)^2)^alpha """
    if (abs(x) > a) : return 0;
    return (1 - (x/a)**2)**alpha

def FrequencyToDensity (F) :
    """ Convert plasma frequency to equivalent density """
    return 4 * (np.pi**2) * (F**2) * me * epslon_o / (e**2)

def FrequencyToX (frac, alpha) :
    """ Convert plasma frequency to equivalent critic position
    return a * np.sqrt( 1 - (FrequencyToDensity(F)/n0)^(1/alpha) ); """
    return a * np.sqrt( 1.0 - frac**(1.0/alpha) );

def Func(x, alpha, r):
    return 1.0 / np.sqrt(1.0 - Profile(x, alpha) / r );

def ProfileInt (xc, alpha, r) :
    """ Make numeric integration based on Density profile with trapezes
        r: ratio of critic density to central density.              """
    if (xc > 0):
        return integ.quad(Func, xc, a, args=(alpha, r),
                epsrel=1.0e-3)[0];
    else:
        return integ.quad(Func, -a, a, args=(alpha, r), 
                epsrel=1.0e-3)[0];

""" 
    Note that for a data set of frequencies in attempt to fitting a curve
    many possible values of center density is tried and all parameters is
    treated as array like.
"""
""" 
    Split the compute of profile integration in two cases:
    1. r <  1 ---> Critical density smaller than Center Density
                   Integrate until equivalent X reflection point.
    2. r >= 1 ---> Critical density greater than Center Density
                   Integrate over plasma diameter.
"""

def OptGroupDelay (F, n0, alpha) :
    """ Function to be optimized - alpha Optional parameter."""
    if (alpha <= 0) :
        groupDelay = np.zeros(F.size);
        groupDelay = groupDelay.reshape(F.shape);
        return groupDelay;
	
    r = FrequencyToDensity(F) / n0;
    groupDelay = np.ones([r.size]) * 1E-9;
    groupDelay = groupDelay.reshape(r.shape);
    for i in range(r.size) : 
        if (r[i] < 1) :
            xc = FrequencyToX(r[i], alpha);
            groupDelay[i] = 2 * ProfileInt(xc, alpha, r[i]) / c;
        else :
            groupDelay[i] = 2 * ProfileInt(False, alpha, r[i]) / c \
            + 2 * (Rwall - a)/c;
    return groupDelay;

def Residues (xdata, ydata, n0, alpha) :
    """ return OptGroupDelay(xdata, n0, alpha) - ydata; """
    return ydata - OptGroupDelay(xdata, n0, alpha);
