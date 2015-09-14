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

import numpy as np, scipy.optimize as opt, scipy.integrate as integ
import parabolic_profile as pp

# Useful constants in IS system
# ------ --------- -- -- ------

a = 0.18; Rwall = 0.22;
epslon_o = 8.854187817E-12; e = 1.6021766E-19; 
c = 299792458; me = 9.109389E-31

# Functions
# ---------

def ChangedSignal(a, b):
    if   (a < 0. and b > 0.): return True
    elif (a > 0. and b < 0.): return True
    else: return False;

def Profile (x, b, s1, s2) :
    """ return (1 - (x/a)**2)**alpha + beta * np.exp(-(x / s)**2); """
    if (abs(x) >= a): return 0;
    if (s1 == 0.): return (b**2) * np.exp(-(x/s2)**2);
    if (s2 == 0.): return np.exp(-(x/s1)**2);
    return np.exp(-(x/s1)**2) + (b**2) * np.exp(-(x/s2)**2);

def Opt_X_root(x, b, s1, s2, frac):
    return Profile(x, b, s1, s2) - frac;

def FrequencyToDensity (F) :
    """ Convert plasma frequency to equivalent density """
    return 4 * (np.pi**2) * (F**2) * me * epslon_o / (e**2)

def FrequencyToX (frac, b, s1, s2) :
    """ Convert plasma frequency to equivalent critic position
    return a * np.sqrt( 1 - (FrequencyToDensity(F)/n0)^(1/alpha) ); """
    f_inf = Opt_X_root(0, b, s1, s2, frac);
    f_sup = Opt_X_root(a, b, s1, s2, frac);
    if (not ChangedSignal(f_sup, f_inf)): return a;
    return opt.brentq(Opt_X_root, 0, a, args=(b, s1, s2, frac));

def Func(x, b, s1, s2, r):
    return 1.0 / np.sqrt(1.0 - Profile(x, b, s1, s2) / r );

def ProfileInt (xc, b, s1, s2, r) :
    """ Make numeric integration based on Density profile with trapezes
        r: ratio of critic density to central density.              """
    #if (xc > (0.995) * a): return 0;
    if (xc > 0):
        return integ.quad(Func, xc, a, args=(b, s1, s2, r),
                epsrel=1.0e-3)[0];
    else:
        return integ.quad(Func, -a, a, args=(b, s1, s2, r), 
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

def OptGroupDelay (F, n0, b, s1, s2) :
    """ Function to be optimized - alpha Optional parameter."""

    r = FrequencyToDensity(F) / n0;
    groupDelay = np.ones([r.size]) * 1E-9;
    groupDelay = groupDelay.reshape(r.shape);
    for i in range(r.size) :
        if (r[i] < 1 + b**2) :
            xc = abs(FrequencyToX(r[i], b, s1, s2));
            groupDelay[i] = 2 * ProfileInt(xc, b, s1, s2, r[i]) / c;
        else :
            groupDelay[i] = 2 * ProfileInt(0, b, s1, s2, r[i]) / c \
                    + 2 * (Rwall - a)/c;
    return groupDelay;
