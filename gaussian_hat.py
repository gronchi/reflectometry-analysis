# -*- coding: utf-8 -*-
"""
    Module for approximated parabolic density profile
    ------ --- ------------ --------- ------- -------
    
    Dependencies : numpy, scipy.optimize, scipy.integrate
    
    Developed by : Alex V. Andriati - USP.
    
    Reflectometry Analysis give us explicit the group delay of the waves
    for a given probe frequency choosen. This approach avoid any transform
    to get directly the density profile and try to fit a curve based on
    least Square method for a alredy declared profile.

    This module try to fit a parabolic profile plus a gaussian, wich means
    that introduce a 'hat' in the center.
"""

import numpy as np, scipy.optimize as opt, scipy.integrate as integ
import parabolic_profile as pp

# Useful constants
# ------ ---------

a = 0.18; Rwall = 0.225; # plasma chamber dimensions (TCABR).
epslon_o = 8.854187817E-12; e = 1.6021766E-19; 
c = 299792458; me = 9.109389E-31

# Functions
# ---------

def ChangedSignal(a, b):
    """ Return True if a * b < 0. """
    if   (a < 0. and b > 0.): return True
    elif (a > 0. and b < 0.): return True
    else: return False;

def Profile(x, alpha, A, s):
    """ return (1 - (x/a)^2)^alpha + A * exp(-(x / s)^2);
    Some conditions required as Profile(|x| > a) = 0 and
    s != 0 to avoid problem with exponential function. """
    if (abs(x) >= a): return 0;
    return (1 - (x/a)**2)**alpha + (A**2) * np.exp(-(x / s)**2);

def Opt_X_root(x, alpha, A, s, frac):
    """ To find numerically the cut off position. """
    return Profile(x, alpha, A, s) - frac;

def FrequencyToDensity(F):
    """ Convert plasma frequency to equivalent density """
    return 4 * (np.pi**2) * (F**2) * me * epslon_o / (e**2)

def FrequencyToX(frac, alpha, A, s):
    """ Convert plasma frequency to equivalent critic position
    that depends on a given profile. """
    f_inf = Opt_X_root(0, alpha, A, s, frac);
    f_sup = Opt_X_root(a, alpha, A, s, frac);
    if (not ChangedSignal(f_sup, f_inf)): return a;
    return opt.brentq(Opt_X_root, 0, a, args=(alpha, A, s, frac));

def Func(x, alpha, A, s, r):
    """ Function to numerically integrate for Group Delay. """
    return 1.0 / np.sqrt(1.0 - Profile(x, alpha, A, s) / r );

def ProfileInt (xc, alpha, A, s, r) :
    """ Make numeric integration based on Density profile with 
    gauss quadrature. If xc < 0 or boolean means that has no
    reflection point in plasma (nc > n_max). """
    if (xc > (0.995) * a): return 0; # tolerance.
    if (xc > 0):
        return integ.quad(Func, xc, a, args=(alpha, A, s, r),
                epsrel=1.0e-3)[0];
    else:
        return integ.quad(Func, -a, a, args=(alpha, A, s, r), 
                epsrel=1.0e-3)[0];

def OptGroupDelay (F, n0, alpha, A, s) :
    """ Function to be optimized, or to theorical curve of
    Group Delay based on the declared Profile above. 
    avoid alpha < 0 and if s = 0 return just parabolic form. """
    if (alpha <= 0) :
        groupDelay = np.zeros(F.size);
        groupDelay = groupDelay.reshape(F.shape);
        return groupDelay;
    if (s == 0): return pp.OptGroupDelay(F, n0, alpha);

    r = FrequencyToDensity(F) / n0;
    groupDelay = np.ones([r.size]) * 1E-9;
    groupDelay = groupDelay.reshape(r.shape);
    for i in range(r.size) :
        if (r[i] < 1 + A**2) :
            xc = abs(FrequencyToX(r[i], alpha, A, s));
            groupDelay[i] = 2 * ProfileInt(xc, alpha, A, s, r[i]) / c;
        else :
            groupDelay[i] = 2 * ProfileInt(False, alpha, A, s, r[i]) / c
            groupDelay[i] += 2 * (Rwall - a)/c; # add vaccum part.
    return groupDelay;

def Residues(xdata, ydata, n0, alpha, A, s):
    """ Return Residuals of each dat point from a given set of params. """
    return ydata - OptGroupDelay(xdata, n0, alpha, A, s);
