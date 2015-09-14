"""
    Python Module of Functions
    ------ ------ -- ---------

    Dependencies: gaussian_hat, matplotlib.pyplot, numpy, scipy.optimize.
    -------------

    Example to Use:
    ------- -- ----
    >>> import gaussian_hat as gh.
    >>> dict = gh.fit_GD(xdata, ydata, path_to_directory, time, show=boolean).
    where just first and second argument are obrigatory.
"""

import scipy.optimize as opt, gaussian_hat as gh;
import numpy as np, matplotlib.pyplot as plt

# Define Figure Custom Options.
from matplotlib import rcParams;
rcParams['font.family'] = 'sans-serif';
rcParams['font.sans-serif'] = ['Computer Modern Sans serif'];
rcParams['ytick.labelsize'] = 'large';
rcParams['xtick.labelsize'] = 'large';
rcParams['text.usetex'] = True;

a = 0.18; Rwall = 0.22; c = 299792458.0;

def CrossCenter(GD):
    """ Find most probably index that cross max density. """
    sentinel = -1;
    for i in range(GD.size) :
        if (GD[i] > 2.2e-9): return i;
    return sentinel;

def fit_GD(pf, gd, path='', time=70.0, show_it=False):
    """ First two arguments takes the data. Others arguments optional.
        Based on data points call curve_fit for the profile defined
        by gaussian_hat, wich can be generalized for others profiles.

        Return a python dictionary with keys:

        params -> Result params of the curve_fit
        mcov ---> Covariance Matrix
        res ----> Residuals for each point to the curve.
    """

    # Find the most possible 'cross center (CC)' density point.
    # Print a message to be patient. Can take more than a min.
    CC = CrossCenter(gd);
    if (CC > 0):
        n0_ref = gh.FrequencyToDensity(pf[CC]);
        print '\nAjustando curvas. Pode demorar... n0_ref = %.2g' % n0_ref
    else:
        n0_ref = gh.FrequencyToDensity(40.5E9);
        print '\nAjustando curvas. Pode demorar... n0_ref = none'

    # Take the best curve by least square method.
    # It depends strongly of initial guess in both quality and time.
    # Besides it, make some statistics of the result.
    guess = [0.87*n0_ref, 1.1, 0.45, a / 4.8];
    p, mcov = opt.curve_fit(gh.OptGroupDelay, pf[:CC], gd[:CC], p0=guess);
    residuals = gh.Residues(pf[:CC], gd[:CC], p[0], p[1], p[2], p[3]);
    # Remake the curve data for more resolution and change units.
    # Calculate some immediately results and standart deviation.
    n_max = p[0] * (1 + p[2]**2);
    Err   = np.std(residuals, ddof=1) / 1e-9;
    PF    = np.linspace(pf[0] - 5E9, pf[-1] + 5E9, 1e3);
    GD    = gh.OptGroupDelay(PF, p[0], p[1], p[2], p[3]);
    # Change units
    pf = pf / 1.0e9;
    PF = PF / 1.0e9;
    gd = gd / 1.0e-9;
    GD = GD / 1.0e-9;

    # Figure Configure
    # ------ ---------

    z = 1; # number of sigmas of filled area.

    fig = plt.figure(time, figsize=(13, 9));
    plt.xlim(pf[0] - 1, pf[CC] + 1);
    plt.ylim(0.5, 2.2);
    plt.xlabel(r'Probe Frequency (GHz)', fontsize=18);
    plt.ylabel(r'Group Delay (ns)', fontsize=18);
    
    # Plot data Points
    # ---- ---- ------

    init_Ka = np.where(pf > 25.5)[0].min();
    plt.plot(pf[:init_Ka], gd[:init_Ka], 'co', markersize=3.0);
    plt.plot(pf[init_Ka:CC], gd[init_Ka:CC], 'bo', markersize=3.0);
    if (CC > 0): plt.plot(pf[CC:], gd[CC:], 'm*', markersize=7.0);
    
    # Optimize parabolic + gaussian profile fit plot
    # -------- --------- - -------- ------- --- ----

    line1 = '$n_{max} =\\ %.2f \\cdot 10^{19} [\\mathrm{m}^{-3}]$\n';
    line2 = '$\\displaystyle \\left(\\frac{\\sigma}{a}\\right)\\ =\\ %.2f$\n';
    line3 = '$\\qquad \\beta  \\ =\\ %.2f$\n';
    line4 = '$\\qquad \\alpha \\ =\\ %.2f$';
    line_name = line1 + line2 + line3 + line4;
    line_name = line_name % (n_max/1E19, p[3]/a, p[2]**2, p[1]);

    plt.fill_between(PF, (GD + z*Err), (GD - z*Err), color='k', alpha=0.15);
    plt.plot(PF, GD, 'k-', linewidth=1.5, label=line_name);
    plt.legend(loc='upper left', fontsize=18);
    plt.grid();
    
    # Save Figure
    # ---- ------

    fig_name = path + 'time_%.2f.png' % time;

    try: fig.savefig(fig_name, dpi=150, bbox_inches='tight');
    except IOError: fig.savefig(fig_name, dpi=150, bbox_inches='tight');

    if(show_it): plt.show(fig);
    else:        plt.close(fig);

    return {'params': p, 'mcov': mcov, 'res': residuals}
