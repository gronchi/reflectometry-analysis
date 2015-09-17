"""

Python Script
------ ------

Dependencies: numpy, signal_analyse, sys.
-------------

Execution: $ python time_evolution shot_number t0 tf dt
----------

    Evaluate Group Delay by optimized curve fitted method. Get just the
    center density from result parameters and estimate error.

    shot_number --> shot to accsess on MDSplus
    to -----------> Initial instant in miliseconds
    tf -----------> final instant in miliseconds
    dt -----------> time interval between calculations

Developed by: Alex Andriati - USP

"""

import sys;
from numpy import arange, sqrt;
import signal_analyse as sa;

folder = raw_input('\nPath of folder to send file results: ');
file_name = raw_input('\nFile name to record results: ');

t1 = float(sys.argv[2]);
t2 = float(sys.argv[3]);
dt = float(sys.argv[4]);
M  = sa.SF_analysis(int(sys.argv[1]), folder);

f = open(M.path + file_name, 'w');

# Warning when do not 
war_msg = '\nParametros para o instante %.3f nao encontrado.'
war_msg = war_msg + ' Maximo de tentativas excedido!'

for t in arange(t1, t2, dt):
    print '\nInstante = %.3fms' %t;
    try: pf, gd, result = M.EvalGD(t, True, saveIm=True);
    except RuntimeError: print war_msg; continue;
    # Extract useful results.
    params = result['params'];
    mcov   = result['mcov'];
    n0_sig = sqrt(mcov[0][0]);
    A_sig  = sqrt(mcov[2][2]);
    n_max  = params[0] * (1 + params[2]**2);
    # Estimate uncertain by covariance matrix.
    der_n0 = (1 + params[2]**2);
    der_A  = params[0] * (2 * params[2]);
    sigma  = sqrt((der_n0 * n0_sig)**2 + (der_A * A_sig)**2);
    # Record result.
    f.write('%f' %t);
    f.write('\t' + str(n_max));
    f.write('\t' + str(sigma));
    f.write('\n');
f.close()
