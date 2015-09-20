REFLECTOMETRY ANALYSIS
------------- --------
------------- --------


An project to extract information from (best) Group Delay analytical curve.
---------------------------------------------------------------------------

Before Anything:

    - Make sure that you have MDSplus installed. If not take it at http://www.mdsplus.org
    - Make sure that you have MDSplus python module. For example try import MDSplus in python shell.
      If not, you can search in http://www.mdsplus.org/dist/ the .egg module file and use easy_install
      to implement it in your python enviroment.
    - Finally, these classes and scripts use several standart python libraries for science like scipy,
      numpy, matplotlib, etc. for best performance. If you dont have them, install before proceed.

How To use:

    - Main module is signal_analyse.py which has the SF_analysis class. I strongly recommend that the
      user type help in python enviroment after import at least once time. To beginners see below:
      >>> import signal_analyse
      >>> help(signal_analyse).

OBS. About code performance:

    - For further improvements on time demanded for run Curve Fit over Group Delay Data, developers
      can take a tour on Cython package to introduce some static data types based on C language.
      
    - A deep investigation can be done in scipy.integrete methods, and in numpy array access order.
    
OBS. About OS:

    - Some manipulations are made in order to write or access files in a folder. Maybe such accesses
      result in bugs in other Operational Systems different from Unix or even others Linux distors.
      All was developed to work in Debian distros.

DEVELOPED BY: Alex Andriati - USP
SUPPORTED BY: Fapesp - Fundação de Amparo a Pesquisa do estado de São Paulo
