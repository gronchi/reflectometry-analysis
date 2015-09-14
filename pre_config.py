# coding: utf-8

"""

Executable   : python pre_config.py
----------

Dependencies : os
------------

    Change pre_config.dat file to configure shots location in current computer
    to run density profile analyses based on Group Delay of Reflectometer.
    If file doesn't exist creat on current directory or from chosen.

    NOTE : The pre_config.dat needs to be located on same directory of
    system archives. Please be sure of that.

"""

import os

ans = raw_input('\nChange folder ? [y|n]: ');
if (ans == 'y' or ans == 'Y'):
    config_file_path = raw_input('\nEnter a path form root folder: ');
    os.chdir(config_file_path);
    f = open('pre_config.dat', 'w');
else:
    f = open('pre_config.dat', 'w');

while True:
    shots_path = raw_input('\nEnter one of shots location on this computer: ');
    if (os.path.exists(shots_path)):
        f.write('shot_path: ' + shots_path + ' ');
        ans = raw_input('\nEnter other location ? [y|n]: ')
        if (ans == 'n' or ans == 'N'): break;
    else:
        print '\ndirectory does not exist.'
        ans = raw_input('Try again ? [y|n]: ')
        if (ans == 'n' or ans == 'N'): break;

f.write('\nextension: .mat');
f.close();
