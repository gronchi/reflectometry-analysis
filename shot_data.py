# coding: utf-8

"""

Executable   : No. Just Python Class.
----------

Dependencies : os, scipy.io, pre_config.dat
------------

    Implements a data structure/class that keep all the data of reflectometer
    from a given shot/discharge in a given interval. Pay attention and assert
    that pre_config.dat is present on the same directory.

    NOTE: If doesn't contains pre_config.dat run pre_config.py file.

"""

import os, scipy.io as sio

class shot_data:
    useless = ['__version__', '__header__', '__globals__'];

    def __init__(self, shot_name):
        self.name = shot_name;
        self.__TakePaths();
        self.__GetNumber();
        if (not hasattr(self, 'gas')): self.gas = 'H'
        self.this_file_path = self.__FindInPaths();
        self.__TakeData(self.this_file_path + self.name + self.ext);

    def __TakePaths(self):
        # From pre_config.dat file take the paths and extension
        # of files in this computer.
        while True:
            try:
                f = open('pre_config.dat', 'r');
                break;
            except IOError:
                print '\nConfigure file is not in current directory.'
                print 'Please make configures now.'
                execfile('pre_config.py');
        self.paths = f.readline().split(' ')[1:-1];
        self.ext   = f.readline().split(' ')[1];

    def __GetNumber(self):
        # Get shot number and possible gas type from name.
        try: self.number = int(self.name);
        except ValueError:
            self.number = int(self.name.split('_')[0]);
            self.gas = self.name.split('_')[1];
            
    def __FindInPaths(self):
        # from a given shot find in suggested paths from pre_config.
        for path in self.paths:
            if os.path.exists(path + self.name + self.ext): return path;
        return path;

    def __TakeData(self, mat_file):
        mat_data = sio.loadmat(mat_file, matlab_compatible=True);
        # Clean up useless fiels returned from scipy.io.loadmat
        for key in mat_data.keys():
            if key in shot_data.useless: continue;
            real_data = mat_data[key];
        # separe band chanels in each field.
        self.time = real_data['t'][0][0].T[0];
        self.trig = real_data['trigger'][0][0].T[0];
        self.Ka   = real_data['Ka'][0][0].T[0];
        self.K    = real_data['K'][0][0].T[0];
        # time data in miliseconds --> Compute Frequency sample
        # Commun sweep time interval.
        self.ts = 8E-6;
        self.fs = round(1.0E3 / (self.time[1] - self.time[0]));
        self.n_sweep = int(round(self.fs * self.ts));

    def time2index(self, t): return int(self.fs * t * 1E-3);

    def NextSweep(self, t):
        # Return index limits of next sweep.
        start_index = self.time2index(t);
        while (self.trig[start_index] > 0): start_index += 1;
        while (self.trig[start_index] < 0): start_index -= 1;
        start_index += 1;
        return (start_index, start_index + self.n_sweep);
