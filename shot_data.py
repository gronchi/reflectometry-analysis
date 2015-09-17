"""

PYTHON CLASS & MODULE
------ ----- - ------

Dependencies : MDSplus -> See www.mdsplus.org
------------

    Implements a data structure/class that keep all the data of reflectometer
    from a given shot/discharge. Need connection with tcabr server or a saved
    file if format like done in SaveShot method.

HOW TO USE :
------------

    >>> from shot_data import shot_data;
    >>> all_accessible_data(s) = shot_data(shot_number);
    >>> s.K, s.KA, s.trig, s.T, ...
    >>> s.SaveShot('some path here')

Developed by: Alex Andriati - USP.

"""

import MDSplus as mds;

class shot_data:

    def __init__(self, shot):
        self.shot_number = shot;
        self.__fromMDSplus();

    def __fromMDSplus(self):
        conn = mds.Connection('tcabrcl.if.usp.br:8000');
        conn.openTree('tcabr_ref', self.shot_number);
        print '\nDonwloading data from the server. Please wait.'
        # Read in MDSplus format first to take time
        # critical proccess, can take minutes depending
        # on connection speed with server.
        self.K  = conn.get('\\KBAND.SIGNAL');
        self.KA = conn.get('\\KABAND.SIGNAL').data();
        self.trig = conn.get('\\TRIGGER.SIGNAL').data();
        # Organize in numpy arrays.
        self.T = self.K.dim_of().data();
        self.K = self.K.data();
        # simple scalars or strings. Common parameters.
        self.rate = conn.get('\\REFPARAMETER.RATE').data();
        self.mode = str(conn.get('\\REFPARAMETER.REFMODE'));
        self.sample = conn.get('\\REFPARAMETER.SAMPLES').data();
        self.angle  = conn.get('\\REFPARAMETER.ANGLE').data();
        # Take correct vector of time (miliseconds)
        self.T = 1E3 * self.T / self.rate;
        # Identify mode and take others parameters.
        if self.mode == 'ff':
            print '\nFixed Frequency mode.'
            self.freq = conn.get('\\FIXEDFREQ.FREQUENCY').data();
        elif self.mode == 'hf':
            print '\nHopping Freq. enable'
            self.freqs = conn.get('\\HOPPINGFREQ.FREQ_TABLE').data();   # GHz
            self.rt = conn.get('\\HOPPINGFREQ.RESTART_TIME').data();    # us
            self.st = conn.get('\\HOPPINGFREQ.TIME_STEP').data();       # us
            self.Nsweep = int(1E-6 * self.rate * st);
        elif self.mode == 'sf':
            print '\nSweep Frequency mode.'
            self.fi = conn.get('\\SWEEPFREQ.FREQ_START').data();    # GHz
            self.ff = conn.get('\\SWEEPFREQ.FREQ_END').data();      # GHz
            self.st = conn.get('\\SWEEPFREQ.SWEEP_TIME').data();    # us
            self.si = conn.get('\\SWEEPFREQ.INTERV_SWEEP').data();  # us
            self.Nsweep = int(1E-6 * self.st * self.rate);
        else: print '\nUnknow mode of operation.'
        conn.closeAllTrees();

    def SaveShot(self, path=''):
        """ Record data of a given shot. Make a file for each channel
            in format .npy, that can easily and rapidly read after. 
            If any path is given save in current directory. """
        if self.mode != 'sf': raise TypeError('Cant save this mode yet');
        # Modules needed to record data.
        import os; from numpy import save;
        if path != '' and not os.path.exists(path): 
            raise IOError('Path dont exist');
        # Correct th path if doesnt finish with '/'
        if path != '' and path[-1] != '/': path = path + '/';
        # Update path to create folder.
        path = path + '#' + str(self.shot_number) + '/';
        if not os.path.exists(path): os.mkdir(path);
        # Start save values
        print '\nSaving data at ' + path;
        save(path + 'K.npy', self.K);
        save(path + 'KA.npy', self.KA);
        save(path + 'trigger.npy', self.trig);
        f = open(path + 'configuration.dat', 'w');
        # write common fixed parameters.
        f.write('rate: %.2g' % self.rate);
        f.write('\nsample: %d' % self.sample);
        f.write('\nangle: %.1f' % self.angle);
        f.write('\nmode: ' + self.mode);
        # Data from sweep frequency mode.
        f.write('\ninit_Freq: %.4f' % self.fi);     # GHz
        f.write('\nfinal_Freq: %4f' % self.ff);     # GHz
        f.write('\nsweep_time: %.4f' % self.st);    # us
        f.write('\nsweep_inter: %.4f' % self.si);   # us
        f.close();
