"""

PYTHON CLASS & MODULE
------ ----- - ------

Dependencies : os, numpy, matplotlib, shot_data, optimization_GD
------------

Description :
-----------

    An object of this class can take beat frequency or group delay at any
    time (in ms) for a given shot passed on constructor. This is a minimal
    version in constant development.

How to use :
--- -- ---
    
    You must have connect with server to take data if tou not have the shot.
    >>> import signal_analyse as sa
    >>> M = sa.SF_analysis(shot_number)

    Other way is to pass shot data folder on this computer with saved format
    of shot_data.py class. Look there for more information.
    >>> import signal_analyse as sa
    >>> M = sa.SF_analysis(path_to_shot)

    At any time you can save Signal data to load after as metioned above.
    This way provides fast loading compared to MDSplus like.
    If you dont give any path on computer will save on current folder.
    >>> M.SaveShot(path(optionally))

    OBS: Some cuts of data are done before extract group delay because of
         instability of boundary and loss of amplitude signal.

DEVELOPED BY: ALEX ANDRIATI - USP

"""

import shot_data as sd
import os; 
import numpy as np;
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import optimization_gd as opt


def time2index(rate, t):
    """ For t in milliseconds
        For rate in Hz. """
    return int(rate * t * 1E-3);

def PreviousPow2(num):
    """ Return an integer 'm' that is
    the floor of number 'num' with 2^m """
    i = 1; 
    while(2**i < num): i += 1;
    return (i - 1);

class SF_analysis:
    """ Data Structure to analysis a sweep frequency signal.
        There are here main methods of interesr to analyse
        Group Delay for any given time in milliseconds. Type
        help of EvalGD and TakeBF methods (main mathods). """
    c = 299792458.0;    # Light speed
    a = 0.18;           # Plasma Radius
    R_wall = 0.22;      # From center of plasma to camara wall.

    def __init__(self, shot, path = '', Tvac = 2):
        """ Constructor. Optionally you can define the path
            to record resulds (second argument) or Time of
            vacuum signal reference (third argument). 
            First Argument takes a integer to connect to
            server by MDSplus and read the data channels
            or a string being a path to saved data from
            SaveShot method like.                   """
        # Critical process to access serve.
        self.SD = sd.shot_data(shot);
        # Validate the results path.
        if path != '': self.DefinePath(path);
        else:          self.path = path;
        # Imaging properties for spectrograms.
        self.nfft = 2 ** PreviousPow2(float(self.SD.Nsweep) / 5);
        self.fft_step = 2;              # 'walk' 2 points to next window
        self.PF = self.__ProbeFreq();   # Units in always GHz.
        self.vacBF = self.__TakeVacuumBF(Tvac);

    def DefinePath(self, folder):
        """ Define a path on current computer to save data analysis. """
        if not os.path.exists(folder): raise IOError('Folder Doesnt exist!');
        if folder[-1] != '/': folder = folder + '/';
        # Make new folder to keep results.
        FullPath = folder + '#' + str(self.SD.shot_number) + '_results/';
        warn_msg = '\nCareful, alredy has an initiated analysis.';
        try: os.mkdir(FullPath);
        except OSError: print warn_msg;
        self.path = FullPath;

    def NearSweep(self, t):
        """ For t im milliseconds
            Return index of nearst 
            sweep start. """
        i = time2index(self.SD.rate, t);
        while (self.SD.trig[i] > 0): i += 1;
        while (self.SD.trig[i] < 0): i -= 1;
        return (i + 1);

    def __ProbeFreq(self):
        # First time centered on first FFT window
        t0 = (self.nfft / 2) / self.SD.rate;
        # Final time centered on last FFT window
        tf = (self.SD.Nsweep - (self.nfft / 2)) / self.SD.rate;
        # Number of points depending on step number points (fft_step)
        length = (self.SD.Nsweep - self.nfft) / self.fft_step + 1;
        # Set the vector of time to one sweep in us
        Tsweep = 1E6 * np.linspace(t0, tf, length);
        sweepRate = (self.SD.ff - self.SD.fi) / self.SD.st;
        # K band ampl. by factor of 2
        PF_K  = 2 * (Tsweep * sweepRate + self.SD.fi);
        # Mark where Ka band initiate
        self.IKA = PF_K.size;
        # KA band ampl. by factor of 3
        PF_KA = 3 * (Tsweep * sweepRate + self.SD.fi);
        return np.concatenate([PF_K, PF_KA]);

    def __TakeVacuumBF(self, t): return self.TakeBF(t);

    def TakeBF(self, time, show=False, saveIm=False):
        """ Call signature example: BF = M.TakeBF(70)
            -----------------------------------------

            Return the beat frequencia of signal to the given instant (time).
            Also is considered a mean of power spectrum of ten sweeps, to
            avoid some resolution problems. """
        S_mean_K  = 0.;
        S_mean_KA = 0.;
        i1 = self.NearSweep(time);
        i2 = i1 + int(self.SD.st * 1E-6 * self.SD.rate);
        # time elapsed to next sweep in milliseconds
        T_elapse = (self.SD.st + self.SD.si) * 1E-3;
        # mean of 10 spectrograms.
        for j in range(10):
            # K band
            S, f, t = mlab.specgram(self.SD.K[i1:i2], NFFT=self.nfft,
                    Fs=self.SD.rate, noverlap=self.nfft-self.fft_step, 
                    pad_to=2**12);
            S_mean_K += S / 10;
            # Ka band
            S, f, t = mlab.specgram(self.SD.KA[i1:i2], NFFT=self.nfft,
                    Fs=self.SD.rate, noverlap=self.nfft-self.fft_step, 
                    pad_to=2**12);
            S_mean_KA += S / 10;
            # update indexes to next sweep.
            i1 = self.NearSweep(time + T_elapse);
            i2 = i1 + int(self.SD.st * 1E-6 * self.SD.rate);
            T_elapse += (self.SD.st + self.SD.si) * 1E-3;
        # avoid useles frequency depends on the case
        if time < 10: 
            limSup = np.where(f > 1.65E7)[0].min();
            limInf = np.where(f < 0.65E7)[0].max();
        else:
            limSup = np.where(f > 1.55E7)[0].min();
            limInf = np.where(f < 0.35E7)[0].max();
        # Take max line in spectrum.
        freqIndexK  = S_mean_K[limInf:limSup].argmax(axis=0) + limInf;
        freqIndexKA = S_mean_KA[limInf:limSup].argmax(axis=0) + limInf;
        # return to default image inferior limite
        limSup = np.where(f > 1.55E7)[0].min();
        limInf = np.where(f < 0.35E7)[0].max();
        # Take frequency of max line.
        BF_K  = f[freqIndexK];
        BF_KA = f[freqIndexKA];
        # remove overlap if it has.
        BF = np.concatenate([BF_K, BF_KA]);
        if (saveIm or show):
            S_K  = S_mean_K[limInf:limSup];
            S_KA = S_mean_KA[limInf:limSup];
            f1 = f[limInf]; # Inferior figure limit
            f2 = f[limSup]; # superior figure limit
            self.__DrawImage(S_K, S_KA, time, f1, f2, BF, show, saveIm);
        return BF;

    def EvalGD(self, t, fit_it=False, show=False, saveIm=False):
        """ Call signature example: PF, GD = M.EvalGD(70, booleans...)
            ----------------------------------------------------------
            
            Use TakeBF method and compare of the given instant with vacuum
            to resolve initialization problem and take group delay.
            If you go to fit the curve, can take some time and signature
            calling changes with the results parameters include. Use: 
            PF, GD, Result_fit = M>EvalGD(70, True, True, True) """
        # Rate of probe frequency fo each band
        sweepRate = (self.SD.ff - self.SD.fi) / self.SD.st;
        # Restrict K band analyses. Avoid plasma boundary.
        pf, vacBF = self.__GoodBounds(self.vacBF, 19, 25.3);
        pf, BF = self.__GoodBounds(self.TakeBF(t), 19, 25.3);
        IKA = np.where(pf > 25.3)[0].min(); # New KA start Index
        GD_K  = 1E-6 * (BF[:IKA] - vacBF[:IKA]) / (2 * sweepRate) + \
                2 * (self.a + self.R_wall) * 1E9 / self.c;
        GD_KA = 1E-6 * (BF[IKA:] - vacBF[IKA:]) / (3 * sweepRate) + \
                2 * (self.a + self.R_wall) * 1E9 / self.c;
        GD = np.concatenate([GD_K, GD_KA]);
        if(fit_it): # Critical step. Can take more than a minute to fit.
            result = opt.fit_GD(pf*1E9, GD*1E-9, t, self.path, show, saveIm);
            return (pf, GD, result);
        return (pf, GD);

    def __GoodBounds(self, bf, f1, f2):
        """ Method to cut problematic data """
        i1_K = np.where(self.PF <= f1)[0].max();
        i2_K = np.where(self.PF <= f2)[0].max();
        pf_K = self.PF[i1_K:i2_K];
        pf_Ka = self.PF[self.IKA:];
        bf_K = bf[i1_K:i2_K];
        bf_Ka = bf[self.IKA:];
        good_pf = np.concatenate([pf_K, pf_Ka]);
        good_bf = np.concatenate([bf_K, bf_Ka]);
        return (good_pf, good_bf);

    def __DrawImage(self, S_K, S_Ka, time, bf_inf, bf_sup, bf, show, saveIm):
        # Scale to simple units.
        pf      = self.PF;
        BF_sup  = bf_sup / 1E6;
        BF_inf  = bf_inf / 1E6;
        BF      = bf / 1E6;
        # image displayed one side for each band.
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(18, 8));
        # Draw Spectrum of K band.
        ax1.imshow(np.log(S_K), origin='lower', aspect='auto', \
                extent=[pf[0], pf[self.IKA-1], BF_inf, BF_sup]);
        # Draw Spectrum of Ka band.
        ax2.imshow(np.log(S_Ka), origin='lower', aspect='auto', \
                extent=[pf[self.IKA], pf[-1], BF_inf, BF_sup]);
        # Draw line of maximum for each spectrum.
        ax1.plot(pf[:self.IKA], BF[:self.IKA], 'k-', linewidth=2);
        ax2.plot(pf[self.IKA:], BF[self.IKA:], 'k-', linewidth=2);
        # Figure limits.
        f.subplots_adjust(wspace=0.02);
        ax1.set_xlim(pf[0], pf[self.IKA-1]);
        ax1.set_ylim(BF_inf, BF_sup);
        ax2.set_xlim(pf[self.IKA], pf[-1]);
        ax2.set_ylim(BF_inf, BF_sup);
        # Figure labels (x comom label)
        ax1.set_ylabel('Beat Frequency (MHz)', fontsize=16);
        f.text(0.5, 0.03, 'Probe Frequency (GHz)', ha='center', fontsize=16);
        ax1.set_title('\"K\" band', fontsize=16); 
        ax2.set_title('\"Ka\" Band', fontsize=16);
        # Switch case if vacuum or plasma to choose a name.
        if time > 10 : 
            figName = '#' + str(self.SD.shot_number) + '_%.2f.jpg' % time
        else : 
            figName = '#' + str(self.SD.shot_number) + '_vacuum.jpg'
        if saveIm: 
            f.savefig(self.path + figName, dpi=150, bbox_inches='tight');
        # if required display image.
        if (show): plt.show(f);
        else:      plt.close(f);

    def SaveShot(self, path=''): self.SD.SaveShot(path);
