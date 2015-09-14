# coding: utf-8

"""

Executable   : No. Just Python Class.
----------

Dependencies : os, numpy, matplotlib, probe_frequency, pre_config
------------

Description :
-----------

    An object of this class can take beat frequency or group delay at any
    time (in ms) of a given shot passed on constructor. This is a minimal
    version in constant development.

How to use :
--- -- ---
    
    For now just '.mat' files are acceptable and for sweep time of
    8 microseconds. Some default importants are the frequency range
    8.4 - 13.5 of reflectometer and some cut-off of problematic
    beat frequency intevals, excluding in K band range f < 19 GHz
    and f > 25.5 GHz. To better this, need be review the line
    maximum method to extract beat frequency of spectrogram.

    Vacuum signal extract of 2 miliseconds by default.

    .TakeBF(time in ms, (optional) save image?, (optional) show ?)
        returns the tuple (probe freq., beat freq) data.
        The optional booleans are to draw image and show

    .EvalGD(time in ms, (optional)fit ?, (optional)show ?)
        returns the tuple (probe freq., group delay) data.
        curve_fit may take until a minute or more to result.

"""

import os, numpy as np, matplotlib.mlab as mlab, pylab
import optimization_gd as opt
from probe_frequency import PF

class signal_analyse:
    # Constants to eval group delay
    c = 299792458.0;
    a = 0.18;
    R_wall = 0.22;

    def __init__(self, shot_data):
        self.data = shot_data;
        self.r_dir = self.MakeDirectory();
        self.pf = PF(self.data.n_sweep, self.data.fs, self.data.ts);
        self.vac_bf = self.__TakeVacuumBF();

    def MakeDirectory(self):
        path = self.data.this_file_path + '#' + self.data.name + '/';
        warn_msg = '\nAnalise de disparo ja iniciada.'  + \
                'Possivel sobreposicao de resultados e figuras.'
        # Make directory to keep image and results.
        try: os.mkdir(path);
        except OSError: print warn_msg;
        return path;

    def __TakeVacuumBF(self): return self.TakeBF(2);

    def TakeBF(self, time, save_im=True, show=False):
        S_mean_K  = 0.
        S_mean_Ka = 0.
        i1, i2 = self.data.NextSweep(time);
        time_elapse = 2 * self.data.ts;
        # mean of 10 spectrograms.
        for j in range(10):
            # K band
            S, f, t = mlab.specgram(self.data.K[i1:i2], NFFT=self.pf.nfft,
                    Fs=self.data.fs, noverlap=self.pf.nfft-2, pad_to=2**12);
            S_mean_K += S / 10;
            # Ka band
            S, f, t = mlab.specgram(self.data.Ka[i1:i2], NFFT=self.pf.nfft,
                    Fs=self.data.fs, noverlap=self.pf.nfft-2, pad_to=2**12);
            S_mean_Ka += S / 10;
            # update indexes to next sweep.
            i1, i2 = self.data.NextSweep(time + 1E3 * time_elapse);
            time_elapse += 2 * self.data.ts;
        # Finish the mean.
        # avoid useles frequency depends on the case
        if time < 10: 
            limSup = np.where(f > 1.65E7)[0].min();
            limInf = np.where(f < 0.65E7)[0].max();
        else:
            limSup = np.where(f > 1.55E7)[0].min();
            limInf = np.where(f < 0.35E7)[0].max();
        # Take max line in spectrum.
        freqIndexK  = S_mean_K[:][limInf:limSup].argmax(axis=0) + limInf;
        freqIndexKa = S_mean_Ka[:][limInf:limSup].argmax(axis=0) + limInf;
        # return to default image inferior limite
        limSup = np.where(f > 1.55E7)[0].min();
        limInf = np.where(f < 0.35E7)[0].max();
        # Take frequency of max line.
        BF_K  = f[freqIndexK];
        BF_Ka = f[freqIndexKa];
        # remove overlap if it has.
        if (self.pf.noverlap_K > 0): BF_K = BF_K[:-self.pf.noverlap_K-1];
        bf = np.concatenate([BF_K, BF_Ka]);
        if (save_im):
            S_K  = S_mean_K[:][limInf:limSup];
            S_Ka = S_mean_Ka[:][limInf:limSup];
            # Concat. Spectrum of two bands.
            self.__DrawImage(S_K, S_Ka, time, f[limInf], f[limSup], bf, show);
        return bf;

    def EvalGD(self, time, fit_it=False, show=False):
        # Rate of probe frequency fo each band
        sweep_rate_K  = 2 * (self.pf.ff - self.pf.f0) / self.data.ts
        sweep_rate_Ka = 3 * (self.pf.ff - self.pf.f0) / self.data.ts
        # Restrict K band analyses.
        pf, vac_bf = self.__GoodBounds(self.vac_bf, 19E9, 25.3E9);
        pf, bf = self.__GoodBounds(self.TakeBF(time), 19E9, 25.3E9);
        i_Ka = np.where(pf > 25.5E9)[0].min();
        GD_K  = (bf[:i_Ka] - vac_bf[:i_Ka]) / sweep_rate_K + 2 * (self.a +
                self.R_wall) / self.c;
        GD_Ka = (bf[i_Ka:] - vac_bf[i_Ka:]) / sweep_rate_Ka + 2 * (self.a +
                self.R_wall) / self.c;
        GD = np.concatenate([GD_K, GD_Ka]);
        if(fit_it): 
            opt_result = opt.fit_GD(pf, GD, self.r_dir, time, show);
            return (pf, GD, opt_result);
        return (pf, GD);

    def __GoodBounds(self, bf, f1, f2):
        i1_K = np.where(self.pf.data < f1)[0].max();
        i2_K = np.where(self.pf.data < f2)[0].max();
        pf_K = self.pf.data[i1_K:i2_K];
        pf_Ka = self.pf.data[self.pf.init_Ka:];
        bf_K = bf[i1_K:i2_K];
        bf_Ka = bf[self.pf.init_Ka:];
        good_pf = np.concatenate([pf_K, pf_Ka]);
        good_bf = np.concatenate([bf_K, bf_Ka]);
        return (good_pf, good_bf);

    def __DrawImage(self, S_K, S_Ka, time, bf_inf, bf_sup, bf, show=False):
        # Scale to simple units.
        pf      = self.pf.data / 1E9;
        bf_sup  = bf_sup / 1E6;
        bf_inf  = bf_inf / 1E6;
        disp_bf = bf / 1E6;
        # image displayed one side for each band.
        f, (ax1, ax2) = pylab.subplots(1, 2, sharey=True, figsize=(18, 8));
        # Draw Spectrum of K band.
        ax1.imshow(np.log(S_K), origin='lower', aspect='auto',
                extent=[pf[0], pf[self.pf.init_Ka], bf_inf, bf_sup]);
        # Draw Spectrum of Ka band.
        ax2.imshow(np.log(S_Ka), origin='lower', aspect='auto',
                extent=[pf[self.pf.init_Ka], pf[-1], bf_inf, bf_sup]);
        # Draw line of maximum for each spectrum.
        ax1.plot(pf[:self.pf.init_Ka], disp_bf[:self.pf.init_Ka], 
                'k-', linewidth=2);
        ax2.plot(pf[self.pf.init_Ka:], disp_bf[self.pf.init_Ka:], 
                'k-', linewidth=2);
        # Figure limits.
        f.subplots_adjust(wspace=0.02);
        ax1.set_xlim(pf[0], pf[self.pf.init_Ka]);
        ax1.set_ylim(bf_inf, bf_sup);
        ax2.set_xlim(pf[self.pf.init_Ka], pf[-1]);
        ax2.set_ylim(bf_inf, bf_sup);
        # Figure labels (x comom label)
        ax1.set_ylabel('Beat Frequency (MHz)', fontsize=16);
        f.text(0.5, 0.03, 'Probe Frequency (GHz)', ha='center', fontsize=16);
        ax1.set_title('\"K\" band', fontsize=16); 
        ax2.set_title('\"Ka\" Band', fontsize=16);
        # Switch case if vacuum or plasma to choose a name.
        if time > 10 : figName = '#' + self.data.name + '_%.2f.jpg' % time
        else : figName = '#' + self.data.name + '_vacuum_%.2f.jpg' % time
        f.savefig(self.r_dir + figName, dpi=150, bbox_inches='tight');
        # if required display image.
        if (show): pylab.show(f);
        else:      pylab.close(f);
