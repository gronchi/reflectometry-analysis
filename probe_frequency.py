import numpy as np

def PreviousPow2(num):
    i = 1;
    while(2**i < num): i += 1;
    return i - 1;

def CutOverlap(f, cut_off):
    if (f[-1] > cut_off):
        points = np.where(f < cut_off)[0];
        return (points.size, f[points]);
    else: return (0, f);

class PF:
    f0 = 8.4E9;
    ff = 13.5E9;
    
    def __init__(self, n_sweep, fs, ts):
        self.ts   = ts;
        self.nfft = 2**PreviousPow2(float(n_sweep) / 5);
        self.data = self.__MakeProbeFrequency(n_sweep, fs);

    def __MakeProbeFrequency(self, n_sweep, fs):
        t0 = (self.nfft / 2) / fs;
        tf = (n_sweep - (self.nfft / 2)) / fs;
        length = (n_sweep - self.nfft) / 2 + 1;
        Time = np.linspace(t0, tf, length);
        pf_K = self.ProbeFrequencyK(Time);
        pf_Ka = self.ProbeFrequencyKa(Time);
        # remove possible overlap between bands (preference to Ka).
        self.noverlap_K, pf_K = CutOverlap(pf_K, pf_Ka[0]);
        self.init_Ka = pf_K.size;
        return np.concatenate([pf_K, pf_Ka]);

    def ProbeFrequencyK(self, t):
        f_sweep_rate = 2 * (self.ff - self.f0) / self.ts;
        return t * f_sweep_rate + 2 * self.f0;

    def ProbeFrequencyKa(self, t):
        f_sweep_rate = 3 * (self.ff - self.f0) / self.ts;
        return t * f_sweep_rate + 3 * self.f0
