import sys, os
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))
import tools.filesystem as fs

import numpy as np
from scipy import sparse
from scipy import fftpack

def mean_corr(fname):
    i,j,v = np.loadtxt(fname,unpack=True)
    m = sparse.csr_matrix((v, (i, j)))
    m = m.todense()
    m = np.nan_to_num(m)
    means = [np.median(np.diag(m, -k)) for k in range(1, m.shape[0])]
    return means

def fft(seq):
    y=fftpack.fft(seq)
    n=len(y)
    power = abs(y[1:(n/2)])**2
    nyquist=1./2
    freq=np.array(range(n/2))/(n/2.0)*nyquist
    return freq[1:-1], power

def smooth(seq):
    """smoothes with fixed 7-day window"""
    seq = [np.mean(seq[i-3:i+4]) for i in range(3, len(seq)-4)]
    return seq

def amplitude(seq):
    """reduce to amplitudes"""
    l = len(seq)
    mx = max(seq[l/3:2*l/3]) #find one of the majors in central third
    idx = seq.index(mx) #index of this item
    idx = idx % 7
    seq = [(i, seq[i]) for i in range(idx, len(seq)-7, 7)]
    return seq


if __name__ == "__main__":
     means = mean_corr(fs.resultsPath("sf_matrix.txt"))

     smth = smooth(means)
     np.savetxt(fs.resultsPath("sf_means.txt"), zip(means, smth))

     ampl = amplitude(means)
     np.savetxt(fs.resultsPath("sf_ampl.txt"), ampl)

     x,y = fft(means)
     np.savetxt(fs.resultsPath("sf_fft.txt"), zip(x,y))

     x,y = fft(smth)
     np.savetxt(fs.resultsPath("sf_smth_fft.txt"), zip(x,y))

     ampl = [s for i,s in ampl]
     x,y = fft(ampl)
     np.savetxt(fs.resultsPath("sf_ampl_fft.txt"), zip(x,y))

