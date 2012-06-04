import numpy as np
import subprocess, sys, os
import scipy
from scipy import sparse
import gc
import tools.progressmeter as pm


def csr_sigma(x,y):
    """
    sigma_xy for sparse matrix
    values need to be 0 or 1
    """
    m,n = x.shape
    size = float(m*n)
    mean_x = x.nnz/size
    mean_y = y.nnz/size
    sigma_xy = (x.multiply(y).nnz - mean_y*x.nnz - mean_x*y.nnz)/size + mean_x*mean_y
    return sigma_xy

def sparse_corr_coeff(x,y):
    """
    corr coeff of two sparse matrices, nach Bronstein
    values considered boolean
    """
    sigma_xy = csr_sigma(x,y)
    sigma_x = np.sqrt(csr_sigma(x,x))
    sigma_y = np.sqrt(csr_sigma(y,y))
    coeff = sigma_xy / (sigma_x * sigma_y)
    return coeff

def matrix_correlation(x,y):
    """
    returns pearson correlation coefficient for 2 matrices
    """
    if x.nnz == 0 or y.nnz == 0:
        return 0
    else:
        x = sparse_reshape(x)
        y = sparse_reshape(y)
        x = sparse.hstack((x,y))
        c = corrcoefcsr(x)
        return c[0,-1]

def corr_matrix(seq, fname):
    """
    calculates all possible pair-wise correlations for matrices in seq

    the matrices in seq are meant to be sparse
    no return value, but saved as i,j,c on-the-fly
    this makes interruption of code possible
    """

    #find the highest value of i already saved
    try:
        i_done = max(np.loadtxt(fname, usecols=[0], dtype=int, unpack=True))
    except IOError:
        i_done = 0

    #start calculation at this value
    p = pm.ProgressMeter(unit="rows", total=(len(seq)-i_done), txt="Calulating correlation matrixes.")
    for i in range(i_done, len(seq)):
        file = open(fname,"a")
        p.update()
        #calculate only lower triangle and save each value directly
        for j in range(0, i):
            c = sparse_corr_coeff(seq[i],seq[j])
            file.write(str(i)+" "+str(j)+" "+str(c)+"\n")
        file.close()
        gc.collect()
    return True

