import sys, os
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))

import tools.filesystem as fs
import matrixlist, corr
import gc

if __name__ == "__main__":

    cm = matrixlist.AdjMatrixSequence(fs.dataPath("D_sf_uvwd_cmpl.txt"))
    gc.collect()
    corr.corr_matrix(cm, fs.resultsPath("sf_slow_matrix.txt"))
    gc.collect()

    cm = matrixlist.AdjMatrixSequence(fs.dataPath("D_sw_uvwd_cmpl.txt"))
    gc.collect()
    corr.corr_matrix(cm, fs.resultsPath("sw_matrix.txt"))
    gc.collect()

    cm = matrixlist.AdjMatrixSequence(fs.dataPath("D_ri_uvwd_cmpl.txt"))
    gc.collect()
    corr.corr_matrix(cm, fs.resultsPath("ri_matrix.txt"))