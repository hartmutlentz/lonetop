#!/usr/bin/env python
#
#
import sys,platform
# Adjust these lines to match the locations of your repositories
if platform.system()=='Darwin': # Mac
    sys.path.append('/Users/lentz/Documents/GitHub_locals/my_py') # optional tools
    sys.path.append('/Users/lentz/Documents/GitHub_locals/lonetop/src/correlations')
    sys.path.append('/Users/lentz/Documents/GitHub_locals/lonetop/src/tools')
elif platform.system()=='Linux': # Linux
    sys.path.append('/users/tsd/lentz/Documents/GitHub_locals/my_py') #optional tools
    sys.path.append('/users/tsd/lentz/Documents/GitHub_locals/lonetop/src/correlations')
    sys.path.append('/users/tsd/lentz/Documents/GitHub_locals/lonetop/src/tools')
else:
    raise NameError, 'unknown operating system.'


import numpy as np, scipy as sc
import Gewindehammer as gwh # optional: for dict2file, cdf2histogram, etc.
from MatrixList import AdjMatrixSequence
from TemporalEdgeList import TemporalEdgeList

#
# ----------------- ----------------- ----------------- -----------------
#

if __name__=="__main__":
    
    the_file='data/T_edgelist.txt'
    At=AdjMatrixSequence(the_file,directed=True,write_label_file=False)
    c=At.unfold_accessibility()
    h=gwh.cdf2histogram(c)
    gwh.dict2file(c,"cumu.txt")
    gwh.dict2file(h,"histo.txt")

    
    # Randomized versions of the temporal network
    """E=TemporalEdgeList(the_file,directed=True)
    print '-> read file: done'
    E.RE()
    print '-> shuffling done. Writing...'
    E.write("Test_RE.txt")
    print '-> file written.' """



