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
    sys.path.append('/users/tsd/lentz/Documents/GitHub_locals/my_py')
    sys.path.append('/users/tsd/lentz/Documents/GitHub_locals/lonetop/src/correlations')
    sys.path.append('/users/tsd/lentz/Documents/GitHub_locals/lonetop/src/tools')
else:
    raise NameError, 'unknown operating system.'


import numpy as np, scipy as sc, Gewindehammer as gwh
from MatrixList import AdjMatrixSequence
from TemporalEdgeList import TemporalEdgeList

#
# ----------------- ----------------- ----------------- -----------------
#

if __name__=="__main__":
    #the_file='/Users/lentz/Desktop/Dissertation-D_sw_Edgelists/Temporal/D_sw_uvd_01JAN2008_31DEC2009.txt'
    the_file='T_edgelist.txt'
    At=AdjMatrixSequence(the_file,directed=True,write_label_file=False)
    #print '-> read file: done'
    #At.time_shuffled()
    #At.time_reversed()
    #print '-> shuffling done. Writing...'
    #At.write("PVM_GST.txt")
    #print '-> file written.'
    #print len(At),At.number_of_nodes
    #c=At.unfold_accessibility()
    #h=gwh.cdf2histogram(c)
    #gwh.dict2file(c,"cumu.txt")
    #gwh.dict2file(h,"histo.txt")

    E=TemporalEdgeList(the_file,directed=True)
    #print '-> read file: done'
    E.RE()
    #print '-> shuffling done. Writing...'
    #E.write("Test_RE.txt")
    #print '-> file written.'



