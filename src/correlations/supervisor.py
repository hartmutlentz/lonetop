import sys, os
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))

import tools.filesystem as fs
import tools.Gewindehammer as gwh
from MatrixProduct import ProductOfAdjacencyMatrices
from MatrixList import AdjMatrixSequence
from scipy.io import mmwrite

########

At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,3))
print 'Matrixsequenz eingelesen'
print At[0]
B=At.time_reversed(True)
print '\n',B[0]
#Z=ProductOfAdjacencyMatrices(At)
#print 'Produkt-Objekt erzeugt'
#    
#D=Z.floyd_warshall_temporal()
#print D
#mmwrite("Test.mtx",D)



