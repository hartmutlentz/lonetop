#! /usr/local/bin/python
#
import sys,os
sys.path.append('/Users/lentz/Documents/Physik/my_py') # Mac
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))
import tools.filesystem as fs
import tools.Gewindehammer as gwh

import numpy as np, networkx as nx
from pprint import pprint
from MatrixList import AdjMatrixSequence
import scipy.sparse as sp
import math,csv,gc,itertools
from scipy.io import mmread,mmwrite,loadmat,savemat



class ProductOfAdjacencyMatrices(list):
    """ Like AdjMatrixSequence, but with focus on products and spectra.
        Also contains creation of random matrices.
        
        Input: AdjMatrixSequence or networkx Graph Generator with parameters in args.
        
        ----------------
        Examples:
        Z=ProductOfAdjacencyMatrices(nx.fast_gnp_random_graph,n=100,p=0.01,directed=True)
        
        At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,3))
        Z=ProductOfAdjacencyMatrices(At)
        
    """
    def __init__(self,matrices,diagonal=1,timespan=100,**args):
        list.__init__(self)
        
        if isinstance(matrices,AdjMatrixSequence): self.extend(matrices)
        else:
            self.matrix_generation(matrices,size=timespan,**args)
        
        self.diagonal=diagonal
        self.number_of_nodes=self[0].shape[0]
        
        if diagonal: self.__add_diagonal()
    
    def __add_diagonal(self):
        # Adds diagonal Matrix to each matrix of self
        D=sp.identity(self.number_of_nodes,dtype=np.int32)
        D*=self.diagonal
        for i in range(len(self)):
            self[i]=self[i]+D
            
    def temporal_path_matrix(self,C):
        """ need cumulated Graph as Matrix C """
        D=sp.identity(self.number_of_nodes,dtype=np.int32)
        return self.full_product_matrix() - D - C
    
    def cumulated(self):
        """ Better: get cumulated from Matrixlist! """
        print "Warning: please compute cumulated matrix from Matrixlist Object!"
        P=self[0].copy()
        for i in range(1,len(self)):
            P=P+self[i]
        D=sp.identity(self.number_of_nodes,dtype=np.int32)
        D=D*self.diagonal *len(self)   
        return P-D
        
    def matrix_transitivity(self,A):
        # Transitivity of a Matrix A
        T=A.multiply(A**2)
        return float(T.nnz)/float((A**2).nnz)
    
    def power_method(self):
        """ Multiplies a random vector with all matrices in self.
        
        """
        x=self.random_vector()
        x/=np.linalg.norm(x)
        results={}
        for i in range(len(self)):
            print i
            #norm1=np.linalg.norm(x)
            x=self[i]*x
            norm2=np.linalg.norm(x)
            
            print norm2
            results[i]=norm2,gwh.the_fle(x)
        return results
            
    
    def full_product_matrix(self):
        P=self[0].copy()
        #for A in self[1:]:
        #    P*=A
        for i in range(1,len(self)):
            print 'full product iteration',i,'non-zeros: ',P.nnz
            P=P*self[i]
        return P
        
    def random_vector(self):
        return np.random.rand(self.number_of_nodes)
        
    def write_large_dense_matrix(self,A,nameoffile='matrix.txt'):
        # Write sparse matrix to textfile
        writer=csv.writer(open(nameoffile,"wb"))
        for i in xrange(n):
            pass
        
        #writer=csv.writer(open(nameoffile,"wb"))
        #indices=zip(A.nonzero()[0],A.nonzero()[1])
        #for i,j in indices:
        #    writer.writerow([i,j,A[i,j]])        
        
        #g=file(nameoffile,'w+')
        #for i,j in A.nonzero():
        #   g.writelines((str(i),'\t',str(j),'\t',str(A[i,j]),'\n'))
        #g.close
        
    def matrix_generation(self,generator,size,**prms):
        for i in range(size):
            self.append(nx.to_scipy_sparse_matrix(generator(**prms)))
        
def write_large_dense_matrix(A,fname='matrix.txt'):
    """ writes dense matrix to text file
    
    """
    n=A.shape[0]
    writer=csv.writer(open(fname,"wb"))
    for i in xrange(n):
        for j in xrange(n):
            writer.writerow([i,j,A[i,j]])
        print i

if __name__=="__main__":
    #Z=ProductOfAdjacencyMatrices(nx.fast_gnp_random_graph,n=100,p=0.01,directed=True)
    At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,2),matr_type='dok')
    #At = AdjMatrixSequence(fs.dataPath("D_sw_uvd_01JAN2009_31MAR2010.txt"),matr_type='dok')
    #C=At.cumulated()
    
    print 'Matrixsequenz eingelesen'
    Z=ProductOfAdjacencyMatrices(At)
    print 'Produkt-Objekt erzeugt'
    gc.collect()
        
    P=Z.full_product_matrix()
    gc.collect()
    
    print 'Schreibe'
    print type(P)
    print P.dtype
    
    writer=csv.writer(open('matrix.txt',"wb"))
    for i,j in itertools.izip(P.nonzero()[0],P.nonzero()[1]):
        writer.writerow([i,j,P[i,j]])




