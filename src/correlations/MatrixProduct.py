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
        
    def matrix_transitivity(self,M):
        # Transitivity of a Matrix M
        T=M.multiply(M**2)
        return float(T.nnz)/float((M**2).nnz)

    def bool_int_matrix(self,M):
        """ Returns matrix with only np.int64: ones. """
        M=M.astype('bool')
        M=M.astype('i')
    
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
                
    def full_product_matrix(self,return_path_matrix=True,return_transitivity=True):
        P=self[0].copy()
        B=sp.csr_matrix((self.number_of_nodes,self.number_of_nodes))
        cumu=[0]
        trans=[0.0]

        if return_transitivity: 
            for i in range(1,len(self)):
                print 'full product iteration',i,'non-zeros: ',P.nnz
                cumu.append(P.nnz)
                try:
                    trans.append(self.matrix_transitivity(P))
                except MemoryError:
                    print 'Excepted Memory Error.'
                    break
                    
                P=P*self[i]
            
            if return_path_matrix:
                P = P.astype('bool')
                P = P.astype('int')    
                return P,cumu,trans
            else:
                return cumu,trans
        else:
            for i in range(1,len(self)):
                print 'full product iteration',i,'non-zeros: ',P.nnz
                cumu.append(P.nnz)
                #trans.append(self.matrix_transitivity(P))
                P=P*self[i]
                B=(P.astype('bool')-B.astype('bool'))*i + B
            
            if return_path_matrix:
                P = P.astype('bool')
                P = P.astype('int')    
                return P,cumu#,trans
            else:
                return B,cumu#,trans
        
    def random_vector(self):
        return np.random.rand(self.number_of_nodes)
        
    def matrix_generation(self,generator,size,**prms):
        for i in range(size):
            self.append(nx.to_scipy_sparse_matrix(generator(**prms)))
 
        
                      
def cdf2histogram(c_in):
    """ Reads cdf and returns histogram. """
    if isinstance(c_in,list):
        c=c_in
    else:
        c=loadtxt("Shortest_Path_cdf.txt",dtype=int,usecols=(1,))
        
    h=[]
    h.append(c[0])
    for i in range(1,len(c)):
        h.append(c[i]-c[i-1])
    return h

if __name__=="__main__":
    #Z=ProductOfAdjacencyMatrices(nx.fast_gnp_random_graph,n=100,p=0.01,directed=True)
    At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,2))
    #At = AdjMatrixSequence(fs.dataPath("nrw_edges_01JAN2008_31DEC2009.txt"))
    #At=AdjMatrixSequence("Data/sociopatterns_hypertext_social_ijt.dat")
    #At=AdjMatrixSequence("Data/sexual_contacts.dat")
    #At.as_undirected()
    #At = AdjMatrixSequence(fs.dataPath("D_sw_uvd_01JAN2009_31MAR2010.txt"),matr_type='dok')
    #C=At.cumulated()
    
    #print len(At),At[0].shape

    print 'Matrixsequenz eingelesen'
    Z=ProductOfAdjacencyMatrices(At)
    print 'Produkt-Objekt erzeugt'
        
    P,c=Z.full_product_matrix(return_path_matrix=True,return_transitivity=False)
    h=cdf2histogram(c)
    
    print 'Schreibe', P.nnz
    gwh.dict2file(c,"Cumu_edges.txt")
    gwh.dict2file(h,"histo.txt")
    #gwh.dict2file(t,"Transitivity.txt")
    
    #out=P.sum(1)
    #inn=P.sum(0)
    #mmwrite("Vir.mtx",out)
    #mmwrite("Vul.mtx",inn)
    mmwrite("sexual_PathMatrix.mtx",P)
    """try:
        mmwrite("Cumulated.mtx",Cumu)
    except:
        savemat("Cumulated.mat",{'C':Cumu})"""




