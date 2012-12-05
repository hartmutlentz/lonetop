#! /usr/bin/python
#
import sys,os
sys.path.append('/Users/lentz/Documents/Physik/my_py') # Mac
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))
import tools.filesystem as fs
import tools.Gewindehammer as gwh

from MatrixList import AdjMatrixSequence
import scipy.sparse as sp
from scipy.io import mmread,mmwrite,loadmat,savemat



class LeslieSequence(list):
    """ Like AdjMatrixSequence, but with focus on products and spectra.
        Also contains creation of random matrices.
        
        Input: AdjMatrixSequence or networkx Graph Generator with parameters in args.
        
        ----------------
        Examples:        
        At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,3))
        Z=LeslieSequence(At)
        
    """
    def __init__(self,matrices,tau=3):
        assert tau<len(matrices), "Infectious period greater than observation time."
        list.__init__(self)
        
        self.extend(matrices)
        self.__transpose()
        
        self.number_of_nodes=self[0].shape[0]
        self.age_classes=tau
        self.max_dimension=self.age_classes*self.number_of_nodes
        
        self.shape_matrices=self.__get_shape_matrices()

    def __transpose(self):
        """ Transpose all matrices in self """
        for i in range(len(self)):
            self[i]=self[i].transpose()
        return
    
    def __get_shape_matrices(self):
        # Leslie age-based matrix shapes
        eye=sp.eye(self.age_classes,self.age_classes,-1)
        top=sp.lil_matrix((self.age_classes,self.age_classes),dtype='bool')
        for j in range(self.age_classes):
            top[0,j]=True
        return top,eye

    def bool_int_matrix(self,M):
        """ Returns matrix with only np.int64: ones. """
        M=M.astype('bool')
        M=M.astype('i')

    def leslie_matrix(self,M):
        """ Returns the Leslie Matrix of a given Matrix M """
        flux,aging=self.shape_matrices
        
        return sp.kron(flux,M)+sp.kron(aging,sp.identity(self.number_of_nodes))

    def initial_state(self,nodes=None):
        #
        if nodes:
            x=sp.lil_matrix((self.max_dimension,1),dtype='int')
            for i in nodes:
                assert i<self.number_of_nodes, 'Initial node not in node set.'
                x[i,0]=1
        else:
            x=sp.lil_matrix((self.max_dimension,1),dtype='int')
            # set the first age class to 1
            for i in range(len(self)):
                x[i,0]=1
                    
        return x

    def oldest_in_x(self,s):
        """ the number of non-zeros for the last self.number_of_nodes elements of s """
        indices=range(self.max_dimension-self.number_of_nodes,self.max_dimension)
        return s[indices,0].nnz
    
    def outbreak(self,initial_nodes=None):
        """ an outbreak """
        state=self.initial_state(initial_nodes)
        outbreak_size=[state.nnz]
        recovered=[0]
        print 'init state\n',state
        
        for i in range(len(self)):
            state=self.leslie_matrix(self[i])*state
            outbreak_size.append(state.nnz)
            print 'Step ',i+1,'\n',state
            new_recovered=self.oldest_in_x(state)
            old_recovered=recovered[-1]
            recovered.append(old_recovered+new_recovered)

            if recovered==self.number_of_nodes: break
            if state.nnz==0: break
        # if no steady state reached
        #recovered.append(state.nnz-self.oldest_in_x(state))
        
        return outbreak_size, recovered

                
    def unfold_accessibility(self,return_accessibility=False):
        #
        assert False, "To be done."
        
        P=self[0].copy()
        cumu=[0]

        for i in range(1,len(self)):
            print 'unfolding accessibility',i,'non-zeros: ',P.nnz
            cumu.append(P.nnz)
            P=P*self[i]
        
        if return_path_matrix:
            P = P.astype('bool')
            P = P.astype('int')    
            return P,cumu
        else:
            return cumu
 
        

if __name__=="__main__":
    At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,3),directed=True)
    #At = AdjMatrixSequence(fs.dataPath("nrw_edges_01JAN2008_31DEC2009.txt"))
    #At=AdjMatrixSequence("Data/sociopatterns_hypertext_social_ijt.dat",directed=False)
    #At=AdjMatrixSequence("Data/sexual_contacts.dat")
    #At.as_undirected()
    #At = AdjMatrixSequence(fs.dataPath("D_sw_uvd_01JAN2009_31MAR2010.txt"),matr_type='dok')
    print 'Matrixsequenz eingelesen',len(At)

    L=LeslieSequence(At,tau=3)
    infected,recovered=L.outbreak(initial_nodes=(1,))
    print infected,recovered
        
    #c=Z.unfold_accessibility(return_accessibility=False)
    #h=gwh.cdf2histogram(c)
    #gwh.dict2file(c,"Randomized-sexual/"+infile+"_Cumu_edges.txt")
    #gwh.dict2file(h,"Randomized-sexual/"+infile+"_histo.txt")






