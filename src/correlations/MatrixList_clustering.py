import sys, os
sys.path.insert(1, os.path.join(sys.path[0], os.pardir))
import tools.filesystem as fs

from scipy.sparse import coo_matrix, csr_matrix, lil_matrix, dok_matrix
import scipy.stats
import numpy as np
from numpy import loadtxt, zeros, savetxt
import scipy
from random import sample
from collections import defaultdict
import gc
import networkx as nx
from scipy.io import mmread,mmwrite,loadmat,savemat



class AdjMatrixSequence(list):
    """
    list of sparse matrixes

    class inherits from list -> has all properties of list
    the constructor expects filename of u,v,w,d-edgelist

    it constructs a list, where each entry is
    a sparse matrix (default: csr)
    the matrixes in the list are ordered by day

    the indices of the matrix represent the nodes, but
    nodes are reindexed to [0..number_of_nodes]
    """

    def __init__(self, edgelist_fname,write_label_file=False,columns=(0,1,2),firstday=None,lastday=None):
        list.__init__(self)
        #if edgelist_fname == fs.dataPath("D_sf_uvwd_cmpl.txt"):
        #    self.first_day = 2555 #2008-2010
        #else:
        #    self.first_day = 0#1825 #2006-2010
        #self.last_day = 9#3650
        self.first_day=firstday
        self.last_day=lastday
        self.fname = edgelist_fname
        self.cols=columns
        self.label_file = write_label_file
        
        self.matricesCreation()
        self.number_of_nodes=scipy.shape(self[0])[0]

    def groupByDays(self, edges):
        """ returns list of tupels: [(d,[(u,v),...]),...] """
        dct = defaultdict(list)
        for u,v,d in edges:
            dct[d].append((u,v))
        dct_s = dict.fromkeys(range(0, self.last_day-self.first_day), [])
        for d in dct:
            dct_s[d-self.first_day] = dct[d]
        return dct_s.items()

    def reindex(self, edges):
        """ for any index in edges returns dict with new_indices
            [0...number_of_unique_indices]
        """
        us, vs, ds = zip(*edges)
        nodes = set(us) | set(vs)
        old_to_new = {}
        for i,n in enumerate(nodes):
            old_to_new[n] = i
        return old_to_new
        
    def all_time_windows(self,max_window=None):
        """ summation over all time windows """
        p_corr={}
        if not max_window: mw=len(self)
        
        for twindow in range(mw):
            print twindow
            p_corr[twindow]=self.single_time_window(twindow)
        return p_corr
        
    def single_time_window(self,windowsize):
        """ summation over a time window """
        prop_corr=[]
        
        for i in range(len(self)-windowsize):
            prop_corr.append(self.two_link_density(i,i+windowsize))
            
        return (scipy.stats.scoreatpercentile(prop_corr,25),\
            scipy.stats.scoreatpercentile(prop_corr,50),\
            scipy.stats.scoreatpercentile(prop_corr,75)\
            )
        
    def two_link_density(self,index1,index2,norm=True):
        """ the link density for 2 step paths """
        C=self[index1]*self[index2]
        
        nlinks=C.sum()
        n=self.number_of_nodes
        
        if norm:
            return float(nlinks)/(float((n-2)*(n**2-n))+float(n**2-n))
        else:
            return float(nlinks)/float((n**2-n))
        
    def matrix_density(self,A):
        # density of a matrix, weights are ignored.
        n=float(scipy.shape(A)[0])
        e=float(A.nnz)
        return e/(n*(n-1.0))        
    
    def density(self):
        # densities of all matrices as list. weights are ignored.
        dens=[self.matrix_density(a) for a in self]
        return dens

    def deep_product(self,twindow=1,start=0):
        """ Product A_1*A_7*A_14... """
        C=self[start].copy()
        links={}
        links[start]=(C*C).sum()
        
        for i in range(start+twindow,len(self)-twindow,twindow):
            C=C*self[i]
            links[i]=C.sum()
            if C.nnz==0: break
        
        return links
        
    def greedy_product_path(self,start=0):
        """ Performs a greedy product path from start time.
            This path is the sequence that maximizes the 2-graph edge density,
            if memory effects are neclected.
        """
        max_days=[]
        pos=start
        day_density={}
        
        while pos<len(self):
            dichte=[0 for i in range(pos+1)]
            for i in range(pos+1,len(self)):
                dichte.append((self[pos]*self[i]).sum())
            
            if max(dichte)==0: break
            max_days.append(dichte.index(max(dichte)))
            pos=max_days[-1]
            
            day_density[pos]=max(dichte)
            
        return day_density
    

    def nx_2correlation_tree(self):
        """ returns acyclic networkx DiGraph with nodes as days and
            2-correlations as edgeweights
        
        """
        G=nx.DiGraph()
        G.add_nodes_from([i for i in range(len(self))])
        
        for i in range(len(self)):
            print i
            for j in range(i,len(self)):
                x=(self[i]*self[j]).sum()
                if x>0:
                    attr={}
                    attr['corr']=x
                    G.add_edge(i,j,attr)
        return G
        
    def cumulated(self):
        """ Returns Cumulted Graph as Matrix """
        C=csr_matrix((self.number_of_nodes, self.number_of_nodes), dtype=np.int64)
        for matrix in self:
            C = C+matrix
        return C
        
    def daily_activity(self):
        """ Dict {day:matrix_density} """
        da={}
        n=float(self.number_of_nodes)
        norma=n*(n-1.0)
        
        for i in range(len(self)):
            da[i]=float(self[i].nnz)/norma
        
        return da      

    def symmetrize_matrix(self,A):
        """ Returns symmetric version of a non-symm Matrix A as bool-int. """
        M = A + A.transpose()
        M = M.astype('bool')
        M = M.astype('float')
        return M
        
    def as_undirected(self):
        """ makes every matrix in self symmetric. """
        for i in range(len(self)):
            self[i]=self.symmetrize_matrix(self[i])

    def clustering_matrix(self,limit=None):
        """ Computes the matrix of clustering coefficients of a matrix sequence.
        
        """
        if limit:
            n=limit
        else:
            n=len(self)
        anzahl=n*(n-1)*(n-2)/6    
        C=lil_matrix((n,n),dtype='float')
        
        # Coefficients
        for i in range(n):
            print "i=",i, " of ",n
            for j in range(i+1,n):
                for k in range(j+1,n):
                    a3=self[i]*self[j]*self[k]
                    clu=(a3.diagonal()).sum()
                    clu_norm=a3.sum()
                    if clu_norm>0.0:
                        #print clu,clu_norm
                        C[j-i,k-j] += clu/clu_norm
        return C

    def matricesCreation(self):
        """ creates list of sparse matrices from input file """
        edges = loadtxt(self.fname, dtype=int, usecols=self.cols)
        
        # first and last days
        dummi, auchdummi, days = loadtxt(self.fname, dtype=int, usecols=self.cols, unpack=True)
        if not self.first_day:
            self.first_day=min(days)
        if not self.last_day:
            self.last_day=max(days)
        
        # use only days between FIRSTDAY and LASTDAY
        edges = [(u,v,d) for u,v,d in edges if (d>=self.first_day) and (d<=self.last_day)]
        
        # get dictionary of new indices and write map-file
        re_dct = self.reindex(edges)
        if self.label_file:
            g=file('oldindex_matrixfriendly.txt','w+')
            for k in re_dct:
                g.writelines(( str(k)+'\t'+str(re_dct[k])+'\n'))
            g.close        
        
        # reindex using this dictionary
        edges = [(re_dct[u],re_dct[v],d) for u,v,d in edges]

        edges = self.groupByDays(edges)

        # the actual construction of the sparse matrices
        mx_index = len(re_dct)
        for d, es in edges:
            us = [u for u,v in es]
            vs = [v for u,v in es]
            bs = [True for i in range(len(es))]

            m = csr_matrix((bs,(us,vs)), shape=(mx_index, mx_index), dtype=np.int32)
            self.append(m)
            
            
if __name__ == "__main__":
    from pprint import pprint
    import tools.Gewindehammer as gwh
    #At = AdjMatrixSequence(fs.dataPath("nrw_edges_01JAN2008_31DEC2009.txt"))
    At=AdjMatrixSequence("Data/sociopatterns_hypertext_social_ijt.dat")
    At.as_undirected()
    #At = AdjMatrixSequence(fs.dataPath("T_edgelist.txt"),columns=(0,1,3))
    #At = AdjMatrixSequence(fs.dataPath("D_sw_uvd_01JAN2009_31MAR2010.txt"))
    print 'Alle: ',len(At)
    C=At.clustering_matrix(limit=2000)
    mmwrite("Clusteringmatrix.mtx",C)
    #Cumu=At.cumulated()
    #try:
    #    mmwrite("Cumulated.mtx",Cumu)
    #except:
    #    savemat("Cumulated.mat",{'C':Cumu})
    
    #x=At.daily_activity()
    #gwh.dict2file(x,"daily_activity.txt")
    
    
