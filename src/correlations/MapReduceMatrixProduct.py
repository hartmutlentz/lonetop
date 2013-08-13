#!/usr/bin/env python
#
#
import scipy.sparse as sp
#
# ----------------- ----------------- ----------------- -----------------
#
def map1(M,N):
    """ first map function for matrix-matrix multiplication.
        True and False in the values are Matrix labels and could be
        alternatively called 'M' and 'N', respectively.
        Considered as Boolean for memory reasons.
    """
    m_hash=dict([(j,(True,i,M[i,j])) for (i,j) in M.keys()])
    n_hash=dict([(j,(False,k,N[j,k])) for (j,k) in N.keys()])
    
    return m_hash,n_hash

def reduce1(m,n):
    """ Reads 2 dicts from map1 and returns {(i,k):[m_ij*n_jk,...]} for different j. """
    new={}

    # Matrix M
    for j in n:
        try:
            _,i,m_val = m[j]
            _,k,n_val = n[j]
            try:
                new[(i,k)].append(m_val*n_val)
            except KeyError:
                new[(i,k)]=[m_val*n_val]
        except KeyError:# do nothing, if entries are 0.
            pass

    # Matrix N
    for j in n:
        _,i,m_val = m[j]
        _,k,n_val = n[j]
        try:
            new[(i,k)].append(m_val*n_val)
        except KeyError:
            new[(i,k)]=[m_val*n_val]

    return new

def map2(x):
    """ Identity """
    return x

def reduce2(d,as_boolean=True):
    """ Reads dict d with keys (i,j) and values [float,float,...].
        Returns dict {(i,j):sum(d[(i,j)]),...}
    """
    if as_boolean: return dict([(k,any(d[k])) for k in d])
    else: return dict([(k,sum(d[k])) for k in d])

def dict_to_dok(d,dim,datatype=bool):
    """ Reads {(i,j):val,...} and converts to square matrix with dimension dim. """
    A=sp.dok_matrix((dim,dim),dtype=datatype)
    for (i,j) in d:
        A[i,j]=d[i,j]

    return A

def mr_matrix_product(M,N):
    """ Reads 2 dok-matrices and returns product matrix. """    
    m,n=map1(M,N)
    r1=reduce1(m,n)
    r2=reduce2(r1)
    
    return r2

if __name__=="__main__":
    A=sp.dok_matrix((3,3),dtype=bool)
    B=sp.dok_matrix((3,3),dtype=bool)
    A[0,1]=1
    A[0,2]=1
    B[1,2]=1

    #print A.keys()
    #a={(1,2):[3,4,5,6],(5,8):[],(12,5):[1,2,3]}
    x= mr_matrix_product(A,B)
    M=dict_to_dok(x,3)
    print type(M)
    

