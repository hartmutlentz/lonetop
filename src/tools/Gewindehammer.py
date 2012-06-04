#!/usr/bin/env python2.3
#
#  untitled.py
#  
#
#  Created by Hartmut Lentz on 12.08.08.
#  Copyright (c) 2008 __MyCompanyName__. All rights reserved.
#
# Packages (Mac OS X): /Library/Python/2.5/site-packages
import math,random,string,pprint,numpy,scipy,scipy.sparse
from scipy.linalg import norm as scipynorm

import networkx as nx
try:
    from my_pool import ListPool
except:
    pass

class MessReihe():
    """ MessReihe is an environment for Multi-CPU computations of a
        function func. The parameter intervall is given at init.
        
        INPUT: function that returns tuple (input,output).
        Output can be anything.
    
    """    
    def __init__(self,func,xmin=0.0001,xmax=0.0009,dx=0.0001):
        assert xmin<=xmax, "x_max has to be greater than x_min."
        
        self.x_min=xmin
        self.x_max=xmax
        self.delta_x=dx
        self.the_function=func
        self.__tested=False
        
        self.input_values=self.parameter_values()
        self.__the_results={}

    def parameter_values(self):
        # get input values
        inputs=[]
        x=self.x_min
        while x<self.x_max:
            inputs.append(x)
            x += self.delta_x
        return inputs

    def single_run(self):
        # run on single CPU
        x=self.x_min
        result={}
        while x<self.x_max:
            result[x]=self.the_function(x)
            x += self.delta_x
        self.__the_results=result
        return
            
    def results(self):
        # getter
        return self.__the_results
        
    def write_results(self,fname):
        # prints results to file
        dict2file(self.__the_results,fname)
                        
    def run(self):
        # run on multi CPU
        listpool=ListPool(self.the_function,self.input_values)
        listpool.run()
        x=listpool.results
        self.__the_results=x
        return
        
    def test(self):
        # Test the input function
        print 'Start testrun...'
        x=self.the_function(self.x_min)
        assert len(x)==2, 'Function must return 2-Tuple.'
        assert type(x)==tuple, 'Function must return 2-Tuple.'
        assert x[0]==self.x_min, 'Function must return input value.'
        assert len(x[1])>=1, 'Function result is None.'
        self.__tested=True
        print '...testrun successfull.\n'

def read_file(filename,delimiter='\t',datatype=numpy.str):
    """
        
    """
    try:
        data=numpy.loadtxt(filename,delimiter=delimiter,dtype=datatype)
        return data.tolist()
    except NameError:
        print 'Numpy.loadtxt not found. Will use other routine.'
        read_file2(filename,delimiter)

def read_file2(filename,delimiter='\t'):
    """ Alternativ: numpy.loadtxt(filename,dtype=numpy.int64) bzw.
        dt={'names':('time','infected'),'formats':(np.int64,np.float32)}
        dann: numpy.loadtxt(filename,dtype=dt)
    """
    
    # Liest die Datei ein
    file_tab=file(filename,"r")
    line=file_tab.readline()
    data=[]
    while line!="":
        row=[]
        for val in string.split(line[:-1],delimiter):
            if string.strip(val)!="":
                row.append(string.strip(val))
        data.append(row)    
        line=file_tab.readline()
    file_tab.close()            
    return data
    
def data_minus_1(data):
    """
    Special function: renames int nodes: (1,2,3,...,n) -> (0,1,2,...,n-1) 
    """
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j]=data[i][j]-1
            if data[i][0]<0:
                print "!!! Nodes with labels < 0 !!!"
                return
            if data[i][1]<0:
                print "!!! Nodes with labels < 0 !!!"
                return
    
def string2int_data(data):
    nr_lines=len(data)
    nr_columns=len(data[0])
    
    for i in range(nr_lines): # Umwandeln in integer
        for j in range(nr_columns):
            data[i][j]=int(data[i][j])

def string2float_data(data):
    nr_lines=len(data)
    nr_columns=len(data[0])
    
    for i in range(nr_lines): 
        for j in range(nr_columns):
            data[i][j]=float(data[i][j])

def extract_col(data,column):
    """
    Extracts column from Data-array and returns it as list
    """
    data2=[]
    for i in range(len(data)):
        data2.append(data[i][column-1])
    return data2
    
def txt_files(searchstring='.txt',sub_dir="Smoothing_edges/05_2/"):
    """ Returns list of strings of all files in subdir.  """
    fil=[]
    for subdir, dirs, files in os.walk(sub_dir):
        fil.append(files)
    
    txt_files=[]
    for i in fil[0]:
        if searchstring in i: txt_files.append(i)
        
    paths=[]
    for i in txt_files:
        paths.append(sub_dir+i)
    
    return paths
    
def convert_scientific2float(inp,outp,sep):
    """
    Reads textfile and changes format from scientific to float, 
    i.e. 1.e-2 -> 0.01.
    Please change lines below for details.
    """
    
    data=read_file(inp,sep)
    string2float_data(data)

    stryx=[]
    for i in range(len(data)):
        strx= "%i" % data[i][0]     # change here
        stry= "%.14f" % (data[i][1])     # change here
        stryx.append((strx,stry))

    m=file(outp,'w+')
    for i in range(len(stryx)):
        m.writelines((stryx[i][0],'\t',stryx[i][1],'\n')) # change here
    m.close

def histogram(seq):
    dmax=max(seq)+1
    freq=[0 for d in range(dmax)]
    for d in seq:
        freq[d] += 1
    return freq
    
def highest_digraph_degree(G,degr='out',rtrn_type='node'):
    return maximum_degree(G,degr='out',rtrn_type='node')
    
def maximum_degree(G,degr='out',rtrn_type='tupel'):
    """
    Input: (Di)Graph G, 
    desired Degreetype (for DiGraphs): 'in' or 'out', 
    rtrn_type: 'dict','list','node'
    Output: Node with highest degree: {node:degree} or [node,degree] or node 
    """
    if not G.is_directed():
        ref=0
        degs=G.degree()
        for node in degs:
            if degs[node] > ref:
                ref=degs[node]
                target=node
    else:
        if degr=='out': outs=G.out_degree()
        elif degr=='in': outs=G.in_degree()
        else: 
            print 'Highest degree: Input incorrect!'
            return None
             
        ref=0
        for node in outs:
            if outs[node] > ref: 
                ref=outs[node]
                target=node
            
    if rtrn_type=='dict': return {target:ref}
    elif rtrn_type=='list': return [target,ref]
    elif rtrn_type=='node': return target
    elif rtrn_type=='tupel': return (target,ref)
    else:
        print 'Highest degree: Input incorrect!'
        return None
        
def longest_range_node(G):
    """ Returns node with the longest range of DiGraph G.
        Note that many nodes are possible, but only one is returned
    """
    x=ranges_single_sources(G)
    y=revert_dictionary(x)
    mm=max(y.keys())
    return y[mm][0]
    
def ranges(G,nodes=False,display=False):
    return ranges_single_sources(G,nodes,display)

def ranges_single_sources(G,nodes=False,display=False):
    """
    Space-saving version of ranges_small_graphs.
    Returns dict.
    """
    if nodes==False:
        nodes=G.nodes()
    if display:
        i=0
        no=G.number_of_nodes()
        rang={}
        for start in nodes:
            laenge=nx.single_source_shortest_path_length(G,start)
            rang[start]=len(laenge)-1
            print i,' von ',no
            i=i+1
    else:
        rang={}
        for start in nodes:
            laenge=nx.single_source_shortest_path_length(G,start)
            rang[start]=len(laenge)-1
    return rang

def node_range(G,node=False):
    if node:
        return range_single_node(G,node)
    else:
        return ranges_single_sources(G)

def range_single_node(G,node):
    """
    Space-saving version of ranges_small_graphs.
    Returns dict. *Use rang.get(i,' ') to avoid key errors*
    """
    laenge=single_source_shortest_path_length(G,node)
    rang=len(laenge)-1

    return rang    

def reachabilities(G,nodes):
    return reachabilities_single_sources(G,nodes)

def reachabilities_single_sources(G,nodes=None):
    """
    Number of nodes that can reach a node i
    Returns dict. 
    """
    if nodes==None:
        the_nodes=G.nodes()
    G1=G.reverse()
    return ranges_single_sources(G1,the_nodes)

def small_world_score(G,startnode,norm_by_range=True):
    """ returns value for histogram number_of_nodes(distance).
        Many close nodes increase this value.
    """
    laenge=shortest_path_length(G,startnode)
    del laenge[startnode]
    l2=revert_dictionary(laenge)
    sc=0.0
    
    for i in range(1,len(l2)):
        sc += len(l2[i])/float(i)

    if norm_by_range:
        if sc==0.0: return sc
        else: return sc/len(laenge)
    else: 
        return sc

def revert_dictionary(di,bijection=False):
    """
    Reverts dictionary. Output as dict of lists
    Input: {1:2, 5:2, 3:8,...}
    Output: {2:[1,5], 8:[3]}
    """
    rev={}
    if bijection:
        for x in di:
            rev[di[x]]=x
        if len(rev)!=len(di): 
            print 'Dictionary not bijective! Cannot revert.'
            return None
    else:
        for x in di:
            if di[x] in rev: rev[di[x]].append(x)
            else: rev[di[x]]=[x]
    
    return rev
    
def sort_list_by_length(lis):
    """
    returns list sorted by the lengths of entries
    [[1,2],[1],[3,4,5,6]] -> [[3,4,5,6],[1,2],[1]]
    """
    def cmp(a,b):
        return len(a)-len(b)
    lis.sort(cmp,reverse=True)
    
def file2dict(file,sep='\t',datatype='int',key_col=1,val_col=2):
    """ reads file and returns dict. 
        Datatypes: int, float or string
    """ 
    d=read_file(file,sep)
    if datatype=='int': string2int_data(d)
    elif datatype=='float': string2float_data(d)
    di={}
    for i in range(len(d)):
        di[d[i][key_col-1]]=d[i][val_col-1]
    
    return di
    
def write_sparse_matrix(A,nameoffile='matrix.txt'):
    # Write sparse matrix to textfile
    indices=zip(A.nonzero()[0],A.nonzero()[1])

    g=file(nameoffile,'w+')
    for k in indices:
        g.writelines((str(k[0]),'\t',str(k[1]),'\t',str(A[k]),'\n'))
    g.close
           
def dict2file(d,nameoffile='dict.txt',sorted=True):
    """ Writes dictionary (or list or tuple) to a textfile
        Sorted by keys, if sorted=True.
    """
    if scipy.sparse.isspmatrix(d):
        write_sparse_matrix(d,nameoffile)
        return
    
    def list2dict(li):
        x={}
        i=0
        for el in li:
            x[i]=el
            i+=1
        return x
    
    if not isinstance(d,dict): d=list2dict(d)
    
    dk=d.keys()
    if sorted: dk.sort()
    
    # if d={ 1: [a,b,c,...], 2:[d,e,f,...],... }
    s=d.values()[0]
    if isinstance(s,dict) or isinstance(s,list) or isinstance(s,tuple):
        laenge=len(d.values()[0])
        g=file(nameoffile,'w+')
        for k in dk:
            wstring=''
            for l in range(laenge): wstring += '\t'+str(d[k][l])
            g.writelines(( str(k)+wstring+'\n' ))
        g.close
        return
    
    g=file(nameoffile,'w+')
    for k in dk:
        g.writelines((str(k),'\t',str(d[k]),'\n'))
    g.close
    
    return

def sort_dict_by_values(d):
    """
    Returns list of tupels ordered by input values (increasing)
    Input: dictionary
    Output: [(a,9),(c,8),(y,2)]
    Note: if values are not unique, output is ordered by keys:
    [(1,9),(17,8),(2,4),(3,4)]
    """
    return sorted(d.items(), key=lambda(k,v):(v,k),reverse=True)

def sort_dict_by_keys(d,rev=False):
    """
    Input: A dictionary
    Output: list of tupels [(key,value),...] ordered by key 
    (increasing, if rev=False)
    """
    def ksort(d,rev):
        keys=d.keys()
        if not rev:
            keys.sort(reverse=False)
        else:
            keys.sort(reverse=True)
        return keys
    
    x=[(k,d[k]) for k in ksort(d,rev)]
    #x=[]
    #for k in ksort(d,rev):
    #    x.append((k,d[k]))
    return x
    
def list_zero_gap(liste):
    """
    returns (length of largest chain of connected 0s, frac(hi values)-frac(low values),single_solution).
    
    The second is one, if all values are large, 0 if 50/50 and -1 if all are small.
    single_solution is True, if there is only one longest chain in liste. 
    """
    if not 0 in liste:
        print 'No 0 in sequence. Will return (0,0.0,True)'
        return (0,0.0,True)
    if len(liste)<2:
        return (0,0.0,True)
    assert liste[-1]!=0, 'Largest value of range histogram cannot be zero.'
    
    def hi_lo_diff(seq,lo,hi):
        """
        for a sequence like [1,2,3,0,0,0,0,0,1,2]
        the difference between higher values and lower values is returned.
        higher and lower are seperated by the largest connected chain of 0s,
        i.e. hi and lo are the indices of the first and last 0 of the chain.
        """
        assert hi<len(seq), 'higher boundary > len(sequence).'
        assert lo<len(seq), 'lower boundary > len(sequence).'
        assert lo<hi, 'lower must be < higher!'
        
        lower=float(sum(seq[0:lo]))/sum(seq)
        higher=float(sum(seq[hi:]))/sum(seq)
        return higher-lower

    laengen={}    
    i=0
    single_maximum=True
    
    while i < len(liste):
        if liste[i]==0.0:
            startpt=i
            while liste[i]==0.0 and i<len(liste)-1: 
                i+=1
            endpt=i
            laengen[(startpt,endpt)]=endpt-startpt
        i+=1
        
    #print laengen
    laengen2=revert_dictionary(laengen)
    the_max_indices = laengen2[max(laengen.values())]
    
    # uniqueness
    thevector=scipy.array(laengen.values())
    uniqueness=angle_to_bisectrix(thevector)

    if len(the_max_indices)>1:
        print '---x--- Gap warning: Multiple maxima!'
        single_maximum=False
        l,h=the_max_indices[0][0], the_max_indices[0][1]
    else:
        l,h=the_max_indices[0][0], the_max_indices[0][1]
    hi_bias=hi_lo_diff(liste,l,h)

    return (max(laengen.values()),hi_bias,single_maximum)

def cdf(seq,rtn_type='dict',norm=True,bridge_values=False):
    """
    Input: Sequence as list. NO STRINGS!!! Convert to float or int first!!!
    Output: cdf (normalized). 
    rtn_type: 'dict', 'list_of_tupels' 
    -> {value,c(value),...} or [(value,frequency),...]
    bridge_values: adds values not present in data
    """
    def extr_tail(list,lower,upper):
        x=[]
        for i in range(lower,upper):
            x.append(list[i][1])
        return x
    
    # find frequencies
    cdf={}
    for i in range(len(seq)):
        if cdf.has_key(seq[i]) == True:
            cdf[seq[i]] += 1.0
        else:
            cdf[seq[i]] = 1.0
        
    # compute cdf
    cdf2=sort_dict_by_keys(cdf)
    
    cdf3={}
    for i in range(len(cdf2)):
        tail=extr_tail(cdf2,i,len(cdf2))
        #print len(tail)
        c_seq=sum(tail)
        wert=cdf2[i][0]
        
        cdf3[wert]=c_seq
            
    #normalize
    if norm:
        summe=sum(cdf3.values())
        for i in cdf3:
            cdf3[i]=cdf3[i]/float(summe)
    
    c=cdf3.values()
    if norm: print 'Conrol sum cdf: ',sum(c)
    
    if bridge_values:
        the_max=int(max(cdf3.keys()))
        the_last_one=1.0
        for i in range(the_max):
            if i in cdf3:
                the_last_one=cdf3[i]
            else:
                cdf3[i]=the_last_one           
    
    if rtn_type=='list_of_tupels':
        return sort_dict_by_keys(cdf3)
    elif rtn_type=='dict':
        return cdf3
    
def graph2dot(G,f_name='Graph.dot'):
    """ Write a (DiGraph) into a dot file <f_name>.dot
    """
    d=file(f_name,'w+')
    if G.is_directed(): d.writelines('digraph G {\n')
    else: d.writelines('graph G {\n')
    d.writelines('    size="6,6";\n')
    d.writelines('    node [shape=circle];\n')
    if G.is_directed():
        for ed in G.edges():
            d.writelines(('    ',str(ed[0]),' -> ',str(ed[1]),';','\n'))
    else:
        for ed in G.edges():
            d.writelines(('    ',str(ed[0]),' -- ',str(ed[1]),';','\n'))
    d.writelines('}\n')
    d.close()    
    

def data2dot_weighted(dat,f_name):
    # schreibt alle drei spalten der daten in eine dot-Datei
    d=file(f_name,'w+')
    d.writelines('digraph G {\n')
    d.writelines('    size="6,6"\n')
    d.writelines('    node [shape=circle]\n')
    for i in range(len(dat)):
        d.writelines(('    ',str(dat[i][0]),' -> ',str(dat[i][1]),' ',\
        '[ label = ',str(dat[i][2]),' ];','\n'))
    d.writelines('}\n')
    d.close()

def data2dot(dat,f_name,digraph=True):
    # schreibt die ersten zwei spalten der daten in eine dot-Datei (digraph oder graph)
    d=file(f_name,'w+')
    if digraph==True: d.writelines('digraph G {\n')
    else: d.writelines('graph G {\n')
    d.writelines('    size="6,6";\n')
    d.writelines('    node [shape=circle];\n')
    if digraph==True:
        for i in range(len(dat)):
            d.writelines(('    ',str(dat[i][0]),' -> ',str(dat[i][1]),';','\n'))
    else:
        for i in range(len(dat)):
            d.writelines(('    ',str(dat[i][0]),' -- ',str(dat[i][1]),';','\n'))
    d.writelines('}\n')
    d.close()
    
def edgefile2dot(edgefile,dotfile,sep,digraph=True):
    """
    Produces dot file from an edgefile.
    """
    komm=read_file(edgefile,sep)
    data2dot(komm,dotfile,digraph)


def hurwitz_zeta(alf,x_m,lim):
    z=[]
    for i in range(lim):
        z.append((i+x_m)**(-alf))
    return sum(z)
    
def riemann_zeta(alf,lim):
    if alf <= 1.0:
        print 'Error: Alpha < 1!'
        return 0
    z=[]
    for i in range(lim):
        z.append(i**(-alf))
    return sum(z)

def mle_alpha(seq,x_m):
    # returns alpha-1 from integer Input data
    # alpha-1 is the slope of the CDF
    fl_seq=[ 0 for i in range (len(seq)) ]
    for i in range(len(seq)):
        fl_seq[i]=float(seq[i])

    temp=[]
    for i in range(len(seq)):
        if fl_seq[i]>x_m:
            temp.append(math.log(fl_seq[i]/(x_m-0.5)))
        
    alpha=1.0+len(temp)/sum(temp)
    print 'Alpha-Abschaetzung: ',len(temp),' Werte wurden beruecksichtigt'

    return alpha-1.0

def print_graph_properties(G):
    print '\nEigenschaften des Graphen:'
    print 'Graph gerichtet?.....................',is_directed(G)
    print 'Anzahl der Knoten....................',G.number_of_nodes()
    print 'Anzahl der Kanten....................',G.number_of_edges()
    print 'Dichte...............................',link_density(G) 
    if G.is_directed():
        outs=G.out_degree()
        outs=outs.values()
        ins=G.in_degree()
        ins=ins.values()
        print 'maximaler In-Grad:...................',max(ins)
        print 'maximaler Aus-Grad:..................',max(outs)
        print 'Graph azyklisch?.....................',is_directed_acyclic_graph(G)
    else:
        grad=G.degree()
        grad=grad.values()
        print 'Maximaler Grad.......................',max(grad)
        #print clustering(G)
        #print 'Mittlerer Clustering-Koeffizient.....',average_clustering(G)
        co=is_connected(G)
        print 'Graph connected?.....................',co,'\n'
        if co and others:
            print 'Exzentrizitaet.......................',eccentricity(G)
            print 'Durchmesser..........................',diameter(G)
            print 'Radius...............................',radius(G)
            print 'Zentrum..............................',center(G)

def print_node_properties(G,node):
    print 'Properties of node #:..',node
    print 'Out-degree:............',G.out_degree(node)
    print 'In-degree:.............',G.in_degree(node)
    print 'Closeness:.............', nx.closeness_centrality(G,node,weighted_edges=True)
    print 'Range:.................', range_single_node(G,node),'\n'

def powerdev_pos_int_sequence(alf,size):
    """
    Power-law random generator from "Computer simulation of 
    attractors in stochastic models with 
    alpha-stable noise". 
    Input: Exponent between 0 and 2. This is the exponent of the CDF. 
    Return type: list with 'size' of positive integers. 
    """
    if alf>2.0:
        print 'exponent greater then 2!!!'
        return 0
    se=[]
    for i in range(size):
        y=random.random()
        v=y*math.pi-math.pi/2.0
        w=expovariate(1.0)
        x=math.sin(alf*v)/(cos(v))**(1.0/alf)*(cos(v-alf*v)/w)**((1.0-alf)/alf)
        x=int(abs(x)) 
        se.append(x)
    return se
    
def powerdev_pos_float(alf):
    """
    Returns a power-law distributed positive float.
    """
    if alf>2.0:
        print 'exponent greater then 2!!!'
        return None
    y=random.random()
    v=y*math.pi-math.pi/2.0
    w=expovariate(1.0)
    x=math.sin(alf*v)/(cos(v))**(1.0/alf)*(cos(v-alf*v)/w)**((1.0-alf)/alf)
    return abs(x) 

def matrix_friendly_node_labels(G,return_mapping=False):
    """ Returns (Di)Graph with nodes labelled 0,1,2,3,...
        Good for matrix representations.
        If return_mapping: returns ((Di)Graph,label_mapping)
    """
    n=G.number_of_nodes()
    new_label={}
    i=0
    
    for node in G.nodes():
        new_label[node]=i
        i += 1
    
    X=relabel_nodes(G,new_label)
    
    if return_mapping:
        return (X,new_label)
    else:
        return X

def w_func_graph2matrix(u,v):
    """ weight function for txt2sparse_matrix
        j is the line of the array dat.
    """ 
    freq=get_edge_data(u,v)['frequenz']
    vol=get_edge_data(u,v)['volumen']
    return vol*freq/912.5    

def graph2sparse_matrix(G,func=None, out='csr',conv_node_labels=False):
    """
    Returns sparse adjacency matrix from weighted or unweighted (Di)Graph G.
    Return type: default csr-Matrix. 'lil' possible.
    Note: Graph interprets weighted edges as distances (probably > 1.0), 
    Trade-Matrix could be < 1.0.
    If nodes of input graph are not 0,1,2,3,..., use conv_node_labels!
    """
    if conv_node_labels: G = matrix_friendly_node_labels(G)
    
    A_sp=scipy.sparse.csr_matrix((G.number_of_nodes(),G.number_of_nodes()), dtype='float32')
    
    if func:
        if G.is_directed():
            for i in G.edges(data=True):
                u, v = i[0], i[1]
                A_sp[int(u),int(v)]= func(u,v)        
        else:
            for i in G.edges(data=True):
                u, v = i[0], i[1]
                A_sp[int(u),int(v)]= func(u,v)
                A_sp[int(v),int(u)]= func(u,v)
    else:
        if G.is_directed():
            for i in G.edges(data=True):
                u, v = i[0], i[1]
                A_sp[int(u),int(v)]= 1.0        
        else:
            for i in G.edges(data=True):
                u, v = i[0], i[1]
                A_sp[int(u),int(v)]= 1.0
                A_sp[int(v),int(u)]= 1.0
    
    if out=='csr':
        return A_sp.tocsr()
    elif out=='lil':
        return A_sp.tolil()
    else: 
        print 'graph2sparse_matrix Input-Error.'
        return None

def graph2laplace_matrix(G,out='csr',conv_node_labels=False):
    """
    Returns sparse (out-)laplace matrix from unweighted (Di)Graph G.
    Return type: default csr-Matrix. 'lil' possible.
    Note: Graph interprets weighted edges as distances (probably > 1.0), 
    Trade-Matrix could be < 1.0.
    If nodes of input graph are not 0,1,2,3,..., use conv_node_labels!
    """
    if conv_node_labels: G = matrix_friendly_node_labels(G)
    
    A_sp=scipy.sparse.csr_matrix((G.number_of_nodes(),G.number_of_nodes()), dtype='int')
    
    if G.is_directed():
        for i in G.edges(data=True):
            u, v = i[0], i[1]
            A_sp[int(u),int(v)]= 1
        outd=A_sp.sum(axis=1)
        for i in range(G.number_of_nodes()):
            A_sp[i,i]=-outd[i,0]
    else:
        for i in G.edges(data=True):
            u, v = i[0], i[1]
            A_sp[int(u),int(v)]= 1
            A_sp[int(v),int(u)]= 1
    
    if out=='csr':
        return A_sp.tocsr()
    elif out=='lil':
        return A_sp.tolil()
    else: 
        print 'graph2sparse_matrix Input-Error.'
        return None

def w_func(dat,j,args):
    """ weight function for txt2sparse_matrix
        j is the line of the array dat.
    """ 
    frequency=float(dat[j][args['frequenz']-1])
    volume=float(dat[j][args['volumen']-1])
    return volume*frequency/912.5        

def txt2sparse_matrix(file,dim,funk=w_func,delimiter='\t',scale=1.0,\
    rtrn_type='lil',symm=False,selfloops=False,**weightcol):
    """
    Converts textfile to sparse matrix. Float scale is multiplied with matrix.
    weightcol: Weight column(s) as dict, counted 1,2,3,...
    example: txt2sparse_matrix(textfile,n,frequenz=3,volumen=4) 
    EXTERNAL FUNCTION w_funk NEEDED!
    """
    
    data=read_file(file,delimiter)
    string2int_data(data)
    
    A_sp=scipy.sparse.csr_matrix((dim,dim), dtype='float32')
    if symm:
        if weightcol:
            for i in range(len(data)):
                A_sp[data[i][0],data[i][1]]=funk(data,i,weightcol)
                A_sp[data[i][1],data[i][0]]=funk(data,i,weightcol)
        else:
            for i in range(len(data)):
                if data[i][0]!=data[i][1]:
                    A_sp[data[i][0],data[i][1]]=1.0
                    A_sp[data[i][1],data[i][0]]=1.0

    else:
        if weightcol:
            for i in range(len(data)):
                A_sp[data[i][0],data[i][1]]=funk(data,i,weightcol)
        else:
            for i in range(len(data)):
                    A_sp[data[i][0],data[i][1]]=1.0
    
    A_sp *= scale
    
    if not selfloops: 
        for i in range(dim): A_sp[i,i]=0.0
    
    if rtrn_type=='csr': return A_sp.tocsr()
    elif rtrn_type=='lil': return A_sp
    else: return None    

def txt2Graph(infile,delimiter='\t',\
    type='DiGraph',selfloops=False,convert2numeric='to_int',\
    attr=False):
    """
    like txt2Graph, but attributes are given to edges.
    Example: 
        if the 3rd column contains the weight of the edges:
        attr={'weight':3}
    Access to edge data: G.get_edge_data(u,v)['weight']
    """
    # read input file as data
    data=read_file(infile,delimiter)
    
    if convert2numeric=='to_int': string2int_data(data)
    elif convert2numeric=='to_float': 
        string2float_data(data)
        for i in range(len(data)):
            data[i][0]=int(data[i][0])
            data[i][1]=int(data[i][1])
    
    # initiate Graph
    if type=='DiGraph': G=nx.DiGraph()
    elif type=='Graph': G=nx.Graph()
    elif type=='MultiDiGraph': G=nx.MultiDiGraph()
    elif type=='MultiGraph': G=nx.MultiGraph()
    else:
        raise Exception, 'txt2Graph_attributed Input Error!'

    # add edges to graph
    if attr:
        for k in range(len(data)):
            u=data[k][0]
            v=data[k][1]
            attr_single_edge={}
            for key in attr: 
                attr_single_edge[key]=data[k][int(attr[key])-1] 
            G.add_edge(u,v,attr_single_edge)
    else:
        for k in range(len(data)):
            u=data[k][0]
            v=data[k][1]
            G.add_edge(u,v)

    if not selfloops: remove_selfloops(G)
    return G
    
def graph2file(G,filen='edgelist.txt'):
    """
    Writes edgelist from graph G into textfile.
    Note that isolated nodes are ignored.
    """
    g=file(filen,'w+')
    for e in G.edges():
        g.writelines((str(e[0]),'\t',str(e[1]),'\n'))
    g.close
    
    return

def graph_with_node_sizes(infile,sep='\t',wcol=3,directn='both'):
    """
    Returns (Di)Graph with edge attributes volume and
    node attributes size(=max(volume)).
    
    node sizes are estimated by returning their maximum
    in or out flow. 
    wcol - weight column (1,2,3,...) of input file.
    Direction for calculation: 'in', 'out' or 'both'.
    'both' returns undirected Graph.
    If node has no flux its size is defined as 1.0.
    """
    # change weight column here
    if directn=='both':
        G=txt2Graph(infile,sep,type='MultiGraph',to_int=False,attr={'vol':wcol})
    else: 
        G=txt2Graph(infile,sep,type='MultiDiGraph',to_int=False,attr={'vol':wcol})
        
    if directn=='in':
        for node in G.nodes():
            temp=[]
            for edge in G.in_edges(node):
                temp.append(G.get_edge_data(edge[0],edge[1])['vol'])
            if len(temp)>0: G.add_node(node,size=max(temp))
            else: G.add_node(node,size=1.0)
    elif directn=='out':
        for node in G.nodes():
            temp=[]
            for edge in G.out_edges(node):
                temp.append(G.get_edge_data(edge[0],edge[1])['vol'])
            if len(temp)>0: G.add_node(node,size=max(temp))    
            else: G.add_node(node,size=1.0)        
    elif directn=='both':
        for node in G.nodes():
            temp=[]
            for edge in G.edges(node):
                temp.append(G.get_edge_data(edge[0],edge[1])['vol'])
            if len(temp)>0: G.add_node(node,size=max(temp))
            else: G.add_node(node,size=1.0)
    else: return None
    
    return G
    
def simple_node_sizes(file,sep='\t',col=3):
    """
    Returns dictionary {node: size}.
    Input: MultiGraph Edgelist and weightcolumn col.
    Node size is computed as maximum of trade volume.
    """
    # read file
    x=read_file(file,sep)
    string2float_data(x)
    n=len(x)
    for i in range(n):
        x[i][0]=int(x[i][0])
        x[i][1]=int(x[i][1])
        
    # compute maximum flow
    nodes={}
    for i in range(n):
        print i
        if x[i][0] in nodes and nodes[x[i][0]]<x[i][col-1]:
            nodes[x[i][0]]=x[i][col-1]
        else:
            nodes[x[i][0]]=x[i][col-1]
            
        if x[i][1] in nodes and nodes[x[i][1]]<x[i][col-1]:
            nodes[x[i][1]]=x[i][col-1]    
        else:
            nodes[x[i][1]]=x[i][col-1]

    return nodes

def log_a(x,a):
    """
        Returns logarithm to real basis a. Log(x) means Ln(x) in Python!
    """
    return (math.log(float(x))/math.log(float(a)))
    
    
def generate_double_pl_digraph(nodes,edges,alpha,beta,isolated_nodes=True,display=False):
    """
    Generates a DiGraph with power-law in-degree (alpha) and power-law 
    out-degree (beta) distribution. Alpha and beta between 0 and 2. 
    First 2 arguments are number of nodes and edges. 
    """
    G=nx.DiGraph()
    if isolated_nodes:
        for j in range(nodes):  # ensures that all nodes (including isolated ones) are callable
            G.add_node(j)
    
    left=[]
    left.append(powerdev_pos_float(alpha))
    right=[]
    right.append(powerdev_pos_float(beta))    
    
    for i in range(1,nodes):
        left.append(left[i-1] + powerdev_pos_float(alpha))  
        right.append(right[i-1] + powerdev_pos_float(beta))

    maxleft=left[-1]
    maxright=right[-1]
                
    while G.number_of_edges() < edges:
        left_pointer=uniform(0,maxleft)
        right_pointer=uniform(0,maxright)
            
        i=0
        while left_pointer > left[i]: i=i+1
        left_pointer=i
        i=0
        while right_pointer > right[i]: i=i+1
        right_pointer=i
        
        G.add_edge(left_pointer,right_pointer)
        if display: print 'Edge ',G.number_of_edges(),' of ',edges
    return G

def angle_to_bisectrix(vect):
    """
    Returns the angle between vect and the bisectrix. Return 
    value is normalized to maximal value for positive half-space.
    """
    if sum(vect)==0.0: 
        return 0.0
    dim=len(vect)
    if dim==1: return 1.0
    phimax=math.acos(math.sqrt(1.0/dim))
    bisect=numpy.ones(dim,dtype='float32')
    term=numpy.dot(vect,bisect)/(scipynorm(vect)*scipynorm(bisect))
    if term <= 1.0:
        return math.acos(term)/phimax
    else:
        print 'angle_to_bisectrix. acos-term: ',term,\
        '. Will treat this as 1.0.'
        return math.acos(1.0)/phimax

def the_fle(vect):
    return the_famous_lentz_entropy(vect)

def the_famous_lentz_entropy(vect):
    return 1.0-angle_to_bisectrix(vect)

def inverse_participation_ratio(vect):
    """
    Returns the inverse participation ratio of a vector (distribution vector). 
    """
    M=sum(vect)**2
    m2=numpy.multiply(vect,vect)
    return float(sum(m2))/float(M)
    
def digraph_bidirectional(G):
    """
    Returns the fraction of edges that are bidirectional. 1 -> undirected 
    """
    #remove_selfloops(G)
    x=G.edges()
    temp=0.0
    for edge in G.edges():
        if G.has_edge(edge[1],edge[0]): temp += 1.0
    
    return temp/float(G.number_of_edges())
    
    
def pylab_histogramm():
    """
    The most important histogramm commands
    """
    import pylab as P
    hihi=P.hist(y.values(),1982)
    #pprint (hihi[0])
    #P.save('hihi.txt',hihi[0])
    #P.xlabel('Closeness')
    #P.ylabel('Anzahl')
    P.show()
    #P.savefig('histo.png')
    """ For matrices (not sparse!)"""
    #P.imshow(A)
    #P.hot()
    #P.colorbar()
    #P.show()
    
def pygraphviz_picture():
    """
    The most important pygraphviz commands
    """
    import pygraphviz as pgv
    A=pgv.AGraph('karate.dot')
    #print A.string() # print to screen
    A.layout('fdp') # layout with default (neato)
    #Layouts: twopi, gvcolor, wc, ccomps, tred, sccmap, fdp, 
    #circo, neato, acyclic, nop, gvpr, dot
    A.draw('simple.png') # draw png


def remove_selfloops(G):
    """
    Removes all selfloops of directed or undirected Graph G.
    """
    x=G.selfloop_edges()
    for edge in G.selfloop_edges():
        G.remove_edge(edge[0],edge[1])
        
def degree_vector(G):
    """
    Returns a numpy float vector containing the degrees of the nodes of G.
    If G is directed: returns tuple: (in_degree vector, out_degree vector).
    Nodes have to be integers.
    """
    if is_directed(G)==False:
        k=numpy.zeros(G.number_of_nodes(),dtype='float32')
        for i in G.nodes():
            k[i]=G.degree(i)
        return k
    else:
        k_in=numpy.zeros(G.number_of_nodes(),dtype='float32')
        k_out=numpy.zeros(G.number_of_nodes(),dtype='float32')
        for i in G.nodes():
            k_in[i]=G.in_degree(i)
            k_out[i]=G.out_degree(i)
        return (k_in,k_out)
        
def power_iteration(M,iterations=0):
    """
    Returns the dominant eigenpair (largest eigenvalue and eigenvector) 
    of Matrix M as a tuple.
    Iterates N steps where N is the dimension of M.
    """
    if iterations==0: iterations=3*shape(M)[0]
    
    def my_norm(x):
        return math.sqrt(add.reduce(real((conjugate(x)*x).ravel())))
        
    r=ones(shape(M)[0],dtype='float32')
    r[0]=0
    if my_norm(M*r)==0.0: return

    x_prev=r/my_norm(r)
    for i in range(1,iterations):
        y=M*x_prev
        x=y/my_norm(y)
        #control=my_norm(x-x_prev)
        eigval=dot(y,x_prev)
        x_prev=x.copy()
        print i
    return (eigval,x)
    
def rk4_step(input_vec,f_matrix,delta_t):
    """
    Runge Kutta 4th order step.
    INPUT:    Position (state) at time t
            ODE system as matrix
    OUTPUT: Position (state) at t+dt
    """
    k1 = f_matrix*input_vec * delta_t
    k2 = f_matrix*(input_vec + 0.5*k1) * delta_t
    k3 = f_matrix*(input_vec + 0.5*k2) * delta_t
    k4 = f_matrix*(input_vec + k3) * delta_t

    return input_vec + 1.0/6.0 *(k1 + 2.0*k2 + 2.0*k3 + k4)    

def rk4_1D(step_range=10,stepsize=0.1,x_init=0.1,t_init=0.0,reverse_output=False):
    """
    Differential equation: x' = f(x,t) -> please define derivs in subroutine derivs
    Initial condition: x(0)=x_0
    Computes a Runge-Kutta step (4th order) for function f with 2 variables and stepsize h.
    INPUT: above
    OUTPUT: list of pairs [(x,t),...]
    example: 
    >>> sol=rk4(10,0.1,0.1,0.0,True)
    >>> print sol
    """
    def derivs(x,t):
        """
        Input function for Runge-Kutta
        """
        dxdt = x- t**2.0 + 1.0 ############### change function here
        return dxdt    

    t=t_init
    x=x_init
    h=stepsize
    number_of_steps=int(step_range/stepsize)
    
    sol1=[]
    sol1.append((x,t))

    for i in range(number_of_steps):
        a=derivs(x,t)
        b=derivs(x+a*h/2.0,t+h/2.0)
        c=derivs(x+b*h/2.0,t+h/2.0)
        d=derivs(x+h*c,t+h)
        x = x + h/6.0*(a + 2.0*b + 2.0*c + d)
        t=t+h
        sol1.append((x,t))
    if reverse_output==True:
        sol2=[]
        for i in range(len(sol1)):
            sol2.append((sol1[i][1],sol1[i][0]))
        return sol2
    else: return sol1

def shift_21(inp,shiftno=1):
    """
    Shifts index from 0,...,n-1 to 1,...,n
    Input: list inp, int shiftno (default=1) for shifting to 2,3,...
    Output: dict
    Note: for-loops in C: for(i=1;i<=n;i++) -> for i in range(1,n+1)
    """
    outp={}
    if (shiftno==1):
        for i in range(len(inp)):
            outp[i+1]=inp[i]
    else:
        for i in range(len(inp)):
            outp[i+shiftno]=inp[i]
    return outp

def peakdetect(y_axis, x_axis = None, lookahead = 500, delta = 0):
    """
    Converted from/based on a MATLAB script at http://billauer.co.il/peakdet.html
    
    Algorithm for detecting local maxima and minmia in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maxima and minima respectively
    
    keyword arguments:
    
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- A x-axis whose values correspond to the 'y_axis' list and is used
        in the return to specify the postion of the peaks. If omitted the index
        of the y_axis is used. (default: None)
    lookahead -- (optional). x-window to consider around peak (measured
        in indices, not in dx!).
        Distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 500) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional). y-detection threshold.
        This specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the algorithm from picking up false peaks towards to end of
        the signal. To work well delta should be set to 'delta >= RMSnoise * 5'.
        (default: 0)
            Delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the algorithm
    
    returns -- two lists [maxtab, mintab] containing the positive and negative
        peaks respectively. Each cell of the lists contains a tupple of:
        (position, peak_value) 
        to get the average peak value do 'np.mean(maxtab, 0)[1]' on the results
    """
    maxtab = []
    mintab = []
    dump = []   # Used to pop the first hit which always if false
       
    length = len(y_axis)
    if x_axis is None:
        x_axis = range(length)
    
    # perform some checks
    if length != len(x_axis):
        raise ValueError, "Input vectors y_axis and x_axis must have same length"
    if lookahead < 1:
        raise ValueError, "Lookahead must be above '1' in value"
    if not (numpy.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"
    assert isinstance(lookahead,int), "lookahead must be integer"
    
    # needs to be a numpy array
    y_axis = numpy.asarray(y_axis)
    
    # maxima and minima candidates are temporarily stored in
    # mx and mn respectively
    mn, mx = numpy.Inf, -numpy.Inf
    
    # Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        #### look for max ####
        if y < mx-delta and mx != numpy.Inf:
            # Maxima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                maxtab.append((mxpos, mx))
                dump.append(True)
                # set algorithm to only find minima now
                mx = numpy.Inf
                mn = numpy.Inf
        
        #### look for min ####
        if y > mn+delta and mn != -numpy.Inf:
            # Minima peak candidate found 
            # look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                mintab.append((mnpos, mn))
                dump.append(False)
                # set algorithm to only find maxima now
                mn = -numpy.Inf
                mx = -numpy.Inf
    
    # Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            maxtab.pop(0)
            #print "pop max"
        else:
            mintab.pop(0)
            #print "pop min"
        del dump
    except IndexError:
        # no peaks were found, should the function return empty lists?
        pass
    
    return maxtab, mintab


def pol_interpolation(xa_in,ya_in,x):
    """
    Polynom interpolation. 
    Input: Data points and value x where to interpolate
    Output: f(x), error(f(x))
    """ 
    xa=shift_21(xa_in)
    ya=shift_21(ya_in)
    ns=1
    n=len(xa)
    n2=len(ya)
    if(n!=n2):
        print "Dimension Error."
        return
    
    dif = abs(x-float(xa[1]))
    c={}
    d={}
    
    for i in range(1,n+1):
        dift=abs(x-float(xa[i]))
        if dift < dif: 
            ns=i
            dif=dift
        c[i]=ya[i]
        d[i]=ya[i]
    y=ya[ns]
    ns=ns-1
    for m in range(1,n):
        for i in range(1,n-m+1):
            ho=xa[i]-x
            hp=xa[i+m]-x
            w=c[i+1]-d[i]
            den=ho-hp
            if (den == 0.0): 
                print "Error"
                #return
            den=w/den
            d[i]=hp*den
            c[i]=ho*den
        if (2*ns < (n-m)): dy=c[ns+1]
        else: 
            dy=d[ns]
            ns=ns-1
        y=y+dy
    return (y,dy)
    
def rat_interpolation(xa_in,ya_in,x):
    """
    Rational function interpolation. 
    Input: Data points and value x where to interpolate
    Output: f(x), error(f(x))
    """
    xa=shift_21(xa_in)
    ya=shift_21(ya_in)
    n=len(xa_in)
    tiny=1.0e-25
    ns=1
    c={} 
    d={}
    hh=abs(x-xa[1])
    for i in range(1,n+1):
        h=abs(x-xa[i])
        if(h==0.0):
            y=ya[i]
            dy=0.0
            return
        elif (h<hh):
            ns=i
            hh=h
        c[i]=ya[i]
        d[i]=ya[i]+tiny
    y=ya[ns]
    ns=ns-1
    for m in range(1,n):
        for i in range(1,n-m+1):
            w=c[i+1]-d[i]
            h=xa[i+m]-x
            t=(xa[i]-x)*d[i]/h
            dd=t-c[i+1]
            if(dd==0.0): 
                print "Pole at interpolation point"
                return
            dd=w/dd
            d[i]=c[i+1]*dd
            c[i]=t*dd
        if (2*ns<(n-m)): dy=c[ns+1]
        else: 
            dy=d[ns]
            ns=ns-1
        y=y+dy
    return (y,dy)
    
    
def spline_interpolation(x_in,y_in,yp1,ypn,x):
    """
    Cubic spline interpolation. 
    Input: 
    Data points (arrays): x_in,y_in. 
    2nd derivatives at first point (yp1) and last point (ypn). 
    x - point where to extrapolate
    Output: f(x)
    """
    def spline(x_in,y_in,yp1,ypn):
        """
        Subroutine for spline interpolation
        """
        n=len(x_in)
        x=shift_21(x_in)
        y=shift_21(y_in)
        u={}
        y2={}
        if (yp1>0.99e30): 
            y2[1]=0.0
            u[1]=0.0
        else:
            y2[1]=-0.5
            u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)    
        for i in range(2,n):
            sig=(x[i]-x[i-1])/(x[i+1]-x[i-1])
            p=sig*y2[i-1]+2.0
            y2[i]=(sig-1.0)/p
            u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
            u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
        if (ypn > 0.99e30): 
            qn=0.0
            un=0.0
        else:     
            qn=0.5
            un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))
        y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0)
        
        k=n-1
        while(k>0):
            y2[k]=y2[k]*y2[k+1]+u[k]
            k=k-1
        return y2

    xa=shift_21(x_in)
    ya=shift_21(y_in)
    y2a=spline(x_in,y_in,yp1,ypn)
    klo=1
    khi=len(x_in)
    while(khi-klo > 1):
        k=(khi+klo)/2
        if (xa[k]>x): khi=k
        else: klo=k
    h=xa[khi]-xa[klo]
    if (h==0.0):
        print "Spline Error"
        return
    a=(xa[khi]-x)/h
    b=(x-xa[klo])/h
    y=a*ya[klo]+b*ya[khi]+((a**3 - a)*y2a[klo]+(b**3 -b)*y2a[khi])*h**2 /6.0
    return y

def simple_statistics(data,display=False):
    """
    Returns average, average dev., standard deviation, 
    variance, skewness and kurtosis of a data list.
    """
    n=len(data)
    ep=0.0
    if (n<2):
        print "Error. Not enough data"
        return
    s=0.0
    for j in range(n): s=s+data[j]
    ave=s/float(n)
    adev=0.0
    var=0.0
    skew=0.0
    curt=0.0
    for j in range(n):
        s=data[j]
        adev=adev+abs(s-ave)
        ep=ep+s
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
    adev=adev/float(n)
    var=(var-ep*ep/float(n))/(float(n-1))
    sdev=math.sqrt(var)
    if display:
        if (var!=0.0):
            skew=skew/(float(n)*var*sdev)
            curt=curt/(float(n)*var*var)-3.0
            print "Average:\t",ave
            print "Av. deviation:\t",adev
            print "Std. deviation:\t",sdev
            print "Variance:\t",var
            print "Skewness:\t",skew
            print "Kurtosis:\t",curt 
            return (ave,adev,sdev,var,skew,curt)
        else:
            print "Average:\t",ave
            print "Av. deviation:\t",adev
            print "Std. deviation:\t",sdev
            print "Variance:\t",var
            print "Variance=0.0 -> no skew or kurtosis."
            return (ave,adev,sdev,var)
    else:
        if (var!=0.0):
            skew=skew/(float(n)*var*sdev)
            curt=curt/(float(n)*var*var)-3.0
            return (ave,adev,sdev,var,skew,curt)
        else:
            #print "Variance=0.0 -> no skew or kurtosis."
            return (ave,adev,sdev,var)

def ks_one(data):
    """
    KS-statistic d and probability prob are computed for given data (list) and
    CDF (subroutine func). (d is maximum distance between both cdf's).
    Returns: (d,prob). 
    Good match for p -> 1.
    """
    def func(x):
        # Cumulative distribution function ranging from 0 (smallest x) to 1 (largest x).
        return x/5.0
    def probks(alam):
        # Kolmogorov-Smirnov probability function
        eps1=0.001
        eps2=1.0e-8
        fac=2.0
        sum=0.0
        termbf=0.0
        a2= -2.0*alam*alam
        for j in range(1,101):
            term=fac*math.exp(a2*float(j)*float(j))
            sum=sum+term
            if (abs(term) <= eps1*termbf or abs(term) <= eps2*sum): return sum
            fac=-fac
            termbf=abs(term)
        return 1.0
        
    fo=0.0
    n=len(data)
    data.sort()
    en=n
    d=0.0
    for j in range(n):
        fn=float(j)/float(en)
        ff=func(data[j])
        dt=max(abs(fo-ff),abs(fn-ff))
        if (dt > d): d=dt
        fo=fn
    en=math.sqrt(en)
    prob=probks((en+0.12+0.11/float(en))*d)
    return (prob,d)
        
def ks_two(data1,data2):
    """
    KS-Test for two data sets.
    Returns: (p-value,D)
    D - maximum distance of the cdf's
    Good match for p->1
    """
    def probks(alam):
        # Kolmogorov-Smirnov probability function
        eps1=0.001
        eps2=1.0e-8
        fac=2.0
        sum=0.0
        termbf=0.0
        a2= -2.0*alam*alam
        for j in range(1,101):
            term=fac*math.exp(a2*float(j)*float(j))
            sum=sum+term
            if (abs(term) <= eps1*termbf or abs(term) <= eps2*sum): return sum
            fac=-fac
            termbf=abs(term)
        return 1.0

    n1=len(data1)
    n2=len(data2)
    j1=0
    j2=0
    fn1=0.0
    fn2=0.0
    data1.sort()
    data2.sort()
    en1=n1
    en2=n2
    d=0.0
    while (j1<n1 and j2<n2):
        d1=data1[j1]
        d2=data2[j2]
        if (d1<=d2): 
            fn1=float(j1)/float(en1)
            j1=j1+1
        if (d2<=d1):
            fn2=float(j2)/float(en2)
            j2=j2+1
        dt=abs(fn2-fn1)
        if (dt > d): d=dt
    en=math.sqrt(float(en1*en2)/float(en1+en2))
    prob=probks((en+0.12+0.11/float(en))*d)
    return (prob,d)
    
def fdf(x):
    # Example function for newton_raphson().
    f=x**2-2.0
    df=2.0*x
    return (f,df)

def newton_raphson(func, x1, x2, acc):
    """
    Returns root of function func within intervall x1 and x2 with accuracy acc.
    Input: x1,x2, func=(f(x),df/dx).
    (Works bad, if f'(x0) = 0)
    """
    jmax=20
    rtrn=0.5*(x1+x2)
    for j in range(jmax):
        (f,df)=func(rtrn)
        dx=f/df
        rtrn=rtrn-dx
        if ((x1-rtrn)*(rtrn-x2)<0.0):
            print "Error. Root not in brackets."
            return
        if (abs(dx)<acc): return rtrn
    print "Error. Maximum number of iterations exceeded."
    return
    
def graph_symmetry(G,subgraph=False):
    """
    Input: a directed graph
    Returns the number of edges of digraph G 
    that have no counterpart (and the subgraph (undirected)).
    """
    Gcp=G.copy()
    x=G.edges()
    lenx=len(x)
    xd=set()
    for i in x:
        xd.add(i)
    for edge in G.edges():
        if (edge[1],edge[0]) not in xd:
            Gcp.remove_edge(edge[0],edge[1])
            #x.remove((edge[0],edge[1]))
    new=Gcp.number_of_edges()
    frac=float(new)/float(lenx)/2.0
    if subgraph==False: return frac
    else: return (Gcp.to_undirected(),frac)

def machine_precision():
    # returns the machine precision
    eps=1.0
    while (1.0+eps !=1.0): eps *= 0.5
    return eps*2.0

def row_sum(M):
    """
    returns numpy vector containing the row sums of M
    """
    co=numpy.zeros(shape(M)[0],dtype='float32')
    ou=M.sum(axis=1)
    out=ou.tolist()
    for i in range(len(co)):
        co[i]=out[i][0]
    return co
    
def col_sum(M):
    """
    returns numpy vector containing the column sums of M
    """    
    co=numpy.zeros(shape(M)[0],dtype='float32')
    ou=M.sum(axis=0)
    out=ou.tolist()
    for i in range(len(co)):
        co[i]=out[0][i]
    return co
    
def txt2laplacian_matrix(file,delimiter,dim,symm=False,outward=True):
    """
    Converts textfile to sparse Laplacian matrix. 
    outward: for directed Laplacians: True for outdegree, False for indegree
    """        
    data=read_file(file,delimiter)
    string2int_data(data)
    
    # negative adjacency matrix
    A_sp=scipy.sparse.lil_matrix((dim,dim), dtype='float32')
    if symm==False:
        for i in range(len(data)):
            if data[i][0]!=data[i][1]: 
                A_sp[data[i][0],data[i][1]]=-1.0
    else:
        for i in range(len(data)):
            if data[i][0]!=data[i][1]:
                A_sp[data[i][0],data[i][1]]=-1.0
                A_sp[data[i][1],data[i][0]]=-1.0

    # degrees on diagonal
    if symm==False:
        out_deg=row_sum(A_sp)
        in_deg=col_sum(A_sp)
        if outward==True:
            for i in range(shape(A_sp)[0]):
                A_sp[i,i]=-out_deg[i]
        else:
            for i in range(shape(A_sp)[0]):
                A_sp[i,i]=-in_deg[i]
    else:
        deg=row_sum(A_sp)
        for i in range(shape(A_sp)[0]):
            A_sp[i,i]=-deg[i]

    return A_sp

        
def deg_av_nei_deg(G,out='deg-av_n_deg-err.txt'):
    """
    Produces txt-file: degree \t average_neighbor degree
    for undirected graphs only!
    """
    if G.is_directed():
        print 'Graph directed!'
        return None
    
    node_deg=G.degree(with_labels=True)
    deg_nodes=revert_dictionary(node_deg)
    deg_av={}
    #deg_er={}
    deg_er1={}
    deg_er2={}
    
    for degree in deg_nodes:
        nodes=deg_nodes[degree]
        temp=[]
        for node in nodes:
            for nei in G.neighbors(node):
                temp.append(G.degree(nei))
        deg_av[degree]=median(temp)
        deg_er1[degree]=quantile(temp,0.25,1)
        deg_er2[degree]=quantile(temp,0.75,1)
        #deg_er[degree]=std(temp)/math.sqrt(float(len(temp)))
        #print degree
    
    if 0 in deg_av: del deg_av[0]
    
    g=file(out,'w+')
    for d in deg_av:
        g.writelines((str(d),'\t',str(deg_av[d]),'\t',str(deg_er1[d]),\
        '\t',str(deg_er2[d]),'\n'))
    g.close        
    return

def directed_deg_av_nei_deg(G,out_d=True,nei_out=False):
    """
    Produces txt-file: degree \t average_neighbor degree
    for directed graphs only!
    out_d=True -> out_degrees on x-axis, out_d=False -> in-degrees on x
    nei_out=False -> neighbors' in-degrees on y; and vice versa 
    """
    if G.is_directed()==False:
        print 'Graph undirected!'
        return

    # output file name
    if out_d:
        if nei_out: outp='outd-av_n_outd.txt'
        else: outp='outd-av_n_ind.txt'
    else:
        if nei_out: outp='ind-av_n_outd.txt'
        else: outp='ind-av_n_ind.txt'
    
    if out_d: node_deg=G.out_degree(with_labels=True)
    else: node_deg=G.in_degree(with_labels=True)
    
    deg_nodes=revert_dictionary(node_deg)
    deg_av={}
    deg_er1={}
    deg_er2={}
    
    # for the right neighbors
    if out_d: G1=G
    else: G1=G.reverse()
    
    for degree in deg_nodes:
        nodes=deg_nodes[degree]
        temp=[]
        for node in nodes:
            if nei_out:
                for nei in G1.neighbors(node):
                    temp.append(G.out_degree(nei))
            else:
                for nei in G1.neighbors(node):
                    temp.append(G.in_degree(nei))                
    
        deg_av[degree]=median(temp) #mean(temp)
        deg_er1[degree]=quantil(temp,0.25,1) #std(temp)/math.sqrt(float(len(temp)))
        deg_er2[degree]=quantil(temp,0.75,1)
        #print degree
        
    if 0 in deg_av: del deg_av[0]
    
    g=file(outp,'w+')
    for d in deg_av:
        g.writelines((str(d),'\t',str(deg_av[d]),'\t',str(deg_er1[d]),\
        '\t',str(deg_er2[d]),'\n'))
    g.close        
    return
    
def degree_clustering_coefficient(G,out='deg_clust_coeff_err.txt'):
    """
    Produces txt file: deg \t clustering coefficient.
    For undirected graphs only!
    """
    cl=nx.clustering(G,with_labels=True)
    
    node_deg=G.degree(with_labels=True)
    deg_nodes=revert_dictionary(node_deg)    
    
    cluster_av={}
    cluster_er1={}
    cluster_er2={}
    
    for degree in deg_nodes:
        nodes=deg_nodes[degree]
        temp=[]
        for node in nodes: temp.append(cl[node])
        cluster_av[degree]=median(temp)
        cluster_er1[degree]=quantil(temp,0.25,1)
        cluster_er2[degree]=quantil(temp,0.75,1)
        #print degree
    
    #if 0 in deg_av: del deg_av[0]
    
    g=file(out,'w+')
    for d in cluster_av:
        g.writelines((str(d),'\t',str(cluster_av[d]),\
        '\t',str(cluster_er1[d]),'\t',str(cluster_er2[d]),'\n'))
    g.close        
    return
        

    
def maximum_weight_subgraph(G1,degree=2):
    """
    Returns a (Di)Graph containing only the highest 
    edges of G. Degree determines the out-degree of the nodes
    of returned Graph.
    """
    def remove_weakest_edges(G,node,weightname,howmany=degree):
        """removes edges with less weight from node.
           Leaves the strongest howmany edges. """
        eddes=[]
        for n in G.neighbors(node):
            wei=G.get_edge_data(node,n)[weightname]
            eddes.append((wei,n))
        eddes.sort(reverse=True)
        while len(eddes)>howmany:
            tokill=eddes.pop()
            G.remove_edge(node,tokill[1])    
    
    G=G1.copy()
    # get name of weight attribute
    sample1, sample2 = G.edges()[0][0], G.edges()[0][1]
    tdic=G.get_edge_data(sample1,sample2)
    weight=tdic.keys()[0]
    
    for node in G.nodes():
        remove_weakest_edges(G,node,weight,degree)
        
    return G


def generate_directed_ER_graph(nnodes,p_edge,rev=0.05,prnt=False):
    """
    Returns directed ER Graph.
    Input: nnodes: number of nodes,p_edge: expected density,
    rev: expected fraction of bidirectional edges
    """
    ER=nx.erdos_renyi_graph(nnodes,p_edge*2.0)
    #ER=fast_gnp_random_graph(nnodes,p_edge*2.0)
    G=ER.to_directed()
    
    # remove some bidirectional edges
    for edge in ER.edges():
        if random.random() > rev/2.0:
            if random.random()<0.5: G.remove_edge(edge[0],edge[1])
            else: G.remove_edge(edge[1],edge[0])

    # print some properties
    if prnt:
        print 'Density: ', link_density(G)
        print 'Fraction of bidirectional edges: ',digraph_bidirectional(G)
    
    return G

def directed_BA_graph(n=2000,m=2,p_r=0.5):
    """ Returns directed Barabasi-Albert-Graph
        with redirection probability p_r.
        n-number of nodes. m-number of old nodes
        a new one is connected to
    """
    G=nx.barabasi_albert_graph(n,m)
    D=G.to_directed()

    for ed in G.edges():
        pass
        if random.random()>p_r:
            x=[ random.choice(( (ed[0],ed[1]),(ed[1],ed[0]) )) ]
            D.remove_edges_from(x)    
    return D
    
def directed_graph(G,p_bi=0.1):
    """    
    Returns directed Graph
    with bidirection probability p_bi.
    """
    D=G.to_directed()

    x=[]
    for ed in G.edges():
        pass
        if random.random()>p_bi:
            x.append( random.choice(( (ed[0],ed[1]),(ed[1],ed[0]) )) )
    D.remove_edges_from(x)    
    return D


def generate_community_test_graph(p_in=0.05,p_out=0.001,mod_size=32,\
    mod_number=4, directed=False,files=False,rtrn_comdict=False,\
    print_info=False):
    """
    returns (un)directed Test Graph.
    Test Graph with 128(mod_size) nodes in 4(mod_number) modules.
    if files=true: 2 txt files: edgelist, {node:comm}
    If rtrn_comdict=true: returns ( (Di)Graph,{node:comm} ).
    Probability for bidirectional inner community edges is set to 0.05.
    
    Guimera, Amaral: Cartography of complex networks
    NOTE: p_out doesn't scale as p_in! It should be p_out<<1!
    """
    graphs=[]
    # node label mapping
    def new_node_labels(G,co):
        nn={}
        for node in G.nodes():
            nn[node]=node+co
        return nn
    
    # generate (directed) ER subgraphs and relabel nodes        
    if directed:
        for i in range(mod_number):
            G=generate_directed_ER_graph(mod_size,p_in)
            while not nx.is_connected(G.to_undirected()): #is_weakly_connected(G):
                print 'module not connected! Will try again. p_in too small?'
                G=generate_directed_ER_graph(mod_size,p_in)
            map=new_node_labels(G,mod_size*i)
            G=nx.relabel_nodes(G,map)
            graphs.append(G)
    else:
        for i in range(mod_number):
            G=nx.erdos_renyi_graph(mod_size,p_in)
            while not nx.is_connected(G):
                print 'module not connected! Will try again. p_in too small?'
                G=nx.erdos_renyi_graph(mod_size,p_in)
            map=new_node_labels(G,mod_size*i)
            G=nx.relabel_nodes(G,map)
            graphs.append(G)
        
    # union subgraphs
    if directed: X=nx.DiGraph()
    else: X=nx.Graph()
    
    for g in graphs:
        X=nx.union(X,g)        

    # fast connect graphs
    # number of edges between subgraphs
    nedges=p_out*float(mod_number)*float(mod_size)*float((X.number_of_nodes()-mod_size))
    if nedges<=float(len(graphs)-1):
        print 'Community TestGraph: p_out -> 0. Will return ring of modules.'
        nedges=float(len(graphs)-1)
            
    if not directed:
        # (weakly) connect minimaly: ring connection
        for i in range(len(graphs)-1):
            u, v = random.choice(graphs[i].nodes()), random.choice(graphs[i+1].nodes())
            X.add_edge(u,v)
        nedges -= float(len(graphs)-1)
        
        while nedges>0.0:
            g1, g2 = 0, 0
            while g1==g2:
                g1=random.choice(graphs)
                g2=random.choice(graphs)
            u,v = random.choice(g1.nodes()), random.choice(g2.nodes())
            X.add_edge(u,v)
            nedges -= 1.0
    else:
        # (weakly) connect minimaly: ring connection
        for i in range(len(graphs)-1):
            u, v = random.choice(graphs[i].nodes()), random.choice(graphs[i+1].nodes())
            if random.choice([True,False]): X.add_edge(u,v)
            else: X.add_edge(v,u)
        nedges -= float(len(graphs)-1)
        
        while nedges>0.0:
            g1, g2 = 0, 0
            while g1==g2:
                g1=random.choice(graphs)
                g2=random.choice(graphs)
            u,v = random.choice(g1.nodes()), random.choice(g2.nodes())
            if random.choice([True,False]): X.add_edge(u,v)
            else: X.add_edge(v,u)
            nedges -= 1.0

    # print information about graph
    ed=0
    for g in graphs:
        ed += g.number_of_edges()
    int_ed=X.number_of_edges()-ed
    if print_info:
        print '\n--- Test Graph ---'
        print 'Within module edges: ', ed
        print 'Inter module edges:  ', int_ed
        print 'Fraction of inter-module edges: ',\
        float(int_ed)/(float(ed)+float(int_ed))
        print 'My Density:', link_density(X),'\n'
    
    # node:community
    nodecom={}
    ii=0
    for g in graphs:
        for node in g.nodes(): nodecom[node]=ii
        ii+=1
    # output files
    if files:
        graph2file(X,'testgraph-edgelist.txt')
            
        y=file('testgraph-node-comm.txt','w+')
        for node in nodecom:
            y.writelines((str(node),'\t',str(nodecom[node]),'\n'))
        y.close    

    if rtrn_comdict: return (X,nodecom)
    else: return X
    
def generate_reciprocity_comm_test_digraph(p_bidirection=0.1,\
    p_in=0.05,p_out=0.005,mod_size=32,mod_number=4,comdict=False):
    """
    returns modular DiGraph. p_bidirection is probability for a link 
    to be bidirectional.
    """
    p_unidirection=1.0-p_bidirection
    
    if comdict:
        H,cd = generate_community_test_graph(p_in,p_out,\
        mod_size,mod_number,rtrn_comdict=comdict)
    else: 
        H=generate_community_test_graph(p_in,p_out,\
        mod_size,mod_number,rtrn_comdict=comdict)
    
    G=H.to_directed()
    #H=G.copy()
    
    for ed in H.edges():
        if random.random()<p_unidirection: 
            if random.random()<0.5: G.remove_edge(ed[0],ed[1])
            else: G.remove_edge(ed[1],ed[0])
        
    #print '\nReciprocity (Testgraph): ', reciprocity(G)[0]
        
    if comdict: return (G,cd)
    return G
    
def shells_single_root(G,root,rtn_stats=False):
    """
    Returns the shells for a startnode root.
    If rtn_stats False: {shell, degrees}
    If rtn_stats True: nodes in shell, 
    av. degree and sts_error, Degree sequence 
    is not returned in latter case.
    """
    
    if G.is_directed(): return None
    
    node_sh={}
    q=[root]
    sh=0
    # add nodes to shells
    while q:
        next=[]
        for node in q:
            if node not in node_sh: node_sh[node]=sh
            for nei in G.neighbors(node):
                if nei not in node_sh:
                    next.append(nei)
        q=[]
        sh += 1
        for node in next: q.append(node)
    # shell: nodes
    x=revert_dictionary(node_sh)
    
    # degrees and average in shell
    if rtn_stats:
        sh_degrees = {}
        sh_av_deg, sh_stderr = {}, {}
        for dist in x:
            sh_degrees[dist] = G.degree(x[dist])
            sh_av_deg[dist] = mean(sh_degrees[dist])
            sh_stderr[dist] = std(sh_degrees[dist])/math.sqrt(len(x[dist]))
        return (x,sh_av_deg,sh_stderr)
    else: 
        sh_degrees = {}
        for dist in x:
            sh_degrees[dist] = G.degree(x[dist])    
        return sh_degrees

def shells_all_roots(G):
    """
    Performs shell computation for all roots in G.
    Returns {shell: average degree}
    """ 
    if G.is_directed(): return None
    
    sh_alldegrees={}
    sh_stderr={}
    for root in G.nodes():
        x=shells_single_root(G,root)
        for shell in x:
            if shell in sh_alldegrees:
                sh_alldegrees[shell].extend(x[shell])
            else:
                sh_alldegrees[shell]=x[shell]
        print 'Root: ', root
        
    for sh in sh_alldegrees:
        sh_stderr[sh]=std(sh_alldegrees[sh])/math.sqrt(float(len(sh_alldegrees[sh])))
        sh_alldegrees[sh]=mean(sh_alldegrees[sh])
        
    return (sh_alldegrees,sh_stderr)
    
def link_density(G):
    # Density of (Di)Graph G
    n=float(G.number_of_nodes())
    e=float(G.number_of_edges())
    if G.is_directed(): return e/(n*(n-1.0))
    else: return 2.0*e/(n*(n-1.0))
    
def is_integer(x):
    # Returns True, if x is an integer
    return type(x)==int
    
def quantile(x, q, qtype = 4, issorted = False): 
    """ 
    NOTE: Alternative: scipy.stats.mstats.mquantiles(x)
    -> [q0.25, q0.5, q0.75]
    
    Given an array x, and a quantile value from from 0 to 1.0, 
    the function will return a value of x which may not be an 
    element of the array such that P(X <= x_q) < q, 
    i.e. q * 100 percent of the data are less than or equal to x_q.

    Args: 
    x - input data 
    q - quantile 
    qtype - algorithm 
    issorted- True if x already sorted.   

    Compute quantiles from input array x given q.
    For median, specify q=0.5.   

    References: 
        http://reference.wolfram.com/mathematica/ref/Quantile.html 
        http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile   

    Author: 
        Ernesto P.Adorio Ph.D. UP Extension Program in Pampanga, Clark Field. 
    """ 
    # one element arrays
    if len(x)<2: return x[0]
    
    if not issorted: y = sorted(x) 
    else: y = x 
    if not (1 <= qtype <= 9): 
        print 'Qtype=',qtype,' is not valid.'
        return None # error! 
    # int or float input
    if is_integer(y[0]) and qtype>3:
        print y[0],'Input must be float for methods 4--9.'
        return None
    elif not is_integer(y[0]) and qtype<4:
        print 'Input must be integer for methods 1--3.'
        return None  

    # Parameters for the Hyndman and Fan algorithm 
    abcd = [ # discontinous
            (0, 0, 1, 0), # inverse empirical distrib. function, R type 1 
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2 
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3   
            # continious
            (0, 0, 0, 1), # California linear interpolation, R type 4 
            (0.5, 0, 0, 1), # hydrologists method, R type 5 
            (0, 1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
            (1, -1, 0, 1), # mode-based method,(S, S-Plus), R type 7 
            (1.0/3, 1.0/3, 0, 1), # median-unbiased , R type 8 
            (3/8.0, 0.25, 0, 1) # normal-unbiased, R type 9. 
            ]   

    a, b, c, d = abcd[qtype-1] 
    n = len(x) 
    g, j = scipy.modf( a + (n+b) * q -1) 
    if j < 0: return y[0] 
    elif j > n: return y[n]   

    j = int(scipy.floor(j)) 
    if g == 0: return y[j] 
    else: return y[j] + (y[j+1]- y[j])* (c + d * g)   

def reciprocity(G):
    """
    See D. Garlaschelli, 2008: Patterns of link reciprocity in directed networks.
    For DiGraphs only!
    """
    remove_selfloops(G)
    r=float(digraph_bidirectional(G)) # fraction of bidir edges
    m=float(G.number_of_edges())
    n=float(G.number_of_nodes())
    den=m/(n*(n-1.0))
    l_bi=r*m # number of bi edges
    
    # rho and rho_min
    rho=(r-den)/(1.0-den)
    rho_min=-den/(1.0-den)
    
    # standard deviation
    rho_bi = ( (l_bi-2.0)/(m-2.0)-(m-2.0)/n/(n-1.0) )/( 1.0-(m-2.0)/n/(n-1.0) )
    rho_mo = ( l_bi/(m-1.0)-(m-1.0)/n/(n-1.0) )/( 1.0-(m-1.0)/n/(n-1.0) )
    s_rho = l_bi*(rho-rho_bi)**2 + (m-l_bi)*(rho-rho_mo)**2
    std_rho = math.sqrt(s_rho)
    
    return (rho, std_rho, rho_min)

def degree_assortativity_coefficient(G):
    """
    Computes assortativity coefficient for an 
    unweighted Graph or DiGraph.
    See: 
    M. Newman, Mixing patterns in networks, 2003
    """
    if not G.is_directed(): 
        G=G.to_directed() 
    
    m=1.0/float(G.number_of_edges())
    def j(G,i):
        return float(G.in_degree(i))-1.0
    def k(G,i):
        return float(G.out_degree(i))-1.0    
    
    a,b,c,d,e = [],[],[],[],[]
    for edge in G.edges():
        s, t = edge
        a.append(j(G,s)*k(G,t))
        b.append(j(G,s))
        c.append(k(G,t))
        d.append(j(G,s)**2)
        e.append(k(G,t)**2)
        
    r1=sum(a)-m*sum(b)*sum(c)
    r2=(sum(d)-m*sum(b)**2) * (sum(e)-m*sum(c)**2)
    
    return r1/math.sqrt(r2)
    
def within_module_degree(edgelist,node_comm_list):
    """
    Returns within module degree as z score for each node in network
    """
    G=txt2Graph(edgelist,type='Graph')
    noco=read_file(node_comm_list,'\t')
    string2int_data(noco)
    
    # { node:comm } dict
    node_comm={}
    comm_nodes={}
    for i in range(len(noco)): 
        node_comm[noco[i][0]]=noco[i][1]
        comm_nodes[noco[i][1]]=[]
    del noco

    # nodes in modules: { module:[nodes] }
    for node in node_comm: 
        comm_nodes[node_comm[node]].append(node)
    
    # within module degree of each node
    z={}
    for node in G.nodes():
        z[node]=0.0
        x=G.neighbors(node)
        for nei in x:
            if node_comm[nei]==node_comm[node]: z[node] +=1.0
    
    # average and std_dev of z in each module
    average={}
    stddev={}
    for commu in comm_nodes:
        seq=[]
        for x in comm_nodes[commu]: seq.append(z[x])
        average[commu]=simple_statistics(seq)[0]
        stddev[commu]=numpy.std(seq) # fuer Grundgesamtheit
        #simple_statistics(seq)[2] # fuer Stichproben-std

    # compute z_score and modify, if std_dev = 0
    z_score={}
    for node in z:
        if stddev[node_comm[node]] > 0.0: 
            z_score[node]= (z[node]-average[node_comm[node]])/stddev[node_comm[node]]
        else: z_score[node]=0.0

    return z_score
    
def participation_coefficient(edgelist,node_comm_list):
    """
    Returns participation coefficient for a modular network.
    For connected graphs only.
    Cluster enumeration has to start with 0 and mustn't have gaps!
    
    """
    G=txt2Graph(edgelist,type='Graph')
    noco=read_file(node_comm_list,'\t')
    string2int_data(noco)
    
    # { node:comm } dict
    node_comm={}
    comm_nodes={}
    for i in range(len(noco)): 
        node_comm[noco[i][0]]=noco[i][1]
        comm_nodes[noco[i][1]]=[]
    del noco

    # nodes in modules: { module:[nodes] }
    for node in node_comm: 
        comm_nodes[node_comm[node]].append(node)

    # inter degree matrix: A[node, comm]=number_of edges
    n_nodes=G.nodes()[-1] #len(node_comm)
    n_modules=len(comm_nodes)

    A_sp={}#scipy.sparse.lil_matrix((n_nodes,n_modules), dtype='float32')
    for node in G.nodes():
        x=G.neighbors(node)
        for nei in x:
            if (node,node_comm[nei]) in A_sp:
                A_sp[(node,node_comm[nei])] += 1.0
            else:
                A_sp[(node,node_comm[nei])] = 1.0
            
    # compute participation coefficient
    k=G.degree()
    p={}
    for i in G.nodes(): #range(n_nodes):
        p[i]=0.0
        #print 'node: ',i
        for j in comm_nodes.keys(): #range(n_modules): 
            if (i,j) in A_sp:
                p[i] += (A_sp[(i,j)]/float(k[i]))**2

        p[i] = 1.0 - p[i]
    
    return p
    
def modularity_Q(G,partition):
    """ Returns modularity for an undirected unweighted 
        graph G given a partition {node:community,...}.
    """
    m=float(G.number_of_edges())
    q=0.0
    
    # nodes in modules: { module:[nodes] }
    comm_nodes={}
    for i in partition:
        comm_nodes[partition[i]]=[]
    for node in partition: 
        comm_nodes[partition[node]].append(node)
    
    if G.is_directed(): 
        print 'for undirected graphs only!'
        return None
    else:
        for module in comm_nodes:
            S=G.subgraph(comm_nodes[module])
            lc=float(S.number_of_edges())
            nodedeg=G.degree(comm_nodes[module])
            dc=float(sum(nodedeg.values()))
            delta_q = lc/m - (dc/2.0/m)**2
            q += delta_q
    
    return q

def newman_modularity_Q(G,partition):
    """
    Returns modularity.
    Input: unweighted (Di)Graph and
    partition - dict { node:community,... }.
    From: Leicht, Newman, PRL 100, 2008
    """
    m=float(G.number_of_edges())
    n=G.number_of_nodes()
    q=0.0

    # nodes in modules: { module:[nodes] }
    comm_nodes={}
    for i in partition:
        comm_nodes[partition[i]]=[]
    for node in partition: 
        comm_nodes[partition[node]].append(node)
        
    if G.is_directed():
        for module in comm_nodes:
            nodes=comm_nodes[module]
            for node in nodes:
                for xnode in nodes:
                    if G.has_edge(node,xnode):
                        q += 1.0-float(G.in_degree(node))\
                        *float(G.out_degree(xnode))/m
                    else:
                        q += 0.0-float(G.in_degree(node))\
                        *float(G.out_degree(xnode))/m
        q = q/m
    else:
        for module in comm_nodes:
            nodes=comm_nodes[module]
            for node in nodes:
                for xnode in nodes:
                    if node in G.nodes() and xnode in G.nodes():
                        if G.has_edge(node,xnode):
                            q += 1.0-float(G.degree(node))\
                            *float(G.degree(xnode))/2.0/m
                        else:
                            q += 0.0-float(G.degree(node))\
                            *float(G.degree(xnode))/2.0/m
        q = q/2.0/m
    
    return q
    
