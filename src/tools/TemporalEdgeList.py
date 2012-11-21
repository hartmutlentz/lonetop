#! /usr/bin/python
#
import sys
try: sys.path.append('/Users/lentz/Documents/GitHub_locals/my_py') # Mac
except: pass
try: sys.path.append('/home/lentz/Documents/GitHub_locals/my_py') # Beta-Cl
except: pass

import Gewindehammer as gwh
import scipy as sc, networkx as nx, numpy as np, random

class TemporalEdgeList():
    """ Class for temporal edgelists as triples (u,v,d),
        Where (u,v) is an edge at time d.
    
    """
    def __init__(self,fname,directed,timecolumn=2):
        #list.__init__(self)
        self.edges=np.loadtxt(fname,usecols=(0,1,timecolumn),dtype='int')
        self.is_directed=directed
        self.times=np.loadtxt(fname,usecols=(timecolumn,),dtype='int',unpack=True)
        self.maxtime=max(self.times)
        self.mintime=min(self.times)
        self.snapshots=self.__get_snapshots()
        self.static_edges=self.__get_static_edges()
    
    def __get_snapshots(self):
        # dict {d:[(u1,v1),(u2,v2)], ...}
        t=dict([(i,[]) for i in range(self.maxtime+1)])
        
        for u,v,d in self.edges:
            t[d].append((u,v))
        return t

    def __get_static_edges(self):
        # all edges present in the static network
        e=[(u,v) for u,v,d in self.edges]
        return list(set(e))

    def edge_occurrence_times(self):
        # dict {(u,v):[t1,t2,...],...}
        et=dict([(se,[]) for se in self.static_edges])

        for u,v,d in self.edges:
            et[(u,v)].append(d)
        return et

    def shuffle_edge_times(self,edge_occs):
        # gives every edge new occurrence times at random. Number of occurrences is conserved.
        new_edge_occs=dict([(se,[]) for se in edge_occs])
        
        for edge in edge_occs:
            for i in range(len(edge_occs[edge])):
                new_edge_occs[edge].append(random.randint(self.mintime,self.maxtime))
        
        return new_edge_occs
        
    def write(self,fname):
        """ writes self to txtfile
            
        """
        gwh.write_array(self.edges,fname)
        return

    def __randomize_graphlet(self,time,maxiterations=False):
        """
            Returns a randomized version of a graph or digraph.
            The degree sequence is conserved.
        """
        
        if maxiterations:
            iterations=maxiterations
        else:
            iterations=len(self.snapshots[time])
        
        def get_legal_edgepair(ed):
            # returns a disjoint pair of edges
            def are_disjoint(fi,se):
                # condition for disjoint edges
                if fi[0]==fi[1] or se[0]==se[1]\
                or fi[0]==se[0] or fi[0]==se[1]\
                or fi[1]==se[0] or fi[1]==se[1]:
                    return False
                else: return True

            def pair_not_in_G(x,y,eds):
                # True, if exchanged edges are not already in G. 
                if (x[0],y[1]) not in eds and (y[0],x[1]) not in eds:
                    return True
                else:
                    return False

            for i in range(10*len(ed)):
                first=random.choice(ed)
                second=random.choice(ed)
                if are_disjoint(first,second) and pair_not_in_G(first,second,ed):
                    return (first,second)
            return False
            #while the_condition(first,second)==True:
            #    first=random.choice(ed)
            #    second=random.choice(ed)
            #return (first,second)

        def legal_graph_condition(e):
            # conditions for useful edgelists
            if len(e)<2: return False
            
            nodeset=set()
            for u,v in e:
                nodeset.add(u)
                nodeset.add(v)
                if len(nodeset)>3: return True

        # switch edges
        edges=self.snapshots[time][:]
        if legal_graph_condition(edges):
            for i in range(iterations):
                erfolg=get_legal_edgepair(edges)
                if erfolg:
                    x, y = erfolg[0], erfolg[1]
                
                    edges.remove((x[0],x[1]))
                    edges.remove((y[0],y[1]))
                    edges.append((x[0],y[1]))
                    edges.append((y[0],x[1]))
                
                #print 'remaining: ', iterations-i

        self.snapshots[time]=edges

    def __update_edges(self):
        # reads snapshots and rewrites edgelist
        new_edges=[]
        for time in self.snapshots:
            for u,v in self.snapshots[time]:
                new_edges.append((u,v,time))
        self.edges=new_edges

    def randomize_edges(self):
        """ Edge randomization for each graphlet
        
        """
        for i in range(self.maxtime):
            self.__randomize_graphlet(i)
        self.__update_edges()
        self.static_edges=self.__get_static_edges()

    def number_of_nodes(self):
        # the number of nodes
        nodes=[]
        for (u,v) in self.static_edges:
            nodes.append(u)
            nodes.append(v)
        nodes=set(nodes)
        return len(nodes)

    def average_size(self):
        """ average edge density """
        all_edges=0
        for time in self.snapshots:
            all_edges += len(self.snapshots[time])
        
        return float(all_edges)/(self.maxtime-self.mintime)

    def random_times(self):
        """ Times at random from uniform distribution
        
        """
        prob=self.average_size()/len(self.static_edges)
        #print prob
        for i in range(self.mintime,self.maxtime):
            edges=[]
            for e in self.static_edges:
                if random.random()<prob:
                    edges.append(e)
            self.snapshots[i]=edges
        self.__update_edges()

    def random_times_fast(self):
        """ Keeps the distribution of graph sizes """
        sizes=[len(self.snapshots[i]) for i in self.snapshots]
        random.shuffle(sizes)
        for i in range(self.maxtime+1):
            edges=[]
            for j in range(sizes[i]):
                edges.append(random.choice(self.static_edges))
            self.snapshots[i]=edges
        self.__update_edges()


if __name__=="__main__":
    from pprint import pprint
    
    #E=TemporalEdgeList("T_edgelist.txt",True,timecolumn=3)
    E=TemporalEdgeList("sociopatterns_113.dat",False)
    #pprint(E.edges)
    #E.randomize_edges()
    print E.average_size(),len(E.snapshots)
    #E.random_times()
    E.random_times_fast()
    print E.average_size(),len(E.snapshots)
    #print E.edges

    #E.write("sociopatterns_113_edge_randomized.txt")

    #print E.edge_occurrence_times()
    #print E.shuffle_edge_times(E.edge_occurrence_times())


