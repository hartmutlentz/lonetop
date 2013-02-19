#! /usr/bin/python
#
import sys
try: sys.path.append('/Users/lentz/Documents/GitHub_locals/my_py') # Mac
except: pass
try: sys.path.append('/home/lentz/Documents/GitHub_locals/my_py') # Beta-Cl
except: pass

import Gewindehammer_lonetop as gwh
import scipy as sc, networkx as nx, numpy as np, random

class TemporalEdgeList():
    """ Class for temporal edgelists as triples (u,v,d),
        Where (u,v) is an edge at time d.
    
    """
    def __init__(self,fname,directed,timecolumn=2):
        #list.__init__(self)
        self.edges=np.loadtxt(fname,usecols=(0,1,timecolumn),dtype='int')
        #self.edges=list(set(self.edges.flat))
        self.is_directed=directed
        self.times=set(np.loadtxt(fname,usecols=(timecolumn,),dtype='int',unpack=True))
        self.maxtime=max(self.times)
        self.mintime=min(self.times)
        self.snapshots=self.__get_snapshots()
        self.static_edges=self.__get_static_edges()
        #assert self.__has_matrix_friendly_node_labels(),"Nodenames must be 0,...,N."

    def __has_matrix_friendly_node_labels(self):
        # check if node labels are matrix friendly
        nodes1,nodes2=zip(*self.static_edges)
        nodes=[]
        nodes.extend(nodes1)
        nodes.extend(nodes2)
        nodes=set(nodes)
        nodes=list(nodes)
        nodes.sort()
        if nodes==range(len(nodes)): return True
        else: return False
    
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

    """def time_edges(self):
        # dict {t1:[(u,v),(x,y),...],...}
        te=dict([(time,[]) for time in self.times])
        
        for u,v,t in self.edges:
            te[t].append((u,v))
        return te
    """
    
    def edge_occurrence_times(self):
        # dict {(u,v):[t1,t2,...],...}
        et=dict([(se,[]) for se in self.static_edges])

        for u,v,d in self.edges:
            et[(u,v)].append(d)
        return et

    def GST(self):
        # alias
        self.shuffle_snapshot_times()
    
    def shuffle_snapshot_times(self):
        # shuffles all snapshots
        t_edges=self.snapshots
    
        new_keys=t_edges.keys()
        random.shuffle(new_keys)
    
        new_t_edges={}
        for i in t_edges:
            new_t_edges[new_keys.pop()]=t_edges[i]

        new_edges=[]
        for t in new_t_edges:
            for (u,v) in new_t_edges[t]:
                new_edges.append((u,v,t))

        self.edges=new_edges
        self.snapshots=self.__get_snapshots()

    def LST(self):
        # alias
        self.shuffle_edge_times()
    
    def shuffle_edge_times(self):
        # gives every edge new occurrence times at random. Number of occurrences is conserved.
        edge_occs=self.edge_occurrence_times()
        new_edge_occs=dict([(se,[]) for se in edge_occs])
        
        for edge in edge_occs:
            for i in range(len(edge_occs[edge])):
                new_edge_occs[edge].append(random.randint(self.mintime,self.maxtime))
        
        new_edges=[]
        for u,v in new_edge_occs:
            for t in new_edge_occs[(u,v)]:
                new_edges.append((u,v,t))
        
        self.edges=new_edges
        self.snapshots=self.__get_snapshots()

    
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
            iterations=len(self.snapshots[time])*100
        
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

            for i in range(100*len(ed)):
                first=random.sample(ed,1)[0]#random.choice(ed)
                second=random.sample(ed,1)[0]#random.choice(ed)
                if are_disjoint(first,second) and pair_not_in_G(first,second,ed):
                    return (first,second)
            return False

        def legal_graph_condition(e):
            # conditions for useful edgelists
            if len(e)<2: return False
            
            nodeset=set()
            for u,v in e:
                nodeset.add(u)
                nodeset.add(v)
                if len(nodeset)>3: return True

        # switch edges
        edges=set(self.snapshots[time][:])
        
        if legal_graph_condition(edges):
            for i in range(iterations):
                erfolg=get_legal_edgepair(edges)
                if erfolg:
                    x, y = erfolg[0], erfolg[1]
                
                    edges.remove((x[0],x[1]))
                    edges.remove((y[0],y[1]))
                    
                    edges.add((x[0],y[1]))
                    edges.add((y[0],x[1]))
                #print 'remaining: ', iterations-i

        self.snapshots[time]=list(edges)

    def __update_edges(self):
        # reads snapshots and rewrites edgelist
        new_edges=[]
        for time in self.snapshots:
            for u,v in self.snapshots[time]:
                new_edges.append((u,v,time))
        self.edges=new_edges

    def RE(self):
        # alias
        self.randomize_edges()
    
    def randomize_edges(self,maxiterations=False):
        """ Edge randomization for each graphlet
        
        """
        for i in range(self.maxtime):
            print "Randomizing ",i," of ",self.maxtime
            self.__randomize_graphlet(i,maxiterations)
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

    def random_times_uniform(self):
        """ Times at random from uniform distribution
            SLOW!
        """
        prob=self.average_size()/len(self.static_edges)
        #print prob
        for i in range(self.mintime,self.maxtime):
            print "Random times uniform. Step ",i," of ",self.maxtime
            edges=[]
            for e in self.static_edges:
                if random.random()<prob:
                    edges.append(e)
            self.snapshots[i]=edges
        self.__update_edges()

    def RT(self):
        # alias
        self.random_times()
    
    def random_times(self):
        """ Keeps the distribution of graph sizes.
            Example: 
            before: time_1: 2 edges, time_2: 5 edges, time_3: 4 edges
              edges are elements of static edges (fixed):
            after: time_1: 5 edges, time_2: 4 edges, time_3: 2 edges
              edges are chosen randomly from static graph
        """
        sizes=[len(self.snapshots[i]) for i in self.snapshots]
        random.shuffle(sizes)
        for i in range(self.maxtime+1):
            #print "Random times. Step ",i," of ",self.maxtime+1
            edges=[]
            for j in range(sizes[i]):
                edges.append(random.choice(self.static_edges))
            self.snapshots[i]=edges
        self.__update_edges()


if __name__=="__main__":
    from pprint import pprint
    the_file='out1.dat'
    E=TemporalEdgeList(the_file,True,timecolumn=2)

    E.RE()
    #print len(E.edges)
    #E=TemporalEdgeList("sociopatterns_113.dat",False)
    #pprint(E.edges)
    #E.randomize_edges()
    #E.random_times()
    #print E.average_size(),len(E.snapshots)
    #print len(E.edges)

    E.write("out1_RE.txt")

    #print E.edge_occurrence_times()
    #print E.shuffle_edge_times(E.edge_occurrence_times())


