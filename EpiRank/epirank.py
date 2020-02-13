# -*- encoding: utf-8 -*-

import pandas as pd    # input DataFrame
import networkx as nx  # main data format
import numpy as np     # making matrix and multiplication

def make_DiGraph(df, origin_col='origin', destination_col='destination', flow_col='flow', largest_connected_component=True, exclude_selfloop=True):
    origins = df[origin_col].tolist()
    destinations = df[destination_col].tolist()
    flows = df[flow_col].tolist()
    g = nx.DiGraph()
    for o,d,f in zip(origins, destinations, flows):
        g.add_edge(o,d, weight=float(f))
    if largest_connected_component:
        g = max(nx.weakly_connected_component_subgraphs(g), key=len)
    if exclude_selfloop:
        g.remove_edges_from(g.selfloop_edges())
    print('graph construction done,'+' no. nodes: '+str( g.number_of_nodes()) +', no. edges: '+str(g.number_of_edges()))
    return g

def run_epirank(g, d=0.95, daytime=0.5, number_of_loops=1000, exfac=None):
    ## exfac should be a matrix with the shape (N, 1)
    Ncount = g.number_of_nodes()
    ### prepare matrix ###
    print('start preparing matrices')
    CN_C = nx.to_numpy_matrix(g, weight = 'weight')
    CN_T = CN_C
    csum = CN_T.sum(axis=0)
    CN = np.zeros((Ncount,Ncount))
    for i in range(CN_T.shape[0]):
        for j in range(CN_T.shape[1]):
            s = float(csum[0,j])
            if s>0:
                v = float(CN_T[i,j])
                CN[i,j] = v/s
    ODt = CN_T.transpose()
    osum = ODt.sum(axis=0)
    CNt = np.zeros((Ncount,Ncount))
    for i in range(ODt.shape[0]):
        for j in range(ODt.shape[1]):
            s = float(osum[0,j])
            if s>0:
                v = float(ODt[i,j])
                CNt[i,j] = v/s

    ### prepare external factor matrix if not specified ###
    exfac_matrix = get_exfac(exfac, g)

    ### initialize epidemic risk value matrix ###
    epidemic_risk = np.matrix(np.ones((1, Ncount))/float(Ncount)).transpose() ## init EpiRank values
    print('preparation done, start iterating')
    ### start running ###
    for i in range(number_of_loops):
        old_epidemic_risk = epidemic_risk.copy()
        epidemic_risk = (1. - d) * exfac_matrix + d * (daytime * CNt * epidemic_risk + (1. - daytime) * CN * epidemic_risk)
        if np.ma.allequal(epidemic_risk, old_epidemic_risk): break
    #print('iteration count:', i)
    print('epirank calculation done after iteration: '+str(i))

    ### prepare result dic ###
    vals = [ i[0] for i in epidemic_risk.tolist() ]
    nodes = list(g.nodes())
    epi_value = {}
    for i in range(len(epidemic_risk)):
        epi_value[nodes[i]] = vals[i]
    return epi_value


def get_exfac(exfackey, g):
    Ncount = g.number_of_nodes()
    if exfackey is None:
        exfac = np.matrix(np.ones((1, Ncount))/float(Ncount)).transpose()
    elif isinstance(exfackey, dict):
        ## please make sure all node is in dict ##
        exlist = []
        for n in g.nodes():
            if n in exfackey:
                v = exfackey[n]
            else:
                print('the dictionary do not contain info for node %s, set as 0'%(n))
                v = 0
            exlist.append(v)
        exsum = float(sum(exlist))
        exlist2 = [ float(ex)/exsum for ex in exlist ]
        exfac = np.matrix(exlist2).transpose()
        #print exfac.shape
    elif isinstance(exfackey, str):
        ## please make sure all nodes have this attribute name ##
        exlist = []
        for n,v in g.nodes(data=True):
            if exfackey in v:
                v2 = v[exfackey]
            else:
                print('the node %s do not have the attribute %s, set as 0'%(n, exfackey))
                v2 = 0
            exlist.append(v2)
        exsum = float(sum(exlist))
        exlist2 = [ float(ex)/exsum for ex in exlist ]
        exfac = np.matrix(exlist2).transpose()
    elif isinstance(exfackey, list):
        ## please make sure the list length match the number of nodes ##
        exlist = exfackey[:Ncount]
        exsum = float(sum(exlist))
        exlist2 = [ float(ex)/exsum for ex in exlist ]
        exfac = np.matrix(exlist2).transpose()
    elif isinstance(exfackey, np.matrix):
        ## please make sure the length of matrix has the same number as nodes ##
        if exfackey.shape==(Ncount, 1):
            exfac = exfackey
        elif exfackey.shape==(1, Ncount):
            exfac = exfackey.transpose()
        else:
            print('the shape is neither (%s, 1) nor (1, %s), force change to default exfac'%(str(Ncount), str(Ncount)))
            exfac = get_exfac(None, g)
    elif isinstance(exfackey, np.ndarray):
        exfackey = np.matrix(exfackey)
        if exfackey.shape==(Ncount, 1):
            exfac = exfackey
        elif exfackey.shape==(1, Ncount):
            exfac = exfackey.transpose()
        else:
            print('the shape is neither (%s, 1) nor (1, %s), force change to default exfac'%(str(Ncount), str(Ncount)))
            exfac = get_exfac(None, g)
    else:
        print('the input is not recognized, force back to default exfac')
        exfac = get_exfac(None, g)
    return exfac
