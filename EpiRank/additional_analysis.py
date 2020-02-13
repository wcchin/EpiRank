# -*- encoding: utf-8 -*-

import pandas as pd            # main data format
import networkx as nx          # main data format
import numpy as np             # log10
from scipy import stats as st  # pearsonr and spearmanr
import copy                    # copy and deepcopy

def htbreak(adic, g=4):
    alist = [ v for k,v in adic.items() ]
    temp = copy.copy(alist)
    breaks = []
    for i in range(g-1):
        avg = sum(temp) / float(len(temp))
        breaks.append(avg)
        temp2 = [ v for v in temp if v>avg ]
        temp = temp2
    #alist2 = []
    adic2 = {}
    for k,v in adic.items():
        lvl = None
        for i in range(len(breaks)):
            if v<=breaks[i]:
                lvl = i
                break
        if lvl is None:
            lvl = len(breaks)
        #alist2.append(lvl)
        adic2[k] = lvl

    #print(alist2)
    return adic2, breaks


def calculate_IOratio(g, exclude_selfloop=True):
    g2 = copy.deepcopy(g)
    if exclude_selfloop:
        g2.remove_edges_from(g2.selfloop_edges())
    indeg = g2.in_degree(weight='weight')
    oudeg = g2.out_degree(weight='weight')
    ioRatio = {}
    for n in g2.nodes():
        if oudeg[n]>0:
            ioRatio[n] = np.log10(indeg[n]/oudeg[n])
        else:
            ioRatio[n] = 0
    return ioRatio


def calculate_metrices(g, return_dic=True, d=0.95, number_of_loops=1000, weighted=True):

    if weighted:
        pagerank = nx.pagerank(g, alpha=d, max_iter=number_of_loops)
        hub_rank, authority_rank = nx.hits(g, max_iter=number_of_loops)
    else:
        g2 = copy.deepcopy(g)
        for n1,n2,wd in g2.edges(data=True):
            g2[n1][n2]['weight'] = float(1.0)
        pagerank = nx.pagerank(g2, alpha=d, max_iter=number_of_loops)
        hub_rank, authority_rank = nx.hits(g2, max_iter=number_of_loops)

    metrices = [pagerank, hub_rank, authority_rank]
    metrices_names = ['pagerank', 'hub_rank', 'authority_rank']
    """
    cal_res = { n:{} for n in g.nodes() }
    for dic, name in zip(metrices, metrices_names):
        for n,v in dic.items():
            cal_res[n].update({name: v})
    """
    #cal_res = { name:res for name,res in zip(metrices_names, metrices) }
    cal_res = list( zip(metrices_names, metrices) )
    if return_dic:
        return cal_res
    else:
        cal_res2 = { a:b for a,b in cal_res }
        df_res = pd.DataFrame.from_dict(cal_res2) #, orient='index')
        df_res = df_res[metrices_names]
        return df_res

"""
def calculate_metrices(g, return_dic=True, d=0.95, number_of_loops=1000, weighted=True):

    windeg = g.in_degree(weight='weight') ## weighted
    woudeg = g.out_degree(weight='weight') ## weighted
    wtodeg = { k:windeg[k]+woudeg[k] for k in windeg.keys() } ## weighted

    if weighted:
        pagerank = nx.pagerank(g, alpha=d, max_iter=number_of_loops)
        hub_rank, authority_rank = nx.hits(g, max_iter=number_of_loops)
    else:
        g2 = copy.deepcopy(g)
        for n1,n2,wd in g2.edges(data=True):
            g2[n1][n2]['weight'] = float(1.0)
        pagerank = nx.pagerank(g2, alpha=d, max_iter=number_of_loops)
        hub_rank, authority_rank = nx.hits(g2, max_iter=number_of_loops)
    #eccentricity = nx.eccentricity(g)
    degree_centr = nx.degree_centrality(g) ## normalized
    indegree_centr = nx.in_degree_centrality(g) ## normalized
    outdegree_centr = nx.out_degree_centrality(g) ## normalized
    #closeness_centr = nx.closeness_centrality(g) ## normalized
    #betweenness_centr = nx.betweenness_centrality(g, weight='weight') ## normalized
    eigenvector_centr = nx.eigenvector_centrality(g, max_iter=number_of_loops) ## normalized

    metrices = [pagerank, hub_rank, authority_rank,
                windeg, woudeg, wtodeg,
                degree_centr, indegree_centr, outdegree_centr,
                eigenvector_centr]
    #, closeness_centr, betweenness_centr, eccentricity
    metrices_names = ['pagerank', 'hub_rank', 'authority_rank',
                      'weighted_indegree', 'weighted_outdegree', 'weighted_total_degree',
                      'degree_centrality', 'indegree_centrality', 'outdegree_centrality',
                      'eigenvector_centrality']
    #, 'closeness_centrality', 'betweenness_centrality', 'eccentricity'
    #cal_res = { name:res for name,res in zip(metrices_names, metrices) }
    cal_res = list( zip(metrices_names, metrices) )
    if return_dic:
        return cal_res
    else:
        df_res = pd.DataFrame.from_dict(cal_res) #, orient='index')
        return df_res
"""

def get_pearson_cor(dic1, dic2, simplify=True):
    keys = dic1.keys()
    numbers1 = [ dic1[k] for k in keys ]
    numbers2 = [ dic2[k] for k in keys ]
    v1, v2 = st.pearsonr( numbers1, numbers2 ) ### cor coef, p-val
    if simplify:
        star = ''
        if v2<=0.001:
            star = '***'
        elif v2<=0.01:
            star = '**'
        elif v2<=0.05:
            star = '*'
        else:
            star = ''
        #v3 = str(v1)+star
        return round(v1, 4), star
    else:
        return v1, v2


def get_spearman_cor(dic1, dic2, simplify=True):
    keys = dic1.keys()
    numbers1 = [ dic1[k] for k in keys ]
    numbers2 = [ dic2[k] for k in keys ]
    v1, v2 = st.spearmanr( numbers1, numbers2 ) ### cor coef, p-val
    if simplify:
        star = ''
        if v2<=0.001:
            star = '***'
        elif v2<=0.01:
            star = '**'
        elif v2<=0.05:
            star = '*'
        else:
            star = ''
        #v3 = str(v1)+star
        return round(v1, 4), star
    else:
        return v1, v2


def validate(validate_d, epi, trad):
    validate_pearson = {}
    validate_spearman = {}
    validate_pearson['epirank'] = get_pearson_cor(epi, validate_d)
    validate_spearman['epirank'] = get_spearman_cor(epi, validate_d)
    orders = ['epirank']
    for n,m in trad:
        validate_pearson[n] = get_pearson_cor(m, validate_d)
        validate_spearman[n] = get_spearman_cor(m, validate_d)
        orders.append(n)
    return validate_pearson, validate_spearman, orders


def prepare_all_table(listtuple_of_dicts, g=None):
    all_dic = {}
    col_order = []
    for n,d in listtuple_of_dicts:
        all_dic[n] = d
        col_order.append(n)
    df_res = pd.DataFrame.from_dict(all_dic)
    df_res = df_res[col_order]
    return df_res
