# -*- encoding: utf-8 -*-

from EpiRank import epirank
from EpiRank import additional_analysis as aa

import pandas as pd
import networkx as nx

df = pd.read_csv('od_flow_data.csv', index_col=0)
print(df.head())

g = epirank.make_DiGraph(df, origin_col='origin', destination_col='destination', flow_col='flow',
                         largest_connected_component=False, exclude_selfloop=False)

epi_vals05 = epirank.run_epirank(g, daytime=0.5, d=0.95)
print('done')

"""
for k,v in epi_vals05.items():
    print(k,v)
"""
