# 
import pandas as pd
import networkx as nx
import random

# 读入化合物名称与靶点蛋白
df_chem = pd.read_csv('../data/PolyphenolProteinInteractions.csv', 
                      dtype={'entrez_id': int})
#   chemical  pubchem_compound_id  entrez_id  symbol
# 0   butein            5281222.0       6718  AKR1D1


# 读入疾病分类与致病
df_disease = pd.read_csv('../data/GenesDisease.csv', dtype={'entrez_id': int})
#            disease  entrez_id
# 0  kidney diseases      79663

# 读入人类相互作用组，并构成网络
df_ppi = pd.read_csv('../data/HumanInteractome_v2017.csv',
                     dtype={'EntrezA': int, 'EntrezB': int})
#    EntrezA  EntrezB
# 0        1      310

G = nx.from_pandas_edgelist(df_ppi, 'EntrezA', 'EntrezB')

# 只考虑相互作用组中最大连通分量(LCC)
largest_cc = max(nx.connected_components(G), key=len)
g = G.subgraph(largest_cc)


# 化合物靶点蛋白与最大连通分量共有的节点
# '(-)-epigallocatechin 3-o-gallate'  == EGCG
chemical = '(-)-epigallocatechin 3-o-gallate'
target_nodes = df_chem[df_chem['chemical'] == chemical]['entrez_id']
nodes_from = set(target_nodes) & set(g)  # set(g) == set(g.nodes())
n_chem = len(nodes_from)

# 疾病致病基因与最大连通分量共有的节点
# 'body weight' 14 个节点， d_s=2.609, d_c=2.020, z_s = 0.278, z_c = -1.527
disease ='nervous system diseases'
disease_nodes = df_disease[df_disease['disease'] == disease]['entrez_id']
nodes_to = set(disease_nodes) & set(g)
n_disease = len(nodes_to)

def cal_lengths(g, nodes_from, nodes_to):
    min_lengths = []
    mean_lengths = []
    for node_from in nodes_from:
        f_length = []
        for node_to in nodes_to:
            if (nx.has_path(g, node_from, node_to)):
                f_length.append(
                    nx.shortest_path_length(g, node_from, node_to))

        f_length = pd.Series(f_length)
        min_lengths.append(f_length.min())
        mean_lengths.append(f_length.mean())

    min_lengths = pd.Series(min_lengths)
    mean_lengths = pd.Series(mean_lengths)
    
    return min_lengths, mean_lengths

min_lengths, mean_lengths = cal_lengths(g, nodes_from, nodes_to)

###
# z_scores 计算相关
n_random = 1000
min_bin_size = 2*max(len(nodes_from), len(nodes_to))
# min_bin_size = 100

degrees = dict(g.degree())
# {1: 25, 310: 135, 368: 7, 1026: 330, 2232: 28, 2886: 58, ...}
# 节点：度

degree_to_nodes = {}
for node, degree in degrees.items():
    degree_to_nodes.setdefault(degree, []).append(node)
# {25: [1, 412,...], 135: [310, ...], ...}
# 度：[节点]

values = list(degree_to_nodes.keys())   # 所有的degree
# [25, 135, 7, ... ]
values.sort()
# [1, 2, 3, 4, ... ]

bin_size = min_bin_size
bins = []   # 若一个度包含的节点较少，则与邻近度合并
i = 0
while i < len(values):
    low = values[i]
    val = degree_to_nodes[values[i]]     # 节点
    while len(val) < bin_size:
        i += 1
        if i == len(values):
            break
        val.extend(degree_to_nodes[values[i]])
    if i == len(values):
        i -= 1
    high = values[i]
    i += 1
    if len(val) < bin_size:
        low_, high_, val_ = bins[-1]
        bins[-1] = (low_, high, val_ + val)
    else:
        bins.append((low, high, val))

# [(1, 2, [10321, 10944, ...]), (3, 5, [79006, ...]), ...]
# [(度，度，节点), ...]

def get_degree_equivalents(degree_holders, bins, g):

    """
    选出和degree_holders中各节点具有相同度的节点
    参数:
        degree_holders: 要匹配度的节点
        bins: 将具有相同（相似）度的节点放到了一起
        g: 图
    """

    nodes_found = {}
    for node in degree_holders:
        d = g.degree(node)   # 获得某一节点的度
        for l, h, nodes in bins:
            if l <= d and h >= d:   # d 位于某一区间内，取出相应的节点
                mod_nodes = list(nodes)
                mod_nodes.remove(node)   # 节点本身和自己具有相同的度，移除
                nodes_found[node] = mod_nodes
                break     #  度位于的区间是唯一的，找到一个就不用往下找了

    return nodes_found

# {1666: [1200, 3553, ...], 84869: [9923, ...], ...}
# {参考节点: 具有相同度的节点, ...}
# nodes_found = get_degree_equivalents(nodes_from, bins, g)

def get_random_nodes(degree_holders, g, bins, n_random,
                     min_bin_size, random_seed):
    random.seed(random_seed)    # 设置随机选取的种子

    values = []
    for i in range(n_random):
        nodes_random = set()
        # 取出具有相同度的节点
        node_to_equivalent_nodes = \
            get_degree_equivalents(degree_holders, bins, g)
        for node, equivalent_nodes in node_to_equivalent_nodes.items():
            chosen = random.choice(equivalent_nodes)  # 从具有相同度的节点中随机选
            for k in range(20):  # 尽量避免选到同一个节点
                if chosen in nodes_random:
                    chosen = random.choice(equivalent_nodes)
            nodes_random.add(chosen)
        nodes_random = list(nodes_random)
        values.append(nodes_random)
    return values

# [[2304, 84230,...], [653440, 797, ...], ...]
# [样本1, 样本2, ...] 共n_random个样本
# random_seed = 452456
# random_nodes_from = get_random_nodes(nodes_from, g, bins, n_random,
#                                      min_bin_size, random_seed)

# seed = 452456 为了验证 
random_seed = 452456
random_nodes_from = get_random_nodes(nodes_from, g=g, bins=bins,
                                     n_random=n_random,
                                     min_bin_size=min_bin_size,
                                     random_seed=random_seed)

random_nodes_to = get_random_nodes(nodes_to, g=g, bins=bins,
                                   n_random=n_random,
                                   min_bin_size=min_bin_size,
                                   random_seed=random_seed)
random_values_list = list(zip(random_nodes_from, random_nodes_to))
# 两列， n_random行

null_min_lengths = []
null_mean_lengths = []
for i, values_random in enumerate(random_values_list):
    nodes_from, nodes_to = values_random
    null_min, null_mean = cal_lengths(g, nodes_from, nodes_to)
    null_min_lengths.append(null_min.mean())
    null_mean_lengths.append(null_mean.mean())

import numpy as np
mu_dc, sigma_dc = np.mean(null_min_lengths), np.std(null_min_lengths)
mu_ds, sigma_ds = np.mean(null_mean_lengths), np.std(null_mean_lengths)

# Python 参数引用传递
# min_lengths, mean_lengths = cal_lengths(g, nodes_from, nodes_to)
dc, ds = min_lengths.mean(), mean_lengths.mean()

z_ds = (ds - mu_ds)/sigma_ds 
z_dc = (dc - mu_dc)/sigma_dc 
p_ds = np.mean(np.array(null_mean_lengths) < ds)
p_dc = np.mean(np.array(null_min_lengths) < dc)

#print(('n_disease: {}, n_chem: {}, ds: {}, z_ds: {}, ' + 
#       'p_ds: {}, dc: {}, z_dc: {}, p_dc: {}')
#      .format(n_disease, n_chem, ds, z_ds, p_ds, dc, z_dc, p_dc))

print(('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}')
      .format(disease, n_disease, n_chem, chemical, ds, z_ds, p_ds, 
             dc, z_dc, p_dc))

