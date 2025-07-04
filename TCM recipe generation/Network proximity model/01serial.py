#!/usr/bin/env python

import pandas as pd
import networkx as nx
import numpy as np
import re, os, itertools, random

# 人类蛋白相互作用网络
def ppi_gen(ppi_file, is_adjacency_list=False, lcc=True, sep=','):
    '''
    ppi_file: 用于生成人类蛋白相互作用的文件，可为相互作用的两列蛋白或邻近表文件
    is_adjacency_list: 是否为邻近表文件
    lcc: 是否取网络的最大连通分量
    '''
    if is_adjacency_list:
        return None
    else:
        df = pd.read_csv(ppi_file,
                         usecols=['EntrezA', 'EntrezB'],
                         dtype={'EntrezA': int,
                                'EntrezB': int},
                         sep=sep)
        G = nx.from_pandas_edgelist(df, 'EntrezA', 'EntrezB')

    if lcc:
        G = G.subgraph(max(nx.connected_components(G), key=len))

    return G

# 网络节点度的分布，生成Zscore时会用
def net_degree_distribution(net_data, min_bin_size=100):
    '''
    net_data: 人类蛋白相互作用网络
    min_bin_size: 将人类蛋白相互作用网络中节点按度进行分组，每组中最少要包含的节点数目
    '''
    node2degrees = dict(net_data.degree())
    degree2nodes = {}
    for node, degree in node2degrees.items():
        degree2nodes.setdefault(degree, []).append(node)

    degrees = sorted(degree2nodes.keys())

    bins = []
    i = 0
    while i < len(degrees):
        degree = degrees[i]
        low_degree = degree
        nodes = degree2nodes[degree]
        while len(nodes) < min_bin_size:
            i = i + 1
            if i == len(degrees):
                break
            nodes.extend(degree2nodes[degrees[i]])
        if i == len(degrees):
            i = i - 1
        high_degree = degrees[i]
        
        i = i + 1
        if len(nodes) < min_bin_size:
            low_degree_, high_degree_, nodes_ = bins[-1]
            bins[-1] = (low_degree_, high_degree, nodes_ + nodes)
        else:
            bins.append((low_degree, high_degree, nodes))
            
    return bins

# 处理药物－蛋白、疾病－基因输入文件，去除名称中特殊字符以便于生成输出文件
def process_assoc_file(in_file, file_type, association, sep=','):
    '''
    in_file: 输入文件，保存有药物－蛋白或疾病－基因之间的对应关系
    file_type: 药物－蛋白文件("chemical")或疾病－基因文件("disease")
    association: 列名称，用于表明药物/疾病、EntrezID存于哪两列，
                 如{"chemical": "entrez_id"}
    '''
    (col_cd, col_gen) = [(k, v) for k, v in association.items()][0]
    df = pd.read_csv(in_file, usecols=[col_cd, col_gen], sep=sep,
                     dtype={col_gen: int}).dropna()
    df = df.rename(columns={col_cd: file_type, col_gen: 'entrez_id'})
    df[file_type + '_orig'] = df[file_type]
    df[file_type] = df[file_type].apply(
        lambda x: re.sub('[^A-Za-z0-9]+', '',x))
    
    return df.drop_duplicates()

# 得到所有药物－疾病组合
def chemical_disease_product(chemical_data, disease_data):
    '''
    chemical_data: 药物－蛋白数据
    disease_data: 疾病－基因数据
    '''
    return set(itertools.product(
        chemical_data['chemical'],
        disease_data['disease']))

# 得到还有多少药物－疾病组合需要计算
def chemical_disease_remained(chemical_disease_pairs, tmp_dir='tmp'):
    '''
    chemical_disease_pairs: 药物－疾病对
    tmp_dir: 已完成药物－疾病对计算的临时文件存放路径
    '''
    finished = set()
    if os.path.isdir(tmp_dir):
        for f in os.listdir(tmp_dir):
            if f.endswith('.csv'): # chemical_disease.csv
                finished.add(
                    tuple(f.replace('.csv', '').split('_')))
    else:
        os.mkdir(tmp_dir)

    chemical_disease_pairs = chemical_disease_pairs - finished
    print(f'There are {len(chemical_disease_pairs)} ' +
          f'chemical_disease pairs to be calcuated')
    return chemical_disease_pairs

# 一对药物－疾病对应的蛋白、基因
def get_chemical_disease_nodes(chemical, disease,
                               chemical_data,
                               disease_data):
    '''
    chemical: 药物－疾病对中药物名称
    disease: 药物－疾病对中疾病名称
    chemical_data: 所有药物－蛋白数据
    disease_data: 所有疾病－基因数据
    '''
    df_chemical = chemical_data[
        chemical_data['chemical'] == chemical ]
    df_disease = disease_data[
        disease_data['disease'] == disease ]
    proteins = set(df_chemical['entrez_id'])
    genes = set(df_disease['entrez_id'])
    return proteins, genes

# 计算邻近度
def cal_proximity(network, nodes_from, nodes_to):
    '''
    network: 人类蛋白相互作用网络
    nodes_from: 一组节点，从该组节点开始计算路径长度
    nodes_to: 一组节点，计算路径长度时到该组节点结束
    '''
    min_lengths = []
    mean_lengths = []
    for node_from in nodes_from:
        lengths = []
        for node_to in nodes_to:
            if (nx.has_path(network, node_from, node_to)):
                lengths.append(
                    nx.shortest_path_length(network,
                                            node_from,
                                            node_to))
        min_lengths.append(np.min(lengths))
        mean_lengths.append(np.mean(lengths))

    ds = np.mean(mean_lengths)
    dc = np.mean(min_lengths)
    return ds, dc

# 取出度相似的节点
def get_equivalent_nodes(network, orig_nodes, degree_bins):
    '''
    network: 人类蛋白相互作用网络
    orig_nodes: 原始节点，取出的节点的度与该组节点的度相似
    degree_bins: 网络中度的分布
    '''
    nodes_found = {}
    for node in orig_nodes:
        degree = network.degree(node)
        nodes_gen = (nodes for l, h, nodes in degree_bins
                     if ((l <= degree) and (h >= degree)))
        nodes = list(next(nodes_gen))
        nodes.remove(node)
        nodes_found[node] = nodes
        
    return nodes_found
    

# 从度相似的节点中随机选取节点
def get_random_nodes(network, orig_nodes, degree_bins,
                     n_random, random_seed):
    '''
    network: 人类蛋白相互作用网络
    orig_nodes: 参考节点组，随机取出的节点与该组内节点度相似
    degree_bins: 网络中度的分布
    n_random: 随机选取的组数
    random_seed: 随机数种子
    '''
    random.seed(random_seed)

    nodes = []
    for i in range(n_random):
        random_nodes = set()
        nodes2eq_nodes = get_equivalent_nodes(network,
                                              orig_nodes,
                                              degree_bins)
        for _, equivalents in nodes2eq_nodes.items():
            chosen = random.choice(equivalents)
            for k in range(20):
                if chosen in random_nodes:
                    chosen = random.choice(equivalents)
            random_nodes.add(chosen)
        random_nodes = list(random_nodes)
        nodes.append(random_nodes)
        
    return nodes

# 邻近度打分
def cal_zscore(network, orig_nodes_from, orig_nodes_to,
               orig_ds, orig_dc, ref_size=1000):
    '''
    network: 人类蛋白相互作用网络
    orig_nodes_from: 原始节点，计算邻近度时从这些节点开始，计算打分时以这些节点为参考
                     随机选取度与之相似的一组节点
    orig_nodes_to: 原始节点，计算邻近度时到这些节点结束，计算打分时以这些节点为参考
                   随机选取度与之相似的一组节点
    orig_ds: 药物－疾病（平均）邻近度
    orig_dc: 药物－疾病（最小）邻近度
    ref_size: 随机选取的节点组的个数，用于构建参考邻近度分布
    '''
    min_bin_size = 2*max(len(orig_nodes_from), len(orig_nodes_to))
    degree_bins = net_degree_distribution(network, min_bin_size)
    
    random_seed = 452456
    random_nodes_from = get_random_nodes(network,
                                         orig_nodes_from,
                                         degree_bins,
                                         n_random=ref_size,
                                         random_seed=random_seed)
    random_nodes_to = get_random_nodes(network,
                                       orig_nodes_to,
                                       degree_bins,
                                       n_random=ref_size,
                                       random_seed=random_seed)
    null_min_lengths = []
    null_mean_lengths = []
    for i in range(ref_size):
        nodes_from = random_nodes_from[i]
        nodes_to = random_nodes_to[i]
        null_mean_length, null_min_length = cal_proximity(network,
                                                          nodes_from,
                                                          nodes_to)
        null_min_lengths.append(null_min_length)
        null_mean_lengths.append(null_mean_length)

    mu_dc, sigma_dc = np.mean(null_min_lengths), np.std(null_min_lengths)
    mu_ds, sigma_ds = np.mean(null_mean_lengths), np.std(null_mean_lengths)

    zs = (orig_ds - mu_ds)/sigma_ds
    zc = (orig_dc - mu_dc)/sigma_dc
    ps = np.mean(np.array(null_mean_lengths) < orig_ds)
    pc = np.mean(np.array(null_mean_lengths) < orig_dc)
    return zs, zc, ps, pc

# 计算Zscore
def get_zscore(chemical_disease_pair, chemical_data, disease_data,
               net_data, tmp_dir='tmp', distance_matrix=None):
    '''
    chemical_disease_pair: 一对药物－疾病
    chemical_data: 药物－蛋白数据
    disease_data: 疾病－基因数据
    net_data: 人类相互作用网络
    tmp_dir: 计算完成的数据存放路径
    distance_matrix: 人类相互作用网络距离矩阵
    '''
    chemical, disease = chemical_disease_pair
    proteins, genes = get_chemical_disease_nodes(
        chemical, disease, chemical_data, disease_data)
    proteins = proteins & set(net_data)
    genes = genes & set(net_data)

    ds, dc = cal_proximity(network=net_data,
                           nodes_from=proteins,
                           nodes_to=genes)
    
    zs, zc, ps, pc = cal_zscore(network=net_data,
                                orig_nodes_from=proteins,
                                orig_nodes_to=genes,
                                orig_ds=ds, orig_dc=dc)

    chemical_name = chemical_data[chemical_data['chemical']
                                  == chemical]['chemical_orig']
    chemical_name = chemical_name.values[0]
    disease_name = disease_data[disease_data['disease']
                                == disease]['disease_orig']
    disease_name = disease_name.values[0]

    n_mapped_chemical = len(proteins)
    n_mapped_disease = len(genes)    
    df = pd.DataFrame([[chemical_name, n_mapped_chemical,
                        disease_name, n_mapped_disease,
                        ds, zs, ps,
                        dc, zc, pc]],
                      columns=['chemical', 'n_mapped_chemical',
                               'disease', 'n_mapped_disease',
                               'shortest', 'z_shortest', 'p_shortest',
                               'closest', 'z_closest', 'p_closest'])
    out_file = '_'.join(chemical_disease_pair) + '.csv'
    out_file = tmp_dir + '/' + out_file
    df.to_csv(out_file, index=None)
    # print(df.iloc[0])
    
        
# 运行
def run(net_data, chemical_data, disease_data,
        tmp_dir='tmp', n_procs=None, distance_matrix=None):
    '''
    net_data: 人类蛋白相互作用组数据
    chemical_data: 药物－蛋白数据
    disease_data: 疾病－基因数据
    tmp_dir: 存储已计算完成了的药物－疾病邻近度结果
    n_procs: 并行核数，若为None则不进行并行计算
    distance_matrix: 人类蛋白相互作用网络距离矩阵
    '''
    # chemical_data = chemical_data[chemical_data['entrez_id']
    #                               .isin(net_data.nodes())]
    # disease_data = disease_data[disease_data['entrez_id']
    #                             .isin(net_data.nodes())]
    # 生成药物－疾病对
    chemical2disease = chemical_disease_product(
        chemical_data = chemical_data,
        disease_data = disease_data)

    # 去除已计算完成的药物－疾病对
    chemical2disease = chemical_disease_remained(
        chemical2disease, tmp_dir=tmp_dir)

    if not n_procs:
        for chemical_disease_pair in chemical2disease:
            get_zscore(chemical_disease_pair,
                       chemical_data, disease_data,
                       net_data, tmp_dir=tmp_dir,
                       distance_matrix=distance_matrix)

if __name__ == '__main__':
    ppi_file = 'HumanInteractome_v2017.csv'
    disease_file = 'PeroxisomalDisorders.csv'
    chem_file = 'Genistein.csv'

    # 生成人类蛋白相互作用组
    ppi_net = ppi_gen(ppi_file)

    # 处理药物、疾病输入文件
    df_disease = process_assoc_file(disease_file,
                                    file_type='disease',
                                    association={'disease': 'entrez_id'})
    df_chemical = process_assoc_file(chem_file,
                                     file_type='chemical',
                                     association={'chemical': 'entrez_id'})
    run(ppi_net, df_chemical, df_disease)
    
    
