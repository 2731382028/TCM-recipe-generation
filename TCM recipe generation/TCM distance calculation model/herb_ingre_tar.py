
import pandas as pd
from proximity_key import *
import numpy as np
from collections import defaultdict
from functools import reduce
import pymysql
import pymysql.cursors



def qc_herb(data):
    """
    清洗草药-成分数据，返回有效的草药-成分字典。
    """
    herb_ingre_dict = data.groupby('herb_id')['ingredients_id'].apply(list).to_dict()
    herb_ingre_dict = {key: [m for m in value if isinstance(m, str)] for key, value in herb_ingre_dict.items() if
                       len(value) != 0}
    herb_ingre_dict_2 = {key: [m for m in value if m[1:].isdigit() == True] for
                         key, value in herb_ingre_dict.items()}
    herb_ingre_dict3 = {key: value for key, value in herb_ingre_dict_2.items() if len(value) != 0}
    print(f"有效草药数量: {len(herb_ingre_dict)}")
    print(f"herb_ingre_dict 示例: {list(herb_ingre_dict.items())[:5]}")  # 打印前5个条目
    return herb_ingre_dict3


class Ingredients:
    def __init__(self, filename):
        """初始化 Ingredients 类，加载数据并初始化相关属性。"""
        self.filename = filename
        self.data = self.read_file()
        self.ingredients = self.data.index.tolist()  # 成分列表
        self.ingre_tar_dict = None  # 成分-目标字典
        self.ingre_id_name_dict_value = self.ingre_id_name_dict()

    def ingre_id_name_dict(self): 
        """成分ID与名称的映射字典（修复版）"""
        try:
            return {
                idx: self.data.loc[idx, 'Ingredient Name'] 
                for idx in self.ingredients
                if 'Ingredient Name' in self.data.columns
            }
        except Exception as e:
            print(f"映射字典创建失败: {str(e)}")
            return {}

    def read_file(self):
        """
        读取数据文件并返回 DataFrame。
        """
        data = pd.read_csv(self.filename, sep=',', encoding='utf-8-sig')
        if 'ingredients_id' not in data.columns:
            raise ValueError("Missing ingredients_id column")
        data = data.dropna(subset=['ingredients_id'])
        data['ingredients_id'] = 'I' + data['ingredients_id'].astype(str)
        # 处理重复索引
        if data.index.duplicated().any():
            print("发现重复成分ID，保留第一个出现项")
            data = data[~data.index.duplicated(keep='first')]
        return data.set_index('ingredients_id')
    
    def ingredients_target_dict(self, G_nodes):
        """优化后的目标字典初始化方法"""
        print(f"正在初始化目标字典，网络节点数: {len(G_nodes)}")
        self.ingre_tar_dict = defaultdict(list)
        g_set = set(G_nodes)
        
        total = len(self.data)
        processed = 0
        
        for idx, row in self.data.iterrows():
            try:
                raw_targets = row.get('targets', None)
                if pd.isna(raw_targets):
                    continue
                
                processed_targets = set()
                for t in str(raw_targets).split(','):
                    clean_t = t.strip().lstrip('T')
                    if clean_t:
                        processed_targets.add(clean_t)
                
                valid_targets = processed_targets & g_set
                if valid_targets:
                    self.ingre_tar_dict[idx] = list(valid_targets)
                
                processed += 1
                if processed % 1000 == 0:
                    print(f"处理进度: {processed}/{total} ({processed/total:.1%})")
                    
            except Exception as e:
                print(f"处理成分 {idx} 时出错: {str(e)}")
                continue
        
        print(f"目标字典初始化完成，有效成分数: {len(self.ingre_tar_dict)}")
        return self.ingre_tar_dict

    
    def ingre_ingre_dis(self, ingre_from, ingre_to, network, distance_method):
        """
        计算两个成分之间的距离。
        """
        if any(ingre not in self.ingre_tar_dict.keys() for ingre in [ingre_from, ingre_to]):
            return None
        else:
            nodes_from = self.ingre_tar_dict[ingre_from]
            nodes_to = self.ingre_tar_dict[ingre_to]
            length_dict = Sets_Lengths(nodes_from, nodes_to).target_lengths(network)
            dis_obj = Network_Distance(nodes_from, nodes_to, length_dict)
            distance = dis_obj.network_distance(distance_method)
            return distance

    def ingre_ingre_dis_all(self, ingre_from, ingre_to, network):
        """
        计算两个成分之间的所有距离方法。
        """
        distance_method_list = ['separation', 'closest', 'shortest', 'kernel', 'center']
        return {method: self.ingre_ingre_dis(ingre_from, ingre_to, network, method)
                for method in distance_method_list}


class Herb:
    def __init__(self, filename):
        """
        初始化 Herb 类，加载数据并初始化相关属性。
        """
        self.filename = filename
        self.data_all = self.read_data2()
        print(f"原始草药数量: {len(self.data_all['herb_id'].unique())}")
        self.data = self.precess_data()
        self.herb_ingre_dict_all = qc_herb(self.data)
        print(f"清洗后有效草药数量: {len(self.herb_ingre_dict_all)}")
        self.herbs = list(self.herb_ingre_dict_all.keys())
        all_ingredients = [ingre for sublist in self.herb_ingre_dict_all.values() for ingre in sublist]
        self.ingredients = list(set(all_ingredients))  # 去重
        self.herb_ingre_dict = self._init_herb_ingre_dict()  # 仅初始化一次
  

    def _init_herb_ingre_dict(self):
        """初始化草药-成分字典（修复版）"""
        return {
            k: [v_2 for v_2 in v if v_2 in self.ingredients]
            for k, v in self.herb_ingre_dict_all.items()
            if len(v) > 0
        }

    def precess_data(self):
        """
        预处理数据，去除无效行。
        """
        data_left = self.data_all.dropna(subset=['ingredients_id'], axis=0)
        return data_left

    def herb_ingre_dict_all(self):
        """
        返回草药-成分字典。
        """
        return qc_herb(self.data)



    def herb_ingre_dict(self, ingredients_target_dict):
        """
        根据成分-目标字典筛选草药-成分字典。
        """
        self.herb_ingre_dict = {
            k: [v_2 for v_2 in v if v_2 in ingredients_target_dict.keys()]
            for k, v in self.herb_ingre_dict_all.items()
        }
        self.herb_ingre_dict = {
            k: v for k, v in self.herb_ingre_dict.items() if len(v) != 0
        }
        return self.herb_ingre_dict
    
    def get_valid_herbs(self):
        """返回所有有效草药的列表"""
        return self.herbs

    def herb_ingretargets(self, herb, ingredients_target_dict):
        """
        返回草药对应的目标列表。
        """
        ingre_target_list = [ingredients_target_dict[ingre]
                             for ingre in self.herb_ingre_dict_all[herb] if ingre in ingredients_target_dict.keys()]
        if len(ingre_target_list) == 0:
            return None
        else:
            return list(set((reduce(lambda x, y: x + y, ingre_target_list))))

    def herb_ingretargets_dic(self, ingredients_target_dict):
        """
        构建草药-目标字典。
        """
        self.herb_ingretargets_dic = defaultdict()
        for herb in list(self.herbs):
            self.herb_ingretargets_dic[herb] = self.herb_ingretargets(herb, ingredients_target_dict)
        return self.herb_ingretargets_dic

    def read_data2(self):
        """
        读取草药-成分数据文件。
        """
        data = pd.read_csv(self.filename, sep=',')
        data['ingredients_id'] = data['ingredients_id'].astype(str).str.replace(r'\.0$', '', regex=True)
        data['ingredients_id'] = 'I' + data['ingredients_id'].astype(str)
        data['herb_id'] = 'H' + data['herb_id'].astype(str)

        return data