
from herb_ingre_tar import *
from herb_ingre_tar import Ingredients
from construct_network import *
from herb_distance_generation import *
import pandas as pd

# generate object

file_name = "data/Interactome2022.csv"
network = Construct_Network(file_name)
g_obj = network.G
# 初始化 Ingredients 对象
G_nodes_with_prefix = list(g_obj.nodes())  # 从网络对象获取原始节点
G_nodes = [node.lstrip('T') for node in G_nodes_with_prefix]  # 去除T前缀

# 初始化成分对象
ingredients_obj = Ingredients('data/stitch_com_entry_ids.csv')

# 显式初始化目标字典（必须步骤）
ingredients_obj.ingredients_target_dict(G_nodes)  # 传入处理后的节点列表


# 再次检查列名
print("Ingredients 对象初始化完成")
print("ingre_tar_dict 是否为空:", ingredients_obj.ingre_tar_dict is None)  # 应该输出 False
print("成分数量:", len(ingredients_obj.ingre_tar_dict))  # 应该输出大于 0 的值
print(f"ingre_tar_dict 是否为空: {len(ingredients_obj.ingre_tar_dict) > 0}")
print(f"ingre_tar_dict 示例: {list(ingredients_obj.ingre_tar_dict.items())[:5]}")

herb_obj = Herb('data/herb_ingredient_pairs.csv')
print("Herb 对象初始化完成")
print("草药数量:", len(herb_obj.herb_ingre_dict))  # 修改点：直接访问属性
print(f"herb_ingre_dict 是否为空: {len(herb_obj.herb_ingre_dict) > 0}")
print(f"herb_ingre_dict 示例: {list(herb_obj.herb_ingre_dict.items())[:5]}")
herb_distance_obj = Herb_Distance(g_obj, ingredients_obj, herb_obj)


from herb_herb_pairs import *
import pandas as pd

herb_info = Herb_Info('data/herbs_data.csv')
fangji = FangJi('data/prescription.txt')
fangji.herbid_frequency_dic(2, herb_info.herb_pinyin_dic)


