
from proximity_key import *
from collections import defaultdict

class Herb_Distance:
    def __init__(self, G_obj, Ingredients_obj, Herb_obj):
        self.G = G_obj
        self.Ingredients = Ingredients_obj
        self.herb_obj = Herb_obj
        # 调用 Ingredients 和 Herb 对象的方法进行初始化
        self.Ingredients.ingredients_target_dict(list(self.G.nodes))
        self.herb_obj.herb_ingretargets_dic(self.Ingredients.ingre_tar_dict)
        print("herb_ingre_dict 示例:", list(self.herb_obj.herb_ingre_dict.items())[:5])


    def herb_herb_length_dict(self, nodes_from, nodes_to, distance_method):
        length_dict = Sets_Lengths(nodes_from, nodes_to).ingre_length(self.length_fuc, distance_method)
        return length_dict

    def herb_herb_dis(self, herb_from, herb_to, distance_method, distance_method_herb_list):
        if any([herb not in self.herb_obj.herb_ingre_dict.keys() for herb in [herb_from, herb_to]]):
            print(f'草药 {herb_from} 或 {herb_to} 不在 herb_ingre_dict 中')
            print(f"当前 herb_ingre_dict 中的草药ID: {list(self.herb_obj.herb_ingre_dict.keys())[:10]}")  # 打印前10个草药ID
            return None
        else:
            nodes_from = self.herb_obj.herb_ingre_dict[herb_from]
            nodes_to = self.herb_obj.herb_ingre_dict[herb_to]
            length_dict = self.herb_herb_length_dict(nodes_from, nodes_to, distance_method)
            saved_lengthes_dict = defaultdict()

            for distance_method_herb in distance_method_herb_list:
                dis_obj = Network_Distance(nodes_from, nodes_to, length_dict)
                distance = dis_obj.network_distance(distance_method_herb)
                saved_lengthes_dict[distance_method_herb] = distance

            distances = {'ingre_method': distance_method,
                         'two_level': {'length_dict': length_dict,
                                       'distances': saved_lengthes_dict}}
            return distances

    def length_fuc(self, nodes_from, nodes_to, distance_method):
        distance = self.Ingredients.ingre_ingre_dis(nodes_from, nodes_to, self.G, distance_method)
        return distance

    def herb_herb_distance_uni(self, herb_from, herb_to, distance_method):
        if any([herb not in self.herb_obj.herb_ingretargets_dic.keys() for herb in [herb_from, herb_to]]):
            print(f'herb {herb_from} or {herb_to} not in herb_ingretarget dictionary')
            return None
        else:
            nodes_from = self.herb_obj.herb_ingretargets_dic[herb_from]
            nodes_to = self.herb_obj.herb_ingretargets_dic[herb_to]
            length_dict = Sets_Lengths(nodes_from, nodes_to).target_lengths(self.G)
            dis_obj = Network_Distance(nodes_from, nodes_to, length_dict)
            distance = dis_obj.network_distance(distance_method)
            distances = {'ingre_method': distance_method,
                         'one_level': {'union': distance}}

            return distances

    def herb_herb_dis_all(self, herb_from, herb_to):
        method_list_ingre = ['separation', 'closest', 'shortest', 'kernel', 'center']
        method_list_herb = ['separation', 'closest', 'shortest', 'kernel', 'center']
        dis_dict = defaultdict()
        for method_ingre in method_list_ingre:
            values_two_level = self.herb_herb_dis(herb_from, herb_to,
                                                 method_ingre, method_list_herb)
            values_one_level = self.herb_herb_distance_uni(herb_from, herb_to, method_ingre)

        # 空值保护
            if values_two_level is None or values_one_level is None:
                print(f"草药对 {herb_from}-{herb_to} 计算失败，跳过")
                continue

            dis_dict[method_ingre] = {'two_level': values_two_level['two_level'],
                                      'one_level': values_one_level['one_level']}

        return dis_dict

    def generator_result(self, herb_pairs_list):
        herb_pairs_distances = defaultdict()
        n = 0
        k = 1
        for herb_pairs in herb_pairs_list:
            herb1, herb1_name, herb2, herb2_name, frequency = herb_pairs
            try:
                distances = self.herb_herb_dis_all(herb1, herb2)
                print(f'yes, herb pairs {herb1} and {herb2} are successful')
                for ingre_method in distances.keys():
                    dict_1 = distances[ingre_method]['two_level']['distances']
                    dict_1.update({
                        'Ingredient-level distance type': ingre_method,
                        'herb1': herb1,
                        'herb1_name': herb1_name,
                        'herb2': herb2,
                        'herb2_name': herb2_name,
                        'frequency': frequency
                    })
                    n += 1
                    herb_pairs_distances[n] = dict_1
                k += 1
                print(f'this is the {k}th successful pairs')
            except:
                continue
        return pd.DataFrame.from_dict(herb_pairs_distances, orient='index')

