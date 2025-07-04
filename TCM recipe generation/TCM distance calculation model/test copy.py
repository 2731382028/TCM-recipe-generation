from generate_objects import herb_obj, herb_info, ingredients_obj, herb_distance_obj, fangji, g_obj

#test the functions
def test_herb_ingre(herb_distance_obj, g_obj):
    dis1=herb_distance_obj.Ingredients.ingre_ingre_dis('I358', 'I98', g_obj, 'separation')
    print("ingre_ingre_dis:", dis1)
    dis_all=herb_distance_obj.Ingredients.ingre_ingre_dis_all('I358', 'I98', g_obj)
    print("ingre_ingre_dis_all:", dis_all)
    #

    dis_uni=herb_distance_obj.herb_herb_distance_uni('H336', 'H346', 'separation')
    print("herb_herb_distance_uni:", dis_uni)
    s = herb_distance_obj.herb_herb_dis_all('H336', 'H346')
    print("herb_herb_dis_all:", s)


def main():
    test_herb_ingre(herb_distance_obj, g_obj)
if __name__ == "__main__":
    main()

