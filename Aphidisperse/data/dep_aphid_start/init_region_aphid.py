import pandas as pd
import geopandas as gpd
import ast
from shapely import Polygon, Point, MultiPolygon
from matplotlib import pyplot as plt
import seaborn as sns 

pd.options.mode.copy_on_write = False

sns.set_theme(rc={'figure.figsize':(11.7,8.27)})

def load_dep(path : str) -> pd.DataFrame:
    """
    Load the French departments with production of prunus spp. species from 
    which Myzus persicae can origin.

    Source : https://www.franceagrimer.fr/content/download/65522/document/CC-FL-2020-Chiffres-cles_FL_2019.pdf
    """
    df_dep = pd.read_csv(path, sep=";")
    liste_dep = ["Maine-et-Loire", "Charente-Maritime", "Dordogne", "Gironde", "Lot-et-Garonne",
                 "Lot", "Tarn-et-Garonne", "Aveyron", "Hérault", "Aude", "Gard", "Bouches-du-Rhône",
                 "Vaucluse", "Drôme", "Ardèche", "Loire", "Rhône", "Isère", "Hautes-Alpes"
                 "Vosges", "Moselle", "Meuse", "Bas-Rhin", "Haut-Rhin", "Meurthe-et-Moselle", "Pyrénées-Orientales"]
    df_dep = df_dep[df_dep["nom"].isin(liste_dep)]

    return df_dep

def df_to_gdf(df : pd.DataFrame) -> pd.DataFrame:
    all_xy = [ast.literal_eval(region) for region in list(df["Geo Shape"])]

    all_polygon = []
    for i in range(0, len(all_xy)):
        if all_xy[i]['type']=='MultiPolygon':
            all_polygon.append(MultiPolygon(all_xy[i]["coordinates"]))
        
        if all_xy[i]['type'] == 'Polygon':
            all_polygon.append(Polygon(all_xy[i]["coordinates"][0]))

    df["Geo Shape"] = all_polygon
    return df

def prunus_df() -> pd.DataFrame:
    """
    Assign the size of cultures in hectares to each region according to 
    https://www.franceagrimer.fr/content/download/65522/document/CC-FL-2020-Chiffres-cles_FL_2019.pdf

    DONE BY HAND : no directly available data.

    Some department did not specify the number of ha of culture. Their ha was estimated from
    substracting culture surface from the most important regions and diving equally the rest
    amongst the departments left. It is therefore only an estimation.

    Example : 
    Total ha of culture = 1000 ha

    Top3 regions : 
        1- Lot (455)
        2- Drôme (200)
        3- Moselle (145)

        Rest = 1000-455-200-135
        Rest = 200

    Let's say there are 10 minor producing regions but no details are given on the amount of ha exploited.
    Then for each of these regions we will estimate the cultivated surface to 200/10 = 20ha.   
    """

    dep = ["Maine-et-Loire", "Charente-Maritime", "Dordogne", "Gironde", "Lot-et-Garonne",
                 "Lot", "Tarn-et-Garonne", "Aveyron", "Hérault", "Aude", "Gard", "Bouches-du-Rhône",
                 "Vaucluse", "Drôme", "Ardèche", "Loire", "Rhône", "Isère", "Hautes-Alpes"
                 "Vosges", "Moselle", "Meuse", "Bas-Rhin", "Haut-Rhin", "Meurthe-et-Moselle"]
    
    culture = [92, 92, 92, 92, 176, 92, 1578, 92, 92, 92, 92, 120, 120, 92, 92, 92, 92, 92, 92, 120, 92, 485, 194, 120, 635]

    dico = dict(zip(dep, culture))
    dico["Tarn-et-Garonne"] += 422
    dico["Gard"] += 1308
    dico["Pyrénées-Orientales"] = 2763
    dico["Bouches-du-Rhône"] += 1823
    dico["Drôme"] = 879
    dico["Lot-et-Garonne"] += 200
    dico["Aude"] += 200
    dico["Hérault"] += 200
    dico["Vaucluse"] += 200
    dico["Ardèche"] += 200
    dico["Isère"] += 200
    dico["Rhône"] += 200

    return dico

def assignation(df_dep : pd.DataFrame, stats : dict) -> pd.DataFrame:
    """
    Link together the previously made dataframes for department and hectares of culture
    """
    df_dep["ha_prunus"] = [0 for i in range(0, 106)]

    for index in range(0, len(df_dep.index)):
        df_dep.iat[index, 7] = stats[df_dep.iat[index, 3]]

    return df_dep

def main():
    df_dep = load_dep(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\EDA_1\Data\contours_dep.csv")
    stats = prunus_df()

    df_dep = assignation(df_dep, stats)
    
    df_dep = df_to_gdf(df_dep)
    print(df_dep)
    return 
    df_dep.to_csv(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Avengers\data\dep_aphid_start\dep_aphid_start.csv")

    return 0

main()