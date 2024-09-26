import pandas as pd
import geopandas as gpd
import numpy as np
from matplotlib import pyplot as plt
from dataclasses import dataclass
from shapely import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable

class Map_background:
    def __init__(self, path):
        self.france = gpd.read_file(path)

class Cases(Map_background):
    def __init__(self, path, path_cases):
        """
        Column in the dataframe :

        id : square identification number
        position : position coordinates of the square's angle
        contain_stations : name and position of the meteorological stations contained
        name_station : list of all station names
        neighbour : id of all neighbouring squares
        args_start_date : arguments of the linear function from which the migration start date of aphid is predicted
        """
        Map_background.__init__(self, path)
        self.df = pd.read_csv(path_cases, index_col="Unnamed: 0")
        self.df['position'] = gpd.GeoSeries.from_wkt(self.df['position'])
        self.df["id"] = self.df["id"].astype(int)
        self.df["neighbour"] = self.df["neighbour"].tolist()

        self.gdf = gpd.GeoDataFrame(self.df, geometry='position')

    def initialize_dataframe(self):
        return self.df, self.gdf

    def show(self):
        """
        Show all the squares on a map
        """

        fig = plt.figure(1)
        ax = fig.subplots(1,1)
        self.france.plot(ax=ax,alpha=0.4, color='grey')
        gpd.GeoDataFrame(self.df, geometry='position').plot(ax=ax, alpha = 0.1, color = "gray", edgecolor="black")

        plt.show()
    
    def show_no_origin_case(self):
        fig = plt.figure(1)
        ax = fig.subplots(1,1)
        self.france.plot(ax=ax,alpha=0.4, color='grey')

        df_temp = self.df.loc[self.df["intercept"]==-1]
        gpd.GeoDataFrame(df_temp, geometry='position').plot(ax=ax, alpha = 0.1, color = "cyan", edgecolor="black")

        plt.show()

    def calc(self, year):
        """
        calculate the day at which aphid migration starts for all squares where the temperature conditions are satisfied
        """
        previsions = []
        for _id in self.df.index:
            a, b, c, slope, intercept = self.df.loc[_id]["a"].astype(float), self.df.loc[_id]["b"].astype(float), self.df.loc[_id]["c"].astype(float), self.df.loc[_id]["slope"].astype(float), self.df.loc[_id]["intercept"].astype(float)
            pred = round(a*np.sin(b*year+c)+(slope*year+intercept))
            if pred == -2023:
                pred = 500
            previsions.append(pred)

        self.df[f"{year}"] = previsions
        self.gdf = gpd.GeoDataFrame(self.df, geometry='position')

        return self.df, self.gdf
    
    def migration_start_distribution(self, gdf, year, ax):
        self.gdf_out = gdf.loc[self.df["intercept"]==-1]
        self.gdf_in = gdf.loc[self.df["intercept"]!=-1]

        self.france.plot(ax=ax,alpha=0.4, color='grey')

        self.gdf_in.plot(ax=ax, column=f"{year}", alpha=1, legend=True, 
                      legend_kwds={"label":"Earlyness of migration in days after 1st January", "orientation":"vertical"}, 
                      cmap="jet", 
                      cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.05),
                      vmin = 110,
                      vmax = 203)
        
        self.gdf_out.plot(ax=ax, column=f"{year}", alpha=0.7, legend=True, color='gray', label='No migration starting from these squares')


@dataclass
class Square(Cases):
    ID : int
    shape : Polygon
    neighbour_id : list
    populations : list
    migration_start_day : int
    general_density : float
    urbanisation : float

    def update_populations(self, all_pop):
        self.populations = [population for population in all_pop if population.position.within(self.shape)]
        return self

    def update_density(self):
        self.general_density = sum([pop.local_density for pop in self.populations])
        if self.general_density != 0:
            self.general_density += self.general_density
        return self
    
    def set_treshold(self):
        for population in self.populations:
            population.treshold_date = self.migration_start_day

        return self

