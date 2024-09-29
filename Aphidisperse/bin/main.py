from init import Cases, Square, Map_background
from init_populations import Population, Culture

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import geopandas as gpd
import random as rd
from shapely import Point, Polygon, MultiPolygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import animation
import itertools
import numpy as np

sns.set_theme(rc={'figure.figsize':(16,9)})

PATH_MAP_FRANCE = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Aphidisperse\data\map\shapefile\FRA_adm0.shp"
PATH_DF_CASES = r'C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Aphidisperse\data\cases_with_weather\cases_wstart.csv'
YEAR = 2023

df, gdf = Cases(PATH_MAP_FRANCE, PATH_DF_CASES).calc(YEAR)

def init_points(path : str) -> list:
    """
    Initialize the starting position of aphid's population based on their
    hibernation behaviour seeking prunus trees.

    In case point visualization is needed : 
    
    fig = plt.figure(1)
    ax= fig.subplots(1)
    Map_background(PATH_MAP_FRANCE).france.plot(ax=ax, alpha=0.4, color='gray')
    gdf_points = gpd.GeoDataFrame(data={"geometry" : all_points})
    gdf_points.plot(ax=ax, alpha=0.8, color='red', markersize=10)
    plt.show()

    """

    all_points = []
    passed_region = []
    df = pd.read_csv(path)
    df["Geo Shape"] = gpd.GeoSeries.from_wkt(df["Geo Shape"])
    gdf = gpd.GeoDataFrame(df, geometry="Geo Shape")

    for index in gdf.index:
        if gdf.iat[index, 4] in passed_region:
            continue

        else:          
            passed_region.append(gdf.iat[index, 4])
            container = gdf.iat[index, 2]
            minx, miny, maxx, maxy = container.bounds

            all_points += [Point(rd.uniform(minx, maxx), rd.uniform(miny, maxy)) for i in range(0, gdf.iat[index, 8]//10)]
    
    return all_points
    
def create_aphid_population_start(_aphid_position_list : list) -> list:
    """
    Initialize the aphids population's positions on the map
    """
    return [Population(position, rd.uniform(0.1, 0.3), 500) for position in _aphid_position_list]

def initialize_cultures(_cultures_position_list : list) -> list:
    """
    Initialize the wanted culture position on the map
    """
    return [Culture(position) for position in _cultures_position_list]

def initialize_squares(_square_df : pd.DataFrame) -> list:
    """
    Initialize the status of the square on the map
    """
    return [Square(ID=df.iloc[line]["id"], shape=df.iloc[line]["position"],
                   neighbour_id=df.iloc[line]["neighbour"], populations=[], 
                   migration_start_day=df.iloc[line][str(YEAR)], general_density=0, urbanisation=0) for line in _square_df.index]

def display_starting_state(all_squares:list, cultures_list:list):
    """
    Display the starting stage of the simulation
    """

    fig = plt.figure(figsize=(16,9))
    ax1= fig.subplots(1, 1)
    Cases(PATH_MAP_FRANCE, PATH_DF_CASES).migration_start_distribution(gdf, 2023, ax1)
    ax1.set_ylabel("latitude")
    ax1.set_xlabel("longitude")
    ax1.set_title(f"Migration starting day prediction of the aphid Myzus Persicae for {YEAR}", fontweight="bold", size=18)
    fig.tight_layout()
    plt.savefig(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Aphidisperse\output\starting_state_temp.png")

    fig = plt.figure(figsize=(16,9))
    ax2= fig.subplots(1, 1)
    Map_background(PATH_MAP_FRANCE).france.plot(ax=ax2, alpha=0.4, color='gray')
    gpd.GeoDataFrame({"geometry":[culture.position for culture in cultures_list]}).plot(ax=ax2, alpha=0.8, color='red', markersize=10)
    gpd.GeoDataFrame({"geometry":[square.shape for square in all_squares], 
                      "density":[square.general_density for square in all_squares]}).plot(ax=ax2, alpha=0.5, column="density", edgecolor="black", legend=True, 
                                                                                          legend_kwds={"label":"Density of aphid's population", 
                                                                                                       "orientation":"vertical"},
                                                                                          cmap="jet", 
                                                                                          cax = make_axes_locatable(ax2).append_axes("right", 
                                                                                                                                     size="5%", 
                                                                                                                                     pad=0.01),
                                                                                          vmin=0, 
                                                                                          vmax=10)

    ax2.set_title("Initial density distribution of aphid populations and culture positions", fontweight="bold", size=18)
    ax2.set_ylabel("latitude")
    ax2.set_xlabel("longitude")

    fig.tight_layout()
    plt.savefig(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Aphidisperse\output\starting_state_aphids.png")

def culture_positions(cultures_list:list):
    """
    Creates a GeoDataframe containing the cultures position
    """

    return gpd.GeoDataFrame({"geometry":[culture.position for culture in cultures_list]})

def update_all(t:int, aphid_population_list:list, all_squares:list):
    """
    Update all the population and square parameters following the new iteration
    """
              
    aphid_population_list = [pop.update() if t >= pop.treshold_date else [pop] for pop in aphid_population_list ]

    aphid_population_list = list(itertools.chain.from_iterable(aphid_population_list))

    all_squares = [square.update_populations(aphid_population_list) for square in all_squares]
    all_squares = [square.update_density() for square in all_squares]
    all_squares = [square.set_treshold() for square in all_squares]

    return aphid_population_list, all_squares

def display(t, ax, all_squares:list, gdf_cultures:gpd.GeoDataFrame):
    """
    Display the background, the cultures and the grid with the relative density of aphids population in each square
    """

    Map_background(PATH_MAP_FRANCE).france.plot(ax=ax, alpha=0.4, color='gray')
    gdf_cultures.plot(ax=ax, alpha=0.8, color='red', markersize=20)
    gpd.GeoDataFrame({"geometry":[square.shape for square in all_squares], 
                      "density":[square.general_density for square in all_squares]}).plot(ax=ax, alpha=0.5, column="density", edgecolor="black", legend=True, 
                                                                                          legend_kwds={"label":"Density of aphid's population", 
                                                                                                       "orientation":"vertical"},
                                                                                          cmap="jet", cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.01),
                                                                                          vmin=0, vmax=10)
    ax.title.set_text(f"Day {t} of simulation")

def loop(t):
    """
    Loop function linked to the animation. 
    Generates the images then compiled in an animation .mp4

    The re-definition of aphid_population_list and all_squares variables as global variable is necessary
    for them to be updated continuously and not go back to their initial definition state.
    """

    print(f"time : {t}")
    global aphid_population_list, all_squares

    _main_fig.clf()
    _main_ax = _main_fig.subplots(1,1)
    
    aphid_population_list, all_squares = update_all(t, aphid_population_list, all_squares)
    display(t, _main_ax, all_squares, gdf_cultures)

    print(len(aphid_population_list))
    return aphid_population_list, all_squares, gdf_cultures

def main(aphid_position_list:list, cultures_position_list:list):
    """
    _aphid_population and _all_squares ==> not global 

    aphid_populations and all_squares ==> global variables
    """
    
    global _main_fig, gdf_cultures, aphid_population_list, all_squares
    
    aphid_population_list = create_aphid_population_start(aphid_position_list)

    cultures_list = initialize_cultures(cultures_position_list)
    all_squares = [square.update_populations(aphid_population_list) for square in initialize_squares(df)]
    all_squares = [square.update_density() for square in all_squares]
    all_squares = [square.set_treshold() for square in all_squares]

    gdf_cultures = culture_positions(cultures_list)

    display_starting_state(all_squares, cultures_list)
    
    _main_fig = plt.figure(1)

    plt.rcParams['animation.ffmpeg_path'] = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\ffmpeg-master-latest-win64-gpl\bin\ffmpeg.exe"   
    FFwriter = animation.FFMpegWriter(10)

    anim = animation.FuncAnimation(_main_fig, loop, frames=270, interval=365, repeat=False)
    anim.save(filename=r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Aphidisperse\output\diffusion.mp4", writer=FFwriter)
    
main(init_points(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Model_prediction\Aphidisperse\data\dep_aphid_start\dep_aphid_start.csv"), [Point(4, 48)])
