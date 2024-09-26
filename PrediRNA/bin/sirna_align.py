from sirna_init import load_seq, sep_seq, LHrna
from sirna_sampling import sample_seq

import pandas as pd
import os
from Bio import Seq
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy import stats
import warnings

sns.set_theme(rc={'figure.figsize':(12,8)})
warnings.filterwarnings('ignore')

def init_siRNA(path : str
               ) -> list:
    """
    See documentation of sirna_init for details.

    Returns a list of the precursors RNA obtainable from a target sequence, 
    and then a list of all the possible siRNA that can be obtained in a plant
    using a sliding window.
    """
    return sep_seq(load_seq(path), None, complement=True, file_out=False)

def select_target(directory : str, 
                  latest = True
                  ) -> str:
    """
    Select the target sequence file to use and create.

    Returns the path of the chosen file.

    The file returned depends on the boolean parameter latest.
    Default value is True, meaning by default the file used is the
    latest available.

    In the following documentation, the file containing the sequences 
    we want to target with siRNA will be called the "target file".
    All the other files will be denominated as "wanted files", even those simulated.
    """
    if latest == True :
        print("================== You are using the latest sequencing data for your analysis. ==================\n")
    elif latest == False :
        print("================== You are using the most wanted sequencing data for your analysis. ==================\n")

    return os.path.join(directory, sorted(os.listdir(directory), key = lambda x:int(x.split('.')[0]), reverse=latest)[0])

def init_df(directory : str, 
            target : str
            ) -> pd.DataFrame:
    """
    Initialize a dataframe based on the number of sequences contained in the target file.

    Returns an empty dataframe with as columns the name of the files, and the sequences id
    of the target file sequences as the index.
    """
    with open(target) as file:
        count = ''.join(file.readlines()).count('>')
    
    return pd.DataFrame(columns=os.listdir(directory), index=[f"seq_{i}" for i in range(0, count)])

def align(attacker : LHrna, 
          seq_b : str, 
          simu=False
          ) -> int:
    """
    Returns the number of siRNA from the target file that aligns (i.e. that can attack) 
    on the wanted file sequence.

    Example input :
    <<< target = ATGCTGCATGCTGATCGTCGATGC
        si = ['TGCTGCATGCTGATCGTCGA', 'TGCTGCATGCTGATCGTCGAT', 'GCTGCATGCTGATCGTCGATG', 'TTGCATGCTGATCGTCGATGT']
    
    >>> n = len(['TGCTGCATGCTGATCGTCGA', 'TGCTGCATGCTGATCGTCGAT', 'GCTGCATGCTGATCGTCGATG'])
    >>> n = 3

    """
    if simu == True:
        return len([si for si in attacker.siARNs if str(si.sequence) in seq_b])
    return len([si for si in attacker.siARNs if str(si.sequence) in Seq.complement(seq_b)])

def count_si(attackers : list, 
             df : pd.DataFrame, 
             col_dir : str, 
             simu=False
             ) -> pd.DataFrame:
    """
    Loop through all possible target file sequences and all possible sequences found in
    the wanted files.
    The siRNA found to work are stored in a list and associated with the target file and 
    wanted file used to count the siRNA efficiency.

    Returns a dataframe filled with the count of siRNA obtained from a target file and efficient
    toward a given wanted file. 

    Example Input :
    <<< df =        wanted_0.fasta
            seq_0   NaN
    
    >>> n = align(seq_0, wanted_0.fasta)
    >>> df =        wanted_0.fasta
            seq_0   [3]
    """
    
    for index_att, lh in enumerate(attackers):
        for index_col, col in enumerate(df.columns):
            sequences = load_seq(os.path.join(col_dir, col))
            values = []
            for read in sequences:
                values.append(align(lh, read.seq, simu))

            df.iat[index_att, index_col] = values
            
    return df

def adjust_df(df : pd.DataFrame, 
              simu = False
              ) -> pd.DataFrame:
    """
    Adjust the format of the dataframe so that there is a column for each siRNA count
    from the target file toward all the sequences in the wanted files.

    This is made to ease further analysis.

    Example input :
    <<< df =        wanted_0.fasta
            seq_0   [3,2,4]
    
    >>> df =        wanted_0_1.fasta    wanted_0_2.fasta    wanted_0_3.fasta
            seq_0   3                  2                  4
    """
    new_frame = pd.DataFrame(index=df.index)
    for index_col, col in enumerate(df.columns):
        for index_ind, elt in enumerate(df.index):
            for index, num in enumerate(df.iat[index_ind, index_col]):
                if index_ind == 0:
                    new_frame[f"{col.split('.')[0]}_{index}"] = num

                else:
                    new_frame.loc[elt,f"{col.split('.')[0]}_{index}"] = num

    if simu == False:
        new_frame.loc["date"] = [int(date.split('_')[0]) for date in new_frame.columns]

    if simu == True:
        new_frame.loc["date"] = [int(date.split('_')[1].split('.')[0]) for date in new_frame.columns]

    return new_frame

def evolution_efficiency(df : pd.DataFrame, 
                         colorset : list,
                         ax1 : plt.Axes,
                         starting_date = 0, 
                         label = 'Reference',
                         ) -> pd.DataFrame:
    """
    We suppose that the decrease in siRNA effect over time is linear.
    Therfore we do a linear regression based on the counts made before to estimate
    the number of siRNA that will still be efficient after a given number of years.

    The program returns a dataframe containing:
    - coef : the slope of the linear regression curve
    - intercept : the intercept of the linear regression curve
    - the color in which it has been drawn to match future plots
    - the starting date of the regression
    - the finishing date of the regression
    - the type (if it is one of the reference curves or a simulated one). Default is reference type.

    By default, two regressions are displayed in all plots for references : the regression based on 
    empirical analysis from NCBI's datasets (blue), and the worst case scenario, where mutations occurs 
    specifically on the worst possible position to hinder siRNA efficiency (red).
    """
    LR = LinearRegression()
    Y_all = df.drop("date", axis=1)
    X = [[date+starting_date] for date in df["date"]]

    dico = {'coef' : [], 'intercept':[], "color":[], "start": [], "end":[], "type":[]}

    for index, col in enumerate(Y_all.columns):
        Y = [[elt] for elt in Y_all[col]]
        LR.fit(X, Y)
        ypred = LR.predict(X)
        
        ax1.plot(X, ypred, label=label, marker='o', color=colorset[index])

        dico['coef'].append(LR.coef_[0][0])
        dico['intercept'].append(LR.intercept_[0])
        dico['color'].append(colorset[index])
        dico["start"].append(df['date'].min()-df['date'].min())
        dico["end"].append(df["date"].max()-df['date'].min())
        dico["type"].append(label)

    return pd.DataFrame(data=dico)

def display_proportion(pred : pd.DataFrame, 
                       ax : plt.Axes, 
                       starting_date : int, 
                       rate = 3
                       ) -> None:
    """
    The linear regression analysis allow to determine the raw number of siRNA that will be efficient depending 
    on time. Here we determine the percentage of active siRNA toward the wanted sequences and in function of
    time.

    The parameters are :
    - pred : the dataframe containing the parameters of the linear regressions
    - ax : the ax on which to plot
    - starting date : the date at which the analysis starts
    - The rate : determines the step at which proportion is displayed (default : 3 years).
    """
    df = pd.DataFrame(data=[[pred.iat[index, 0]*(year+starting_date)+pred.iat[index,1] 
                                   for year in range(pred["start"].max(), pred["end"].min()+rate, rate)] 
                                   for index in pred.index],
                            columns = [year+starting_date for year in range(pred["start"].max(), pred["end"].min()+rate, rate)])
    
    df_value = df.copy(True)
    for index in df_value.index:
        df_value.iloc[index, :] = (df_value.iloc[index, :]/df_value.iloc[index, :].max())*100

    value = 2/len(df_value.index)

    dico_significance = {index : None for index in df_value.index if pred.iat[index, 5] == "Reference"}
    for index in df.index:
        if pred.iat[index, 5] == "Reference":
            for num in df.index:
                if pred.iat[num, 5] != "Reference":
                    stat, pvalue = stats.chisquare((df.iloc[num, :]/df.iloc[num, :].sum())*df.iloc[index, :].sum(), df.iloc[index, :])
                    
                    if pvalue > 0.05:
                        dico_significance[index] = False

    for index in df_value.index:
        ax.bar([x+value for x in df_value.columns], df_value.iloc[index, :], 
               width=1/pred.shape[0], color=pred.iloc[index, 2], align='center', alpha=1)
        
        value+=1/len(df_value.index)
    
    if list(dico_significance.values()) == [False]*len(dico_significance.values()):
        ax.plot([min(df.columns), min(df.columns), max(df.columns)+2, max(df.columns)+2], [100, 102, 102, 100], lw=1.5, c='k')
        ax.text((min(df.columns) + max(df.columns)+2)/2, 102, "Non-significant (pvalue >> 0.05)", ha='center', va='bottom', color='k')

    ax.set_ylim((0,110))
    ax.set_ylabel("% of siRNA efficient ")
    ax.set_xlabel("Time (y)")
    ax.set_title("Evolution of the efficiency of siRNA generated from p21\ntargeting the BYV sequence depending on time.")
    
    plt.show()
    

def calc_intercept(low_bound : int, 
                   f_val : float,
                   rate = 2.987
                   ) -> float:
    """
    Calculates the intercept of a linear function based on a known value at a given point.
    """
    return f_val + low_bound*21*rate

def show_wcs(df : pd.DataFrame, 
             ax1 : plt.Axes,
             rate : float = 2.987,
             ) -> pd.DataFrame :
    """
    Calculates and display the worst case scenario of mutations for interfering RNA based solutions.
    
    One mutation occuring can hinder at most 21 different sequences of siRNA. The worst possible happening
    are mutations occuring 21 nucleotides appart from one another.
    """
    _max = max([df.iloc[i, 0][0] for i, val in enumerate(df.index)])
    _bounds = (min([int(col.split('.')[0]) for col in df.columns]), max([int(col.split('.')[0]) for col in df.columns]))
    
    w_values = []
    for i in range(_bounds[0], _bounds[1]+1):
        w_values.append(_max)
        _max-=21

        if _max < 0:
            _max = 0
    
    ax1.plot([i for i in range(_bounds[0], _bounds[1]+1)], w_values, color='r', label="Worst case scenario")
    ax1.legend()

    dico = {'coef' : [-21*rate], 
            'intercept':[calc_intercept(_bounds[0], w_values[0])], 
            "color":['red'], 
            "start": [_bounds[0]-_bounds[0]], 
            "end":[_bounds[1]-_bounds[0]], 
            "type":["Reference"]
            }
    
    print(dico['coef'], _bounds[0], w_values[0])
    
    return pd.DataFrame(data=dico)

def display_base_evolution(df : pd.DataFrame, 
                           ax1 : plt.Axes, 
                           ) -> pd.DataFrame:
    """
    Calculate and send to the `evolution_efficiency` function the count of the siRNA
    efficient toward a given sequence depending on time.

    Display the scatterplot of the count of siRNA efficient toward a sequence, for the references.
    """
    color_set = []
    for index_col, col in enumerate(df.columns):
        for index_ind, elt in enumerate(df.index):
            for num in df.iat[index_ind, index_col]:
                ax1.scatter(int(col.split('.')[0]), num, color=f"C{index_ind}", marker='o')

                if f"C{index_ind}" not in color_set:
                    color_set.append(f"C{index_ind}")

    adjusted_df = adjust_df(df).T

    prop = evolution_efficiency(adjusted_df, list(set(color_set)), ax1)

    ax1.set_xlabel("time (y)")
    ax1.set_ylabel('Number of siRNA efficient')
    ax1.set_title("Evolution of the number of siRNA generated from p21\nefficient on the BYV between 1993 and 2023")
    return prop

def sort_frame(df : pd.DataFrame) -> pd.DataFrame:
    """
    Sort the columns of a dataframe in increasing order.
    """
    return pd.DataFrame(data=[df[col] for col in sorted(df.columns, key = lambda x: int(x.split('_')[1].split('.')[0]))]).T

def display_simulation(df : pd.DataFrame, 
                       ax1 : plt.Axes, 
                       starting_date = 0, 
                       number_attackers = 1
                       ) -> pd.DataFrame:
    """
    Display the simulated informations which are the result of stochastic sampling based on empirical observations.
    """
    color_set = []
    for index_col, col in enumerate(df.columns):
        for index_ind, elt in enumerate(df.index):
            for num in df.iat[index_ind, index_col]:    

                ax1.scatter(int(col.split('_')[1].split('.')[0]) + starting_date, 
                            num, color=f"C{index_ind + number_attackers}", alpha=0.05, marker='o')
                if f"C{index_ind+2}" not in color_set:
                    color_set.append(f"C{index_ind + number_attackers}")
               
    adjusted_df = adjust_df(df, simu=True).T
    prop = evolution_efficiency(adjusted_df, color_set, ax1, starting_date=starting_date, label="Simulated data")

    ax1.legend()
    return prop

def main_simulation(db_dir : str, 
                    ax1 : plt.Axes,
                    latest_target = True, 
                    simu_path = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\output\simulation",
                    starting_date = 1993,
                    number_attackers = 1
                    ) -> None:
    """
    The idea of simulation data is to give an idea of whether the sequence sampled are relevant.
    To do so the process is the following :

    From the target file (available on NCBI), we randomly generate a given number of sequences with
    a number of substitutions based on the substitution rate determined thanks to BEAST analysis.

    Typically the sequences are generated for the **same time as for the empirical analysis** to allow
    comparison.

    Then, **even if the file chosen is the latest**!! We compare the number and proportions of efficient 
    siRNA between empirical an simulated sequences, **as if they started on the same date, even though it
    is not the case if the latest sequence were chosen**.
    This is done to ease the distribution analysis. 
    """
    
    target_file = select_target(db_dir, latest=latest_target)
    
    all_attackers = init_siRNA(target_file)

    df_simu_target = init_df(simu_path, target_file)

    ref = str(all_attackers[0].mother_seq)
    sample_seq([ref[i:i+3]for i in range(0, len(ref), 3) if len(ref[i:i+3]) >2])
    df_simu_target = count_si(all_attackers, df_simu_target, simu_path, simu=True)

    sorted_df_simu_target = sort_frame(df_simu_target)
    print(sorted_df_simu_target)

    return display_simulation(sorted_df_simu_target,
                                ax1, 
                                starting_date=starting_date, 
                                number_attackers = number_attackers)

def main(db_dir : str, 
         latest_target : bool, 
         simu = True
         ):
    """
    Activating function for the processes.
    The simu parameter allows to simulate siRNA efficiency over time 
    when set on True (default False).
    """
    proportions = []
    target_file = select_target(db_dir, latest=latest_target)
    all_attackers = init_siRNA(target_file)
    df_ali = init_df(db_dir, target_file)

    df_ali = count_si(all_attackers, df_ali, db_dir)
    
    fig = plt.figure()
    ax1, ax2 = fig.subplots(2,1)
    fig.tight_layout(pad=4)

    proportions.append(display_base_evolution(df_ali, ax1))
    proportions.append(show_wcs(df_ali, ax1))

    if simu == True:
        proportions.append(main_simulation(db_dir, 
                        ax1,
                        latest_target= True,
                        starting_date = int(target_file.split('\\')[-1].split('.')[0]),
                        number_attackers=len(df_ali.index)))

    starting_date = int(target_file.split('\\')[-1].split('.')[0])
    df_param = pd.concat(proportions, ignore_index=True)
    display_proportion(df_param, ax2, starting_date)
    
    return 0

main(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\alignment_all_seq\db_p21", False, False)
