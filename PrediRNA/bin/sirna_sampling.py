import pandas as pd
import os
import random as rd

BASE_DIR = '\\'.join(os.path.dirname(os.path.abspath(__file__)).split('\\')[:-1])

def init_frames(path :str) -> list:
    return [pd.read_csv(os.path.join(path, file), index_col="Unnamed: 0") for file in os.listdir(path) if file.endswith('.csv') == True]

def binary_search_codon(df : pd.DataFrame, 
                        value : float, 
                        cumsum : float
                        ) -> pd.DataFrame:
    while df.shape[0] != 1:
        if df.iloc[df.shape[0]//2, 9] > value:
            df = df.iloc[0:df.shape[0]//2, :]
            
        else:
            df = df.iloc[df.shape[0]//2:df.shape[0], :]
            
    return df

def gen_sampling(all_df : 
                 list, batch = 10
                 ) -> list :
    
    liste_all_samples = []
    for index, df in enumerate(all_df):
        temp = []
        for i in range(0, batch):
            cumsum = df.iloc[:, 8].sum()
            liste_rand = [rd.uniform(0, cumsum) for i in range(0, 1+round(cumsum))]

            seqrow_list = [binary_search_codon(df, rand, cumsum) for rand in liste_rand]

            temp.append(pd.concat(seqrow_list))

        liste_all_samples.append(temp)
    
    return liste_all_samples

def write_year_0(ref, batch):
    path = os.path.join(BASE_DIR, "output","simulation")
    
    with open(os.path.join(path, "year_0.fasta"), 'w', encoding='utf-8') as file:
        for index in range(0, batch):
            file.write(f">initial_seq\n{ref}\n")

def write_year_1(df : pd.DataFrame) -> None:
    path = os.path.join(BASE_DIR, "output","simulation")
    
    with open(os.path.join(path, "year_1.fasta"), 'w', encoding='utf-8') as file:
        for index in range(0, len(df.index)):
            file.write(f">{df.iat[index, 0]}|{df.iat[index, 2]}|{df.iat[index, 3]}\n{''.join(df.iat[index, 1])}\n")

def fasta_forge(year : int, 
                df : pd.DataFrame, 
                batch : int, 
                seed_df : pd.DataFrame, 
                refseq : str
                ) -> __file__:
    """
    outputing fasta files with following id columns : 
    > <id batch sampling>|<id of the merged sequences>|<probability of sequence>
    <sequence>

    """
    path = os.path.join(BASE_DIR, "output","simulation")
    name = f"year_{year}.fasta"
    with open(os.path.join(path, name), 'w', encoding='utf-8') as file:
        for index in df.index:
            file.write(f">{df.iat[index, 0]}|{df.iat[index, 2]}|{df.iat[index, 3]}\n{''.join(df.iat[index, 1])}\n")
    return 0

def merge_seq(sample_id : int, 
              df : pd.DataFrame, 
              refseq : str
              ) -> pd.DataFrame:
    merged = refseq
    merged_id = df.iat[0,0]
    merged_proba = df.iloc[:, 8].prod(axis=0)

    for index_ind in range(1, len(df.index)):
        merged[int(df.iat[index_ind, 2][0])] = df.iat[index_ind, 1][int(df.iat[index_ind, 2][0])]
        merged_id += f"_{df.iat[index_ind, 0].split('_')[-1]}"
        
    res_df = pd.DataFrame(data = [[sample_id, merged, merged_id, merged_proba]], 
                          columns=["sample_id", "merged_seq", "merged_id", "merged_proba"], index=[sample_id])
    return res_df

def cross_df(all_samples : list, 
             refseq : str, batch, 
             seed_df : pd.DataFrame, 
             str_ref :str
             ) -> int:
    
    for index_year, liste_sample in enumerate(all_samples):
        temp = []
        for index_sample, sample in enumerate(liste_sample):
            if sample.shape[0] > 1:
                temp.append(merge_seq(index_sample, sample, refseq))
            else:
                df_1sub = pd.DataFrame(data = [[index_sample, sample.iat[0,1], sample.iat[0,0], sample.iat[0,8]]], 
                          columns=["sample_id", "merged_seq", "merged_id", "merged_proba"], index=[index_sample])
                temp.append(df_1sub)

            temp_df = pd.concat(temp)

        fasta_forge(index_year, temp_df, batch, seed_df, str_ref)
    return 0       

def adjust_df(liste_df : list) -> list:
    for df in liste_df:
        df["cumsum_proba"] = [df.iloc[0:index, 8].sum()+value for index, value in enumerate(df.iloc[:, 8])]
        df["codon_position"] = df["codon_position"].apply(lambda x: x.strip("[]").replace("'", "").split(", ") if x != '[]' else list())
        df["Seq"] = df["Seq"].apply(lambda x: x.strip("[]").replace("'", "").split(", ") if x != '[]' else list())

    return liste_df

def sample_seq(refseq : str, 
               path = os.path.join(BASE_DIR, "output","simulation"), 
               batch = 5
               ):
    """
    The sampling is done as follow :
    - A sequence
    """
    ref_refseq = "".join([i for i in refseq])
    liste_df = init_frames(path)
    liste_df = adjust_df(liste_df)

    liste_df = sorted(liste_df, key= lambda x:x.iloc[:, 8].sum())
    
    list_all_samples = gen_sampling(liste_df, batch)
    cross_df(list_all_samples, refseq, batch, liste_df[0], ref_refseq)
    
    return 0

