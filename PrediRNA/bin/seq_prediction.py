from seed_main import init_alignement

import re, math, warnings, sys
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import math
from collections import Counter

sns.set_theme(rc={'figure.figsize':(12,8)})
warnings.filterwarnings('ignore')

def parser(path : str) -> pd.DataFrame:
    """
    Parse fasta file with Regex.

    re.findall(r">(.*?)(?=.<|\n)", seeds) ==> Finds all elements between a > and \n
    re.findall(r"[ATCG]+", seeds) ==> Finds all string containing ATCG elements.

    returns a dataframe with following columns :
    >>> index   id   Seq    codon_position   ref_codon   Shannon_entropy

    index : number of the sequence
    id : identification of the sequence
    Seq : Sequence
    codon_position : position of the mutation
    ref_codon : Codon in the consensus sequence
    shannon_entropy : Entropy at the mutation's position.
    
    """
    with open(path, 'r', encoding='utf-8') as seeds:
        seeds = ''.join(seeds.readlines())
    
    all_id = re.findall(r">(.*?)(?=.<|\n)", seeds)
    seeds = re.sub(r">(.*?)(?=.<|\n)", "", seeds)
    all_seq = re.findall(r"[ATCG]+", seeds)

    dico_all = {'id_prot' : [], 'Seq': [], 'codon_position':[], 'ref_codon': [], 'Shannon_entropy':[]}
    for index, value in enumerate(all_seq):
        id_prot, codon_position, ref_codon, shannon_entropy = all_id[index].split('|')

        dico_all['id_prot'].append(id_prot)
        dico_all['codon_position'].append([int(codon_position)])
        dico_all['ref_codon'].append(ref_codon)
        dico_all['Shannon_entropy'].append(float(shannon_entropy))
        
        dico_all['Seq'].append([str(value[i:i+3]) for i in range(0, len(value)-2, 3)])

    return pd.DataFrame(data=dico_all)

def poisson_conserved_codon(S:list, 
                            y : int, 
                            l = 0.0224, 
                            multiple_S = False,
                            sub = 10
                            ) -> list:
    """
    Model a Poisson probability law of parameter l calculated from the molecular rate of evolution
    from a BEAST analysis, giving the rate of substitutions per site per year.

    S is the entropy at a given position.
    y is the number of years passed from the starting year.
    l is the lambda parameter of the Poisson probability law.

    sub is the last number of substitutions we want to know the probability of.
        It is set by default at 10, meaning by default the probability calculated 
        ranges for 0 to 9 substitutions.
    """

    l = l*y
    all_s = []
    s = 0
    if type(S) != list:
        S = [S]

    for k in range(0, sub):
        temp = []
        for s in S:
            ls = l*s
            
            temp.append((ls**k/(math.factorial(k)))*math.exp(-ls))      

        all_s.append(math.prod(temp))
    return all_s

def show_multiple_poisson(df : pd.DataFrame, 
                          YEARS : int
                          ) -> None:
    """
    Display the distribution Poisson Law for sequences with more than
    one substitution that were predicted from the seeds.

    The X axis corresponds to the number of substitutions.
    The Y axis corresponds to the probability of k substitutions happening.
    """
    passed = []
    values = {}
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    for index, score in enumerate(df["coef_conservation"]):
        if score not in passed :
            passed.append(score)
        else:
            values[df.iat[index, 0]] = poisson_conserved_codon(score, YEARS, multiple_S=True)

    for key, value in values.items():
        ax1.plot(range(0, 10), value, label=key)

    ax1.set_xlabel('number of double substitutions')
    ax1.set_ylabel('Probability')
    ax1.set_title(f"Probability of a number of double substitutions happening in {YEARS} year")

    plt.savefig(fr"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\figures\event_substitution\event_sub_poisson_distribution_cycle2_{YEARS}y.png")
    

def show_poisson(df : pd.DataFrame, 
                 dico_num_cons : dict, 
                 YEARS : int
                 ) -> None:
    
    """
    Display the distribution Poisson Law for sequences with more than
    one substitution that were predicted from the seeds.

    The X axis corresponds to the number of substitutions.
    The Y axis corresponds to the probability of k substitutions happening.
    """

    dico_plot_S = {dico_num_cons[s] : poisson_conserved_codon(s, YEARS) for s in df['coef_conservation'].unique()}
    df[f"proba_seq_>=1_sub_{YEARS}y"] = pd.Series()

    #print(sum([sum(i[1:len(i)]) for key, i in dico_plot_S.items()]))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    for key, value in dico_plot_S.items():
        ax1.plot(range(0, 10), value)
        df[f"proba_seq_>=1_sub_{YEARS}y"].loc[df["num_conservation"] == key] = sum(value[1:len(value)])#/df.loc[df['num_conservation']==key, 'nS_per_position'].iloc[0]

    ax1.plot(range(0, 10), poisson_conserved_codon(1, YEARS, l=0.504), color='r', label='Non-conserved substitution in all the sequence')
    ax1.set_xlabel('number of substitutions')
    ax1.set_ylabel('Probability')
    ax1.set_title(f"Probability of a number of substitutions happening in {YEARS} year in one sequence for the entropy value in the sequence")
    ax1.legend(ncols=2)
    
    plt.savefig(fr"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\figures\event_substitution\event_sub_poisson_distribution_cycle1_{YEARS}y.png")

def coef_calc(df : pd.DataFrame, max_S : float) -> tuple:
    """
    Determines several important values for the rest of the program.
    
    dico_nS_per_position is a dictionnary counting the occurences of each
                         Shanon's entropy value.
    
    df["nS_per_position"] : number of occurences of each Shanon's entropy value.
    df["coef_conservation"] : Conservation coefficient given by the relative entropy ratio 
                              multiplied by the number of times this entropy value appears.
    """
    alt_df = df.drop_duplicates("codon_position", keep='first')
    dico_nS_per_position = dict(Counter(alt_df["Shannon_entropy"]))

    df["nS_per_position"] = [dico_nS_per_position[df.iat[index, 4]] for index in range(0, df.shape[0])]
    df["coef_conservation"] = [(max_S/s)*dico_nS_per_position[s] for s in df['Shannon_entropy']]

    dico = {elt : (index) for index, elt in enumerate(sorted(df["coef_conservation"].unique()))}
    df["num_conservation"] = [dico[s] for s in df["coef_conservation"]]

    return df, dico

def merge_seq(df : pd.DataFrame, 
              YEARS  : int
              ) -> pd.DataFrame:
    """
    Merge the sequences together based on one baseline sequence.
    Outputs a dataframe containing all possible sequence's substitutions,
    with all the relevant data.

    Example :

    seq_ref = ['AGC', 'TTG', 'GGA', 'AAT'] //This sequence has a substitution at index 1
    seq_merge = ['AGC', 'TTT', 'GGC', 'AAT'] //This sequence has a substitution at index 2

    >>> ns = ['AGC', 'TTG', 'GGC', 'AAT'] //List uncommon element have been combined.
    The output is therefore :

    >>> new_df = id_prot   Seq    codon_positions    ref_codon    Shannon_entropy   coef_conservation   num_conservation   proba_seq_>=1_sub_{YEARS}y
                 0         ns     [1,2]              [TTT, GGA]   [s1, s2]          [c1, c2]            [nc1, nc2]         P(s1)*P(s2)

    This program is VERY heavy. Time constraints forced a naive approach to offer the possibility of predicition.
    Time complexity is O(n^4) and space complexity is O(n^n).

    Improvement is needed to improve efficiency, despite the program working.
    """
    dico_all = {'id_prot' : [], 
                'Seq': [], 
                'codon_positions':[], 
                'ref_codon': [], 
                'Shannon_entropy':[], 
                "coef_conservation":[],
                "num_conservation":[],
                f"proba_seq_>=1_sub_{YEARS}y":[],
                }
    
    for index_ref, seq_ref in enumerate(df["Seq"]):
        for index_merge, seq_merge in enumerate(df["Seq"]):
            for pos in df.iat[index_merge, 2]:
                new_seq = [i for i in seq_ref]
                new_seq[pos] = seq_merge[pos]

            if len(list(set([df.iat[index_ref, 0], df.iat[index_merge, 0]]))) != 1:
                dico_all['Seq'].append(new_seq)
                dico_all['id_prot'].append(''.join([df.iat[index_ref, 0], df.iat[index_merge, 0]]))
                dico_all['codon_positions'].append(df.iat[index_ref, 2] + df.iat[index_merge, 2])
                dico_all['ref_codon'].append([df.iat[index_ref, 3], df.iat[index_merge, 3]])
                dico_all['Shannon_entropy'].append([df.iat[index_ref, 4], df.iat[index_merge, 4]])
                dico_all['coef_conservation'].append([df.iat[index_ref, 5], df.iat[index_merge, 5]])
                dico_all['num_conservation'].append([df.iat[index_ref, 6], df.iat[index_merge, 6]])
                dico_all[f"proba_seq_>=1_sub_{YEARS}y"].append(math.prod([df.iat[index_ref, 7], df.iat[index_merge, 7]]))

    new_df = pd.DataFrame(data = dico_all)
    return new_df

def fasta_forge(list_df : list, 
                YEARS : int
                ) -> __file__:
    """
    Creates a fasta file from the dataframe containing all the sequences.
    The format is as follows :

    >id_seq|proba sequence
    ATGCTGATCGATCGCGTAGCTAGCTAGC...

    The file is named after the number of years and number of cycles undergone.
    """
    for index, df in enumerate(list_df):
        with open(fr"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\output\{index+1}sub_{YEARS}y_RNA_pol.fasta", 'w', encoding='utf-8') as fastafile:
            for elt in range(0, len(df.index)):
                fastafile.write(f">{df.iat[elt,0]}|{df.iat[elt,7]}\n{''.join(df.iat[elt,1])}\n")
    return 0

def csv_forge(liste_df : list, 
              YEARS : int
              ) -> __file__ :
    """
    Creates a csv file from the dataframe containing all the sequences.
    The file is named after the number of years and number of cycles undergone.
    """
    for index, df in enumerate(liste_df):
        df.to_csv(fr"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\output\{index+1}sub_{YEARS}y_RNA_pol.csv")
    return 0

def sort_all_df(liste_df : list, 
                YEARS : int
                ) -> list:
    """
    Sort the order of elements in the dataframe by their probability of happening.
    """
    return [df.sort_values(by=f"proba_seq_>=1_sub_{YEARS}y", ascending=False) for df in liste_df]

def main(MSA_PATH : str, cycles : int, y : int, forge = True):
    all_df =[]
    max_S = init_alignement(MSA_PATH)[1].max(axis=1).iloc[0]
    df = parser(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\temp\seeds_RNA_pol.fasta")
    df = df[df['Shannon_entropy'] != df['Shannon_entropy'].max(axis=0)] # take out the non conserved substitutions
    df, dico_num_cons = coef_calc(df, max_S)

    show_poisson(df, dico_num_cons, y)

    all_df.append(df)
    
    for i in range(0, cycles):
        new_df = merge_seq(df, y)
        show_multiple_poisson(new_df, y)
        all_df.append(new_df)

    all_df = sort_all_df(all_df, y)

    if forge == True :
        fasta_forge(all_df, y)
        csv_forge(all_df, y)
        return 0
    else:
        return all_df

if __name__ == "__main__":
    path, cycles, n_years = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    for n in range(1, n_years+1):
        main(path, cycles, n)
        
"""main(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\alignment_all_seq\ALL_seq_RNAPOL_ali.fasta", 0, 30)"""