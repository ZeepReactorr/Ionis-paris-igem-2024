from seed_init_MSA import Alignment, AA_LIST, CODONS_LIST, NUC_LiST
from seed_init_codon_graph import CodonGraph

import numpy as np
from collections import Counter
import pandas as pd
import re

def init_alignement(path : str) -> tuple:
    """
    Initialize the MSA of sequences to be studied.
    Return the MSA under the form of a matrix, and then determine Shanon's entropy 
    for each column of the matrix.

    Example : 

    input : 
    Seq1 ACCTGC
    Seq2 ACTCGC
    Seq3 ACCTAC

    A possition specific scoring matrix is then determined for each codon with the formula 
    weight(codon) = (freq(codon) + 1)/(Number_seq + len_alphabet) :
    
        0           1           
    ACC 0.375       0.125
    ACT 0.25        0.125  
    TGC 0.125       0.25
    CGC 0.125       0.25
    TAC 0.125       0.25

    Then Shanon's entropy is calculated for each position with the formula:
    S = log2(len(alphabet) + Sum(weight_at_position_i * log2(weight_at_position_i))

        0           1
    S   0.16629     0.57193

    output : 
    >>> mat = 
    [['ACC', 'TGC'],
    ['ACT', 'CGC'],
    ['ACC', 'TAC']]

    df_shanon =             0           1
                Entropy     0.16629     0.07193
    """

    ali = Alignment(path)
    mat, dicomat = ali.init_codon()
    count_codon = ali.counter(mat)
    df = ali.calc_PSSM(count_codon, CODONS_LIST)
    df_shanon = ali.Shanon(df, CODONS_LIST)
    return mat, df_shanon

def init_consensus(mat : np.matrix) -> list :
    """
    Initialize the reference sequence that will be used. Can be either :
    - The consensus sequence of the MSA
    - The most recent sequence of the MSA

    return a list of the most conserved elements in the sequence. Example :
    
    input : 
    Seq1 ACCTGC
    Seq2 ACTCGC
    Seq3 ACCTGC

    output : 
    >>> ['ACC', 'TGC']
    """
    return [max(Counter(line), key=Counter(line).get) for line in mat.tolist()]

def init_newest(path : str):
    ali = Alignment(path)
    mat, dicomat = ali.init_codon()
    return [max(Counter(line), key=Counter(line).get) for line in mat.tolist()]


def liste_codon_calc(df_entropy : pd.DataFrame, consensus : list) -> tuple:
    """
    Generates the synonym codon list at the positions where they are not the most conserved.
    
    Example : 
    
    Input :
    df_entropy =
                0           1
    Entropy     0.16629     0.07193

    consensus = 'ACCTGC'

    output : 
    >>> ([['ACC', 'ACT'], ['TGC', 'CGC', 'TAC']], {0 : {'codon_consensus': 'ACC', 'codon_entropy' : 0.16629},
                                                   1 : {'codon_consensus': 'TGC', 'codon_entropy' : 0.07193}})

    """

    Graph = CodonGraph().graph
    gen_seq = []
    trace = {}

    for index, codon in enumerate(consensus):
        if df_entropy.iat[0, index] != df_entropy.max(axis=1).values[0]:
            gen_seq += [Graph[codon]]
            trace[len(gen_seq)-1] = {"codon_consensus" : codon, "codon_entropy" : df_entropy.iat[0, index]}

        else:
            gen_seq.append(codon)

    return gen_seq, trace

def seed_sequence_generation(all_seq : list, consensus : list, trace : dict) -> list:
    """
    Generates all the sequences from the consensus sequence with on each, one mutated
    codon at a position where the entropy indicates a lack of conservation.

    Example : 

    Input :
    consensus = 'ACCTGC' 
    all_seq = [['ACC', 'ACT'], ['TGC', 'CGC', 'TAC']]
    trace = {0 : {'codon_consensus': 'ACC', 'codon_entropy' : 0.16629},
             1 : {'codon_consensus': 'TGC', 'codon_entropy' : 0.07193}})

    output : [(sequence, index codon altered, trace from the previous dictionnary)]
    >>> seeds =[('ACCTGC', 0, {'codon_consensus': 'ACC', 'codon_entropy' : 0.16629}),
                ('ACTGC', 0, {'codon_consensus': 'ACC', 'codon_entropy' : 0.16629}), 
                ('ACCTGC', 1, {'codon_consensus': 'TGC', 'codon_entropy' : 0.07193}),
                ('ACCCGC', 1, {'codon_consensus': 'TGC', 'codon_entropy' : 0.07193}),
                ('ACCTAC', {'codon_consensus': 'TGC', 'codon_entropy' : 0.07193})]
    """
    
    seeds = []
    temp = []
    for index, codon in enumerate(all_seq):
        if type(codon) == list and codon != []:
            for elt in codon:
                temp.append(elt)
                temp += consensus[index+1:]

                if temp != consensus :
                    seeds.append((temp, index, trace[index]))

                temp = consensus[:index]

            temp.append(consensus[index])

        else:
            temp.append(codon)

    return seeds

def temp_output(seeds : list, seq_annot : str) -> __file__:
    with open(fr"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\temp\seeds_{seq_annot}.fasta", 'w', encoding='utf-8') as seed_file:
        for index, seed in enumerate(seeds):
            #columns = id protein | codon position | ref codon | Shannon's entropy
            seed_file.write(f">{seq_annot}_pred_{index}|{seed[1]}|{seed[2]['codon_consensus']}|{seed[2]['codon_entropy']}\n\
{''.join(seed[0])}\n")

def main(MSA_PATH : str,  seq_annot = "", newest_path = "", param = 'consensus') -> int:
    """
    Activates all previously mentionned functions.

    The parameter "consensus" (default) makes the baseline sequence on which substitutions occur be
    the consensus sequence of the MSA, while the parameter "newest" uses the most 
    recent sequence to do so. An additional path to the file containing the most recent
    sequence is then required.

    """
    mat, df_entropy = init_alignement(MSA_PATH)

    if param == "consensus" :
        consensus_codons = init_consensus(mat)
    elif param == "newest":
        consensus_codons = init_newest(newest_path)

    all_seq, trace = liste_codon_calc(df_entropy, consensus_codons)
    sequences = seed_sequence_generation(all_seq, consensus_codons, trace)
    temp_output(sequences, seq_annot)
    return 0

main(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\alignment_all_seq\ALL_seq_RNAPOL_ali.fasta", 
        "RNA_pol", 
        r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\alignment_all_seq\most_recent_BYV.fst",
        "newest")
    