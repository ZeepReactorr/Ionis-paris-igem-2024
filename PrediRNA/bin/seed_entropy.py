from seed_init_MSA import Alignment

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_theme(rc={'figure.figsize':(12,8)})

PATH_MSA_BYV = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\alignment_all_seq\ALL_seq_RNAPOL_ali.fasta"
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'] 
NUC_LiST = ['A', 'T', 'C', 'G', '-']
CODONS_LIST = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

def init_MSA_RNApol(path):
    """
    Initialize Alignment object from fasta file.
    """
    return Alignment(path)

def PSSM(ali : Alignment, matrix  : np.matrix, liste_elt : list) -> pd.DataFrame:
    """
    Initialize the Position Specific Scoring Matrix for the loaded MSA.
    Specific documentation available in Alignment
    """

    dico = ali.counter(matrix)
    return ali.calc_PSSM(dico, liste_elt)

def calc_Shanon_entropy(ali : Alignment, df : pd.DataFrame, letter_list : list) -> pd.DataFrame:
    """
    Determine the Shanon's entropy for each position in the MSA.
    Specific documentation available in Alignment.
    """
    return ali.Shanon(df, letter_list)

def display_entropy(df_nuc : pd.DataFrame, df_pro : pd.DataFrame, df_cod : pd.DataFrame) -> None:
    """
    Display a heatmap of all the Shanon's Entropy for all positions.
    Graph are superposed one on top of another in a triple subplot along the Y axis.

    Top plot is the heatmap of Amino acids.
    Middle plot is the heatmap of Codons.
    Bottom plot is the heatmap of Nucleotides.

    Colorscales are normalized for better visualization.
    """
    fig = plt.figure(1)
    ax1, ax2, ax3 = fig.subplots(3, 1)
    
    sns.heatmap(data = df_pro, ax=ax1, robust=True, cmap = 'plasma_r', 
                cbar_kws={"label":"Shanon's Entropy\nat each position in MSA", "orientation":"vertical"}, 
                vmin=df_pro.min(axis=1), vmax=df_pro.max(axis=1))
    
    ax1.set_title("Heatmap of Shanon's Entropy for the BYV's MSA showing the relative conservation of amino acids\nat each position for the RdRP conserved sequence")

    sns.heatmap(data = df_cod, ax=ax2, robust=True, cmap = 'plasma_r', 
                cbar_kws={"label":"Shanon's Entropy\nat each position in MSA", "orientation":"vertical"}, 
                vmin=df_cod.min(axis=1), vmax=df_cod.max(axis=1))
    
    ax2.set_title("Heatmap of Shanon's Entropy for the BYV's MSA showing the relative conservation of codons\nat each position for the RdRP conserved sequence")

    sns.heatmap(data = df_nuc, ax=ax3, robust=True, cmap = 'plasma_r', 
                cbar_kws={"label":"Shanon's Entropy\nat each position in MSA", "orientation":"vertical"}, 
                vmin=df_nuc.min(axis=1), vmax=df_nuc.max(axis=1))
    
    ax3.set_title("Heatmap of Shanon's Entropy for the BYV's MSA showing the relative conservation of nucleotides\nat each position for the RdRP conserved sequence")

    fig.tight_layout()
    plt.savefig(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\mutation_prediction\figures\full_BYV_genome_shanon_entropy_and_correlation_heatmap.png")

def main_nuc(Ali : Alignment) -> pd.DataFrame:
    """
    Load and initialize the Dataframe of Nucleotide's entropy.
    """
    _nuc_Ali = Ali.init_nuc()    
    _PSSM = PSSM(Ali, _nuc_Ali[0], NUC_LiST)
    Sh_entropy = calc_Shanon_entropy(Ali, _PSSM, NUC_LiST)    
    return Sh_entropy

def main_pro(Ali : Alignment) -> pd.DataFrame:
    """
    Load and initialize the Dataframe of Amino acids's entropy.
    """
    _prot_Ali = Ali.init_prot()    
    _PSSM = PSSM(Ali, _prot_Ali[0], AA_LIST)
    Sh_entropy = calc_Shanon_entropy(Ali, _PSSM, AA_LIST)    
    return Sh_entropy

def main_codon(Ali : Alignment) -> pd.DataFrame:
    """
    Load and initialize the Dataframe of Codon's entropy.
    """
    _codon_Ali = Ali.init_codon()    
    _PSSM = PSSM(Ali, _codon_Ali[0], CODONS_LIST)
    Sh_entropy = calc_Shanon_entropy(Ali, _PSSM, CODONS_LIST)    
    return Sh_entropy

def main():
    """
    main activating function.

    Displays all plots.
    """
    _Ali = init_MSA_RNApol(PATH_MSA_BYV)

    nuc_entropy = main_nuc(_Ali)
    pro_entropy = main_pro(_Ali)
    cod_entropy = main_codon(_Ali)
    
    display_entropy(nuc_entropy, pro_entropy, cod_entropy)
    return 0

main()
