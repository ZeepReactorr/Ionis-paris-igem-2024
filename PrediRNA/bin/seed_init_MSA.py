from Bio import SeqIO, Seq
import pandas as pd
import numpy as np
from collections import Counter

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

class Alignment:
    """
    Parse the MSA of proteins and provides valuable attributes from implemented methods.
    Works for Nucleotides, Codons and Amino acids.

    Feature can be added as seen fit.
    """
    def __init__(self, 
                 filepath : str) -> None :
        """
        input : file in .fasta format ressembling this :

        >identification_of_nucleotide_sequence
        ATCTCGCTGATCGCTGTACGATCTGCTCGATGCTA...
        or
        >identification_of_amino_acids_sequence
        AERGCRGAIRNFGRGRDNPPGLYRAILESSTY

        Creates a list of SeqRecord instances from BioPython, with the
        associated variables id and seq allowing yo access the file content.
        """
        self.sequences = [elt for elt in SeqIO.parse(filepath, "fasta")]

    def init_nuc(self) -> tuple:
        """
        Initialize a matrix of sequences based on the alignment of the file.
        The matrix is transposed to reduce time complexity of operations, which
        consists in counting each element's occurence in the MSA.

        Returns a Numpy matrix of the MSA and a dictionnary linking directly the sequences with their ID.

        Example input :
        <<< self.sequences = [
                            Seq(id='seq_1', seq=Seq('ATGCC')),
                            Seq(id='seq_2', seq=Seq('TTGCC'))]

                            
        Output :
        >>> self.mat = A T
                       T T
                       G G
                       C C
                       C C
            self.dicomat = {'seq_1':'ATGCC',  'seq_2' = 'TTGCC'}
        """

        self.mat = np.matrix([list(elt.seq) for elt in self.sequences]).T
        self.dicomat = dict(zip([elt.id for elt in self.sequences], [list(elt.seq) for elt in self.sequences])) 
        return self.mat, self.dicomat

    def init_prot(self) -> tuple:
        """
        Initialize a matrix of sequences based on the alignment of the file.
        The matrix is transposed to reduce time complexity of operations, which
        consists in counting each element's occurence in the MSA.

        Returns a Numpy matrix of the MSA and a dictionnary linking directly the sequences with their ID.

        Example input :
        <<< self.sequences = [
                            Seq(id='seq_1', seq=Seq('ETAIL')),
                            Seq(id='seq_2', seq=Seq('ETALL'))]

                            
        Output :
        >>> self.mat = E E
                       T T
                       A A
                       I L
                       L L

            self.dicomat = {'seq_1':'ETAIL',  'seq_2' = 'ETALL'}
        """

        self.mat = np.matrix([list(Seq.translate(elt.seq, gap='-')) for elt in self.sequences]).T
        self.dicomat = dict(zip([elt.id for elt in self.sequences], [list(Seq.translate(elt.seq, gap='-')) for elt in self.sequences])) 
        return self.mat, self.dicomat
    
    def init_codon(self) -> tuple:
        """
        Initialize a matrix of sequences based on the alignment of the file.
        The matrix is transposed to reduce time complexity of operations, which
        consists in counting each element's occurence in the MSA.

        Returns a Numpy matrix of the MSA and a dictionnary linking directly the sequences with their ID.

        Example input :
        <<< self.sequences = [
                            Seq(id='seq_1', seq=Seq('ATGCCT')),
                            Seq(id='seq_2', seq=Seq('TTGCCT'))
                            ]

                            
        Output :
        >>> self.mat = ATG TTG
                       CCT CCT
            self.dicomat = {'seq_1': ['ATG', 'CCT'],  'seq_2' = ['TTG', 'CCT']}
        """

        self.mat = np.matrix([[str(elt.seq[i:i+3]) for i in range(0, len(elt.seq)-2, 3)] for elt in self.sequences]).T
        self.dicomat = dict(zip([elt.id for elt in self.sequences], [[str(elt.seq[i:i+3]) for i in range(0, len(elt.seq)-2)] for elt in self.sequences])) 
        return self.mat, self.dicomat

    def counter(self, matrix:np.matrix) -> dict:
        """
        Count the occurences of each elements given a position in the MSA.
        
        Returns a dictionnary of the number of each element at their position.

        Example input :
        nucleotides:
        <<< self.mat = A T
                       T T
                       G G
                       C C
                       C C

        >>> self.count_pro = {0 : {A:1, T:1, C:0, G:0},
                              1 : {A:0, T:2, C:0, G:0},
                              2 : {A:0, T:0, C:0, G:2},
                              3 : {A:0, T:0, C:2, G:0},
                              4 : {A:0, T:0, C:2, G:0}}

        The idea for amino acid count is the same as nucleotide's.
        Codons :
        <<< self.mat = ATG TTG
                       CCT CCA      

        >>> self.count_pro = {0 : {ATG:1, CCT:0, ..., TTG:1},
                              1 : {ATG:0, CCT:2, ..., TTG:0}}  
        """
        self.count_pro = {index : dict(Counter(seq)) for index, seq in enumerate(matrix.tolist())}

        for key, value in self.count_pro.items():            
            self.count_pro[key] = value

        return self.count_pro
    
    def calc_PSSM(self, dico_count : dict, liste_elt : list) -> pd.DataFrame:
        """
        Determines a Position Specific Scoring Matrix (PSSM) for the MSA given in input.
        Associate a score to the chances of another coding element appearing instead of the 
        consensus nucleotide/codon/amino-acid.

        The formula used to determine this PSSM is the following :
        weight(elt) = (freq(elt) + 1)/(Number_seq + len_alphabet)

        Returns a Dataframe with the alphanet as index and the positions as columns,
        associating both to give the element position specific score for the MSA.

        Example input :
        Codons:
        <<< Seq1 [ACC, TGC]
            Seq2 [ACT, CGC]
            Seq3 [ACC, TAC]

        >>>         0           1           
                ACC 0.375       0.125
                ACT 0.25        0.125  
                TGC 0.125       0.25
                CGC 0.125       0.25
                TAC 0.125       0.25
        """

        self.df = pd.DataFrame(data=np.zeros((len(liste_elt), len(dico_count))), index=liste_elt, columns=range(0, len(dico_count)))

        for key, value in dico_count.items():                       
            for nuc, num in value.items():
                self.df[key][nuc] += num
        
        for col in self.df.columns:
            self.df[col] = (self.df[col] + 1)/(sum(self.df[col]) + len(liste_elt))

        return self.df
    
    def Shanon(self, df:pd.DataFrame, letter_list : list) -> pd.DataFrame:
        """
        Calculates Shanon's entropy for the MSA, i.e. the quality of information transfer
        overtime.
        Good information transfer means high value of Entropy and no substitution.
        Bad information transfer means low value of Entropy and depends on the number of substitutions 
        at a given position.

        Shanon's entropy is calculated for each position in the sequence with the formula:
        S = log2(len(alphabet) + Sum(weight_at_position_i * log2(weight_at_position_i))

        Returns a dataframe containing the value of Shanon's entropy for each position in the sequence.
        Highlight the positions where substitutions occur.

        Example input:
        Codons:
        <<<         0           1           
                ACC 0.375       0.125
                ACT 0.25        0.125  
                TGC 0.125       0.25
                CGC 0.125       0.25
                TAC 0.125       0.25        

        >>>     df_shanon =             0           1
                            Entropy     0.16629     0.07193

        Here for example, there is one codon substitution in position 0 and there are two in position 1.
        Accordingly the entropy is higher in position 0 and lower in position 1.
        """
        sh_list = {col : np.log2(len(letter_list)) + np.sum(df[col]*np.log2(df[col])) for col in df.columns}

        return pd.DataFrame(data=sh_list, index=["Entropy"])
