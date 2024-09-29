import os
from Bio import SeqIO, Seq
from dataclasses import dataclass

#change into the desired directory to store the results of the sorting
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

@dataclass
class LHrna:
    mother_seq:str
    frame:int
    mother_seqid:str
    siARNs : list

@dataclass
class SIrna:
    sequence:str
    ID : str
    mother_id : str

def load_seq(file:__file__, 
             description="Load the file's sequences in a list"
             ) -> list:
    """
    Parse the fasta file and create Seq objects into a list of the double stranded RNA precursor
    """
    return [seq for seq in SeqIO.parse(file, "fasta")]

def sep_seq(seq_list : list, 
            out_dir : str, 
            complement : bool = False, 
            file_out : bool = True, 
            ) -> list:
    
    """
    Create a list of objects containing all possible siRNA sequence from a double stranded RNA precursor.
    Uses a sliding window approach to decompose the origin sequence into as many 21 nucleotides sequences
    as possible.

    The length of the siRNA is fixed at 21.
    The parameter complement is an optional parameter (default False) creating the complement siRNA to the original sequence.
    The parameter file_out is an optional parameter (default True). When set on False, instead of a .fasta file
    the output is a list of LHrna objects.

    The fasta format is as follows :
    >{name_seq}_{number_siRNA}
    ATCTGCTAG...

    Example:
    <<< SeqRecord(id=Seq_0, seq(ATCGATCGTAGCTAGCGTGGCT), ...)

    >>> LHrna(mother_seq=seq(ATCGATCGTAGCTAGCGTGGCT), mother_seqid=Seq_0,
                    siARNs =[SIrna(sequence = 'ATCGATCGTAGCTAGCGTGGC', ID = 'Seq_0_0', mother_id = 'Seq_0'),
                             SIrna(sequence = 'TCGATCGTAGCTAGCGTGGCT', ID = 'Seq_0_1', mother_id = 'Seq_0')
                            ]
             )
    """
    #

    if complement == False:
        liste_rna = [LHrna(mother_seq=read.seq, frame=21, mother_seqid=read.id,
                        siARNs=[SIrna(sequence=read.seq[index:index+21],
                                    ID=f"{read.id}_{index}", mother_id=read.id) 
                                    for index in range(0, len(read.seq)) if len(read.seq[index:index+21]) > 20])
                    for read in seq_list]
    
    if complement == True :
        liste_rna = [LHrna(mother_seq=Seq.complement(read.seq), frame=21, mother_seqid=read.id,
                        siARNs=[SIrna(sequence=Seq.complement(read.seq[index:index+21]),
                                    ID=f"{read.id}_{index}", mother_id=read.id) 
                                    for index in range(0, len(read.seq)) if len(read.seq[index:index+21]) > 20])
                    for read in seq_list]

    file_list = []

    #create a fasta file per sequence targeted containing all siRNA sequences possible
    if file_out == True:
        for LH in liste_rna:
            all_rna = ''.join([f">{rna.ID}\n{rna.sequence}\n" for rna in LH.siARNs])
            file_list.append(f"siRNA_from_{LH.mother_seqid.split('/')[0]}.fasta")
            with open(os.path.join(out_dir, f"siRNA_from_{LH.mother_seqid.split('/')[0]}.fasta"), 'w') as file:
                file.write(f">{LH.mother_seqid}\n{LH.mother_seq}\n{all_rna}")
                
        return 0
    else:
        return liste_rna

