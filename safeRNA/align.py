import os, subprocess
from dataclasses import dataclass
import numpy as np

#change into the desired directory to store the results of the sorting
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

def blast(query_dir, target_dir, outdir): 
    """
    siRNA ==> query
    ref transcriptome ==> subject
    align each siRNA to desired sequences
    """
    #loop through all siRNA files
    for target in os.listdir(target_dir):
        #loop through all target sequence files
        for query in os.listdir(query_dir):
            out = f"{query[:-6]}_vs_{target[:-6]}.blast"
            #run blast on the fasta selected
            subprocess.run(f"blastn -task blastn-short -word_size 21 -query {os.path.join(query_dir, query)} -subject {os.path.join(target_dir, target)} -out {os.path.join(outdir, out)}")
    return True    

"""blast(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\siRNA_seq",
    r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\blast_db",
    r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\out_dir")"""
 
@dataclass
class Summarhit:
    rna_ID : str
    total: int
    target_hit : list

def report(blast_dir, description = "Analyzes blast report results")-> list:
    hit_count = []
    for index, reports in enumerate(os.listdir(blast_dir)):
        
        #add the id of the RNA sequence analyzed to the list
        hit_count.append(Summarhit(rna_ID=reports[:-6], total=None, target_hit=[]))

        #open blast reports
        with open(os.path.join(blast_dir, reports)) as blast_report:
            blast_report = ''.join(blast_report.readlines())
            #split blast report by query
            sep_rep = blast_report.split("Query= ")
            #count the absence of hits
            hit_count[index].total = len(sep_rep)-1 - blast_report.count("***** No hits found *****")
        
        #if hits are found, puts the id of sequences hitted in a dictionnary
        if hit_count[index].total != 0:
            file_hit = reports.split('vs')[1][:-6].strip('_') 
            all_hit = [elt for elt in sep_rep if "***** No hits found *****" not in sep_rep]
            hit_count[index].target_hit.append(
                {file_hit : [hit.split('\n')[0] for hit in all_hit[1:]]
                 })
    
    #The idea is to output a report file, but so far I've procrastinated this part of the code :D
    print([i for i in hit_count], sep='\n')        
    return 0    
            
'report(r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\out_dir")'
