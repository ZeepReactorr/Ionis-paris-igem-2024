import os, subprocess, sys
from dataclasses import dataclass

@dataclass
class Summarhit:
    rna_ID : str
    total: int
    target_hit : list
    seq_dir: str

#change into the desired directory to store the results of the sorting
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

def report(blast_dir, db, rep_dir, seq_dir, description = "Analyzes blast report results")-> list:
    hit_count = []
    for index, reports in enumerate(os.listdir(blast_dir)):
        if reports.endswith(".blast"):
            #add the id of the RNA sequence analyzed to the list
            hit_count.append(Summarhit(rna_ID=reports[:-6], total=None, target_hit=[], seq_dir=seq_dir))

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
                all_hit = [elt for elt in sep_rep if "***** No hits found *****" not in elt]
                hit_count[index].target_hit += [hit.split('\n')[0] for hit in all_hit[1:]]
    
    #report file writing
    filename = f"{db}_report.txt"
    with open(os.path.join(rep_dir, filename), "a", encoding='utf-8') as report:
        for hit in hit_count:
            si_ori, seq_target = hit.rna_ID.split('_vs_')[0].split('siRNA_from_')[1], hit.rna_ID.split('_vs_')[1]
            link = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{seq_target}/"
            report.write(f"{hit.rna_ID}\t{si_ori}\t{seq_target}\t{db}\t{hit.total}\t{hit.target_hit}\t{hit.seq_dir}\t{link}\n")
    
    return 0

def blast(sirna_dir, db_dir, outdir, rep_dir): 
    """
    siRNA ==> query
    ref transcriptome ==> subject
    align each siRNA to desired sequences
    """

    #loop through all the databases
    roots = [root for root, dirs, files in os.walk(db_dir) if 'GCF' in root]
    fna = "cds_from_genomic.fna"

    for root in roots:
        db = root.split("\\")[3]

        for query in os.listdir(sirna_dir):
            out = f"{query[:-6]}_vs_{root[-15:]}.blast"

            subprocess.run(f"blastn -task blastn-short -word_size 21 -query {os.path.join(sirna_dir, query)} \
                        -subject {os.path.join(root, fna)} -out {os.path.join(outdir, out)}")
            
        report(outdir, db, rep_dir, root)

        for file in os.listdir(outdir):
            if file.endswith(".blast"):
                os.remove(os.path.join(outdir, file))

    return 0    

def analyzer(rep_dir):
    hits = []
    for file in os.listdir(rep_dir):
        with open(os.path.join(rep_dir, file)) as report:
            report_lines = [line for line in report.readlines() if line.split('\t')[4] != '0']
        
        hits += report_lines

    filename = '_general_report.txt'
    with open(os.path.join(rep_dir, filename), 'w', encoding='utf-8') as final:
        if hits != []:
            for hit in hits:
                final.write(hit)
        else:
            return 1
    
    return 0   
                
def main(siRNA_seq, genome_database, out_dir, reports):
    blast(siRNA_seq, genome_database, out_dir, reports)

    if analyzer(reports) == 0:
        print("hit(s) found in database")
    else:
        print("No hit found in database !")
    return 0

if __name__ == "__main__":
    siRNA_seq, genome_database, out_dir, reports = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    main(siRNA_seq, genome_database, out_dir, reports)

"""
siRNA_seq = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\siRNA_seq"
genome_database = r"D:\iGEM\genome_database"
out_dir = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\out_dir"
reports = r"C:\Subpbiotech_cours\BT4\iGEM\Dry_lab\Génomique\predict_siRNA\reports"

main(siRNA_seq, genome_database, out_dir, reports)
"""