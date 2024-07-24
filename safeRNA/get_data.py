import sys, os, subprocess, re, zipfile

DATA_TO_GET_DIR = r"D:\iGEM\tool"
DIR = r"D:\iGEM\genome_database"

def report(taxon, status, info):
    if status == "successful":
        n_genomes = re.findall("[0-9]+", str(info).split("genome")[0])
    else :
        n_genomes = "NaN"

    with open("report.txt", 'a', encoding='utf-8') as report:
        report.write(f"{taxon}\t{status}\t{n_genomes[0]}\n")

def getter(taxon):
    os.chdir(DIR)
    test = subprocess.run(
            f"datasets download genome taxon {taxon} --assembly-level complete --assembly-version all \
            --exclude-atypical --filename {taxon}_db.zip --include cds --released-after 2010",
            capture_output=True, shell=True
        )

    if "Error" in str(test.stderr) :
        report(taxon, "failed", test.stderr)    

    else:
        report(taxon, "successful", test.stderr)

def main(taxon_file):
    with open(taxon_file, 'r', encoding='utf-8') as taxon_list:
        taxon_list = ''.join(taxon_list.readlines()).split('\n')
    for taxon in taxon_list:
        getter(taxon)

        try : 
            os.makedirs(os.path.join(DIR, f"{taxon}_db"))
            with zipfile.ZipFile(f"{taxon}_db.zip", 'r') as zip_ref:
                zip_ref.extractall(os.path.join(DIR, f"{taxon}_db"))     
        except Exception :
            continue

        os.remove(f"{taxon}_db.zip")
        print(f"{taxon} done")
    
if __name__ == '__main__':
    if "report.txt" in os.listdir():
        os.remove(os.path.join(DIR, "report.txt"))

    taxon_file = sys.argv[1]
    main(taxon_file=taxon_file)

