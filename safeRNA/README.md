# safeRNA role and user guide 

## The problem

The siRNAs produced are short, averaging between 21 and 24 nucleotides. They are highly specific, but their small size could lead to the targeting of undesired sequences present in transcriptomes by chance. The tool therefore aims to check that the sequences produced will target only the virus and no others.

## Solution principle

To address this issue and preserve the soil's diversity and especially the long neglected microbiome, we collected every available genomes about the taxa found in the French soil microbiome associated with *Beta vulgaris* culture. Then all those genomes would be aligned with the interfering RNA sequence planned to be engineered to identify which can be kept and which have to go. Only the coding sequences from organisms are retrieved as siRNA only target RNA. To limit the risk of sequencing errors interefering, multiple sequencing from the same strain/species is allowed.

## The programs

The tool combines three main programs, all written in Python 3.11.7 : 
- main.py : from precursor sequences, calculate all the possible siRNA sequences that can be obtained using a sliding window.
- get_data.py : allow to download all coding sequences from a given taxon in NCBI
- align.py : makes the alignment between siRNA sequences and every genome available in the genome database gained using get_data.py and generates the report

### main.py

Written in python, works with terminal command lines. 2 arguments :
- `dir_in` : The directory with the sequences fasta files of precursor RNA sequences 
- `dir_out` : The directory in which all the possible siRNA fasta files will be generated

```BASH
python main.py ~/DIR/IN ~/DIR/OUT
```

### get_data.py

Written in python, works with terminal command lines and a side file containing the list of taxon's name for which genomes have to be downloaded. Two arguments are required : 
- the PATH to the .txt file containing the list of taxon's name.
- The PATH to the directory where you want the database to be downloaded

```BASH
python get_data.py ~/PATH/TO/data_to_get.txt ~/PATH/TO/DATABASE/DIR
```

A report summarizing the number of genomes dowloaded by taxa is also generated.

### align.py

Written in python and is the core of the program. Takes four arguments in the command line :
- siRNA_seq : The directory where siRNA sequences fasta file can be found
- genome_database : the directory where the NCBI genomes have been downloaded
- out_dir : the directory in which the program will work (**keep it empty**)
- reports : the directory in which reports will be generated per taxa and **where the final report named _general_report.txt will be generated**

**The file _general_report.txt contains only information about the HITS that were detected**

```BASH
python align.py ~/PATH/TO/SIRNA/DIR ~/PATH/TO/GENOME/DOWNLOADED ~/PATH/TO/OUT_DIR ~/PATH/TO/REPORT/DIR
```

**Format of the file _general_report.txt** <br>

The report contains 8 different columns

Column name | Column meaning
--- | ---
precursor vs target | Give the name of the precursor from which siRNA are obtained versus the target
precursor | name of the precursor protein 
target | ID of the genome targeted
taxa | the taxon database in which the target was classified when downloaded
number of hits | number of siRNA sequences with a hit in the targeted genome
list of siRNA ID | the ID of the siRNA sequence that hit the targeted genome
Path to genome file | The path to the genome file in the database locally built from NCBI
link | URL directing toward the Genome Assembly Description

## References

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., and Madden, T.L. 2009. BLAST+: architecture and applications. BMC Bioinformatics, 10, 421.