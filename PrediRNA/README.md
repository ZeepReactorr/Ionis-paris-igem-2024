# PrediRNA - User guide

## The problem

RNA interference only works when the small RNA sequence is perfectly stringent toward its target. If even one base amongst the 21 is not paired correctly, the RISC protein complex will not recognize the binding, and the virus genome will not be degraded. Consequently, the two parameters on which depends the long-lasting of our solution for farmers are the virus' sequences we aim, and the rate of mutation of the virus, which are known to be very high. 

## The solution

PrediRNA is a <b>from scratch</b> software aiming to identify the mutational pattern of targetted sequences, resolve it with the right method, and forecast the efficiency of siRNAs targetting a virus over a given period of time (30 years here). This multifunctionnality and versatility makes PrediRNA a complete tool allowing us to target appropriately the BYV sequences, but could also be applied on any other virus ambitionned to be targetted with interfering RNA.

## Getting started

First clone the repository
```
git clone https://gitlab.igem.org/2024/software-tools/ionis-paris/
```

Move into the SafeRNA directory
```
cd ionis-paris/PrediRNA
```

Install the necessary packages
```
pip install -r requirements.txt
```

### Prerequisite

We left all our data on BYV in the repository to use as example. However should you want to use the progam on other viruses, follow the next steps to ensure proper working :
- Compute an MSA of your virus sequences in nucleotide, codons and amino-acid mode.
- Select the most recent sequence and place it in a FASTA file in `alignment_all_seq` named `most_recent.fasta` (needed for seed and sequence prediction, ignore if you just want to run the pattern identification).
- Order you sequences by their time of sequencing and each sequences of the same year in the same FASTA file, but make a folder for each year. Then create a new folder in `alignment_all_seq` and place all your folders in it (needed for the efficiency forecast, ignore if you just want to identify a pattern or predict sequences).

## User guide

### Pattern identification

To run the pattern identification tool:
```
python ~/PATH/TO/seed_entropy.py ~/PATH/TO/MSA.fasta 
```

This program outputs a triple plot figure which you can find in the folder `figures`, from which the pattern of substitution can be identified. Some examples of what the figure should look like are already available on the BYV.

### Sequences prediction

Two steps are required to run the sequence prediction tool. Four arguments are needed for the first step where the sequences' seeds are generated.
```
python ~/PATH/TO/seed_main.py ~/PATH/TO/MSA.fasta protein-annotation newest ~/PATH/TO/most_recent.fasta 
```

Lets break down the arguments: 
- `~/PATH/TO/MSA.fasta`: Path to the MSA (Multiple Sequence Alignment) from which the seeds are generated. <b>Obligatory argument</b>
- `protein-annotation`: name of the studied protein. You can name it however you'd like if it is unknown. <b>Obligatory argument</b>
- `newest`: Determines on which sequence the seeds are calculated. Can either be `newest` or `consensus`. `newest` can only be used if a paht associated to the most recent sequence is added.
- `~/PATH/TO/most_recent.fasta`: Path to the most recent sequence from the MSA.

The seeds can be found in the `temp` folder under .fasta and .csv format for more exhaustive informations regarding the substitution of the seed, its probability, etc.

The second step for sequence predictions is runned as follows:

```
python ~/PATH/TO/seq_prediction.py ~/PATH/TO/MSA.fasta cycles duration rate-sub
```

Let's break down the arguments:
- `~/PATH/TO/MSA.fasta`: Similar to the first step, the path to your MSA.
- `cycles`: number of mutations cycle (we advise to keep under 2 to avoid too long runtime depending on the sample)
- `duration`: how many years would you like to study
- `rate-sub`: The rate of substitution of your sequence

Finally for the efficiency prediction, make sure you followed the instructions for the database creation in the Prerequisite. The program is ran as follow:
```
python ~/PATH/TO/sirna_align.py ~/PATH/TO/db_folder
```

Breakdown of the arguments:
- `~/PATH/TO/db_folder`: the path to the folders where your sequences ordered by year are stored.

The program outpus a double plot figure which can be found in the `figures` folder.

!! Important !! <br>
For the calculation of the worst case scenario, we blace ourselves in the same time-frame as the one where the sequences were sequenced. However the worst case scenario is time-frame independent assuming a constant rate of substitution. Thus for example the results obtained in the time frame 1993-2021 would be identical to those of the time frame 2025-2053.

> The simulation step which was supose to emulate the most recent sequences' evolution was deprecated at the last minute and we are currently working on a fix.
