# Ionis-Paris-iGEM2024

## Project description 

A quarter of the worlds’ sugar is obtained from sugar beets, which are cultivated all around the world. It is especially important in France, and one of our biggest agricultural sector, from which is obtained the sugar we pour in our tea, coffee, or our cars, under the form of biofuel.

The aphids Myzus persicae and Aphis fabae are the most threatening pests toward B. vulgaris. They can resist most pesticides, and are the main vectors for the Beet Yellow Virus (BYV) [2]. More than a third of the crops are affected yearly, and infected plants suffer a yield loss that can reach 50%, which is why agricultors are wary of these pests as they represent a major threat toward the crops [3].

Neonicotinoids (NEOs) prevail as the sole pesticides able to repel efficiently enough the aphids to protect the crops [3]. As a result, in 2016, about 100% of sugar beet seeds were undergoing preventive treatment of NEOs [4]. However, in 2018 it was confirmed that this pesticide had adverse effects on human and pollinators health. However there usage was not totally restricted in the EU until 2023, as years of derogation were granted especially for sugar beets [5]. Now deprived of their most efficient defense, european beet and beet farmers are most vulnerable to pests.

Our iGEM project was born to adress this issue and save the Sugar Beets ! This project is carried out by 12 students of SupBiotech in France, keen on challenges, adventure and agroecology !

Our idea is to engineer precursors of small interfering RNAs (siRNA) that naturally exists in in plants as their main defense mechanism, and to target and destroy specifically viral RNA[6]. This biomimetic approach will be combined with a second one, hijacking the properties of viral capisdes and use them as vector for the RNAs to penetrate the plant.

## Software description

Cap'siRNA's dry lab developped three softwares to support the project. Programming language is Python 3.11.7 for all three programs, and the Kernel for the jupyter notebooks is Python 3.11.7 as well. Aphidisperse and PrediRNA were fully developped by our Dry lab and are not based on any pre-existing programs. SafeRNA leverages NCBI's tools, datasets.exe and local BLAST to do its office. Aside from that, SafeRNA developpment was also done from scratch.

Documentation and last code version will be available on October 2nd.
- Aphidisperse : Aphid dispersion model <br>
- PrediRNA: Prediction of BYV mutation and effect on siRNA efficiency. <br>
- SafeRNA : Tool leveraging NCBI's datasets tool to ensure specificity of siRNAs toward BYV.

## References

[1] Hossain R, Menzel W, Lachmann C, Varrelmann M. New insights into virus yellows distribution in Europe and effects of beet yellows virus, beet mild yellowing virus, and beet chlorosis virus on sugar beet yield following field inoculation. Plant Pathol. 2021;70(3):584-593. doi:10.1111/ppa.13306

[2] Margaritopoulos JT, Kati AN, Voudouris CC, Skouras PJ, Tsitsipis JA. Long-term studies on the evolution of resistance of Myzus persicae (Hemiptera: Aphididae) to insecticides in Greece. Bull Entomol Res. 2021;111(1):1-16. doi:10.1017/S0007485320000334

[3] Grimmer MK, Bean KMR, Qi A, Stevens M, Asher MJC. The action of three Beet yellows virus resistance QTLs depends on alleles at a novel genetic locus that controls symptom development. Plant Breed. 2008;127(4):391-397. doi:10.1111/j.1439-0523.2008.01515.x

[4] Hauer M, Hansen AL, Manderyck B, et al. Neonicotinoids in sugar beet cultivation in Central and Northern Europe: Efficacy and environmental impact of neonicotinoid seed treatments and alternative measures. Crop Prot. 2017;93:132-142. doi:10.1016/j.cropro.2016.11.034

[5] Zhang D, Lu S. Human exposure to neonicotinoids and the associated health risks: A review. Environ Int. 2022;163(March):107201. doi:10.1016/j.envint.2022.107201

[6] Zhan J, Meyers BC. Plant Small RNAs: Their Biogenesis, Regulatory Roles, and Functions. Annu Rev Plant Biol. 2023;74:21-51. doi:10.1146/annurev-arplant-070122-035226
