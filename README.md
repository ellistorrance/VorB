# VorB
"Virus or Bacteria?" (VorB) is a tool for population level surveillance of putative prophages & virus-like bacterial machinery. 

VorB compares:(1) a user-submitted GeNomad output (https://github.com/apcamargo/genomad) with (2) a user-submitted population of bacterial or archaea annotated genomes (.fna) to help the user differentiate between true prophages and host machinery with phage homology.

VorB is written with Python v3.7 and requires Python libraries: NumPy, Matplotlib, Seaborn, and Pandas as well as the bioinformatic tools Diamond (https://github.com/bbuchfink/diamond), MAFFT (https://mafft.cbrc.jp/alignment/software/), and HHsuite (https://github.com/soedinglab/hh-suite).

VorB Requires:
1) a path to a GeNomad output folder.
2) a path to a folder containing annotated genomes in fasta format (.fna) from related taxa. NOTE: make sure no other files are in this folder other than the annotated genomes you wish to query. 

To run:
          
          python start.py -in “path to genomad output” -query “path to folder containing annotated host seqs”

Optional Flags:
          
          -out : specify where VorB will output to. Default is within the genomad output folder.
          -restart : restarts run, deletes any files made by VorB in previous runs.

Test environment is correct with included files:
          
          python start.py -in <path>/Staphylococcus_hominis/genomad/ -query <path>/Staphylococcus_hominis/genes/ -out <path>/Staphylococcus_hominis/
          
          This test should run in ~10 minutes and output to S. hominis folder.


VorB File Output:

![Screenshot 2023-06-21 at 12 59 02 PM](https://github.com/ellistorrance/VorB/assets/60077187/ddef5ed4-cd5a-464d-bc7b-c8a6b5db8d69)



VorB Graphical Output:
![Screenshot 2023-06-21 at 12 59 13 PM](https://github.com/ellistorrance/VorB/assets/60077187/0b1b0c14-1867-4f66-a1c3-c5e74135e807)

^ The example above shows a bacteria with 2 putative prophages as annotated by GeNomad (V1 (green) & V2 (orange)). Host (bacterial/archaeal) universal COGs are included (in blue) to act as reference in the case the user chooses to use more distantly related taxon in their query group. Here, the scatterplot on the left compares the number of strains found in the query group (total N=22) containing a gene with homology to the putative prophage genes identified by GeNomad versus the % average pairwise identity (nucleotide) of the gene across query strains it was identified in.  We see V1 genes are almost universaly present in all strains queried at high average pairwise identity. On the right is barchart of the PHROG annotatation profile of the genes from the putative prophages identified with genomad (y) and the number of genes from each putative propahge which fell within each category. The high level of conservation in V1 evidenced by the graph on the left as well as its lack of head and packaging genes indicate the V1 is likely a vertically transmitted host gene set with phage homology. In this case, a Tailocin. 

Citations:

You can move, but you can't hide: identification of mobile genetic elements with geNomad
Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. — bioRxiv (2023), DOI: 10.1101/2020.11.01.361691.

HH-suite3 for fast remote homology detection and deep protein annotation
Steinegger M, Meier M, Mirdita M, Vöhringer H, Haunsberger S J, and Söding J, BMC Bioinformatics (2019), 473. doi: 10.1186/s12859-019-3019-7

PHROG : families of prokaryotic virus proteins clustered using remote homology.
Terzian P*, Olo Ndela E*, Galiez C, Lossouarn J, Pérez Bucio RE, Mom R, Toussaint A, Petit MA, Enault F.
NAR Genomics and Bioinformatics, Volume 3, Issue 3, September 2021, lqab067, https://doi.org/10.1093/nargab/lqab067

Sensitive protein alignments at tree-of-life scale using DIAMOND
Buchfink B, Reuter K, Drost HG, Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

MAFFT multiple sequence alignment software version 7: improvements in performance and usability.
Katoh, Standley. Molecular Biology and Evolution (2013) 30:772-780



