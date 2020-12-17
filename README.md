# hombaampC
Homoplasy-based analysis on FOX resistant ESBL negative Escherichia coli.

# note
These scripts are the actual scripts used in the studie and are not polished.

Homoplasy_ref_based.R
  - extract regions
  - MSA visualization
  - Treespace (not used)
  - HomoplasyFinder
    - Homoplasy
  - Homoplasy windows
  - SNP to methylation correlation (not used)
  - pyseer renaming and parsing
  - compare results with pyseer
  - homoplasy INDEL positions promotor and attenuator


Calculate_Fisher_nt_trait.R
  - mutation to trait
  - CIRCOS plot
  - INDEL promotor and attenuator and fisher exact
  - 24 contingency table/ Fisher exact
  - gubbins VCF file parse
  - Table 1 piecharts
  - combine all files (consistency index + fisher + FDR)

Folder:
Alternative_analysis
Contains Rscripts only used during preliminary analysis and alternative analysis stategies.
Code not used for obtaining the final results.

# dependencies
R packages:
ape
insect
readxl
msa
homoplasyFinder
ggplot2
phangorn
stringr
VennDiagram
reshape2
RColorBrewer
scales

R packages optional:
treespace
phytools
adegraphics
evobiR
ggpubr
grid
gridExtra
seqinr


# best practise
See how the steps are performed and reproduce and reuse for your own study.

I would recommend to use a conda environment and new environment with
R and Rstudio to perform the analysis.

# preprint (version 2,  July 07, 2020)
Genome-wide analysis in Escherichia coli unravels an unprecedented level of genetic homoplasy associated with cefotaxime resistance

Jordy P.M. Coolen, Evert P.M. den Drijver, Jaco J. Verweij, Jodie A. Schildkraut, Kornelia Neveling, Willem J.G. Melchers, 
Eva Kolwijck, Heiman F.L. Wertheim, Jan A.J.W. Kluytmans,  Martijn A. Huynen

https://doi.org/10.1101/2020.06.01.128843

# authors
J.P.M. Coolen and E.P.M. den Drijver

# license
[MIT](https://choosealicense.com/licenses/mit/)
