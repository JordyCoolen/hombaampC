library(phylotools)

# rename fasta files to ampC names

names <- read.delim("/Users/jordycoolen/Downloads/20190621_ampC_genes_Voorbereiding_Checkpoints/ampC_200nt_and_ampCgenes/20190621_ampC_genes_naming.txt")

rename.fasta("/Users/jordycoolen/Downloads/20190621_ampC_genes_Voorbereiding_Checkpoints/ampC_200nt_and_ampCgenes/ampC_all.aln",
             names,
             "/Users/jordycoolen/Downloads/20190621_ampC_genes_Voorbereiding_Checkpoints/ampC_200nt_and_ampCgenes/test.aln")
