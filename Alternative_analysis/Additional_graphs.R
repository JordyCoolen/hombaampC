## read matrix and create heatmap

m <- read.csv('Core_genes/Matrix.csv', row.names = 1)
heatmap(as.matrix(m))
plot(hclust(as.dist(m), method="average"), main='UPGMA', cex=0.5)
## create table of results

require(openxlsx)
library("dplyr")
library(circlepackeR) #devtools::install_github("jeromefroe/circlepackeR")
library(data.tree) #install.packages("data.tree")
setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/")
data <- read.xlsx('20190606_Lijst_AmpC_project_data_promoter_attenuator_correctie_mulvey_typen.xlsx', sheet=1)
data <- cbind(data, 'species'=rep("Ecoli",172))
# subdata <- select(data, "Sets.(1:Radboudumc,2:External,3:Amphia)", 
#                   MLST, promoter.variant, attenuator.variant)
subdata <- select(data, MLST, AmpC.gen.nucleotide.variant)
subdata <- select(data, MLST, AmpC.gen.amino.acid.variant)
subdata$MLST[subdata$MLST == "-"] <- "Unknown"
subdata <- plyr::count(subdata)
# subdata <- subdata[subdata$freq != 1,]

subdata$pathString <- paste("Ecoli", subdata$MLST, subdata$AmpC.gen.nucleotide.variant, sep = "/")
subdata$pathString <- paste("Ecoli", subdata$AmpC.gen.nucleotide.variant, subdata$MLST, sep = "/")
subdata$pathString <- paste("Ecoli", subdata$AmpC.gen.amino.acid.variant, subdata$MLST, sep = "/")
subdata$pathString <- paste("Ecoli", subdata$MLST, subdata$AmpC.gen.amino.acid.variant, sep = "/")

# subdata$pathString <- paste("Ecoli", subdata$Sets..1.Radboudumc.2.External.3.Amphia., 
#                             subdata$MLST, subdata$promoter.variant, subdata$attenuator.variant, sep = "/")
#subdata$pathString <- paste("Ecoli", subdata$attenuator.variant, subdata$MLST, sep = "/")
population <- as.Node(subdata)

# Make the plot
circlepackeR(population, size = "freq")

# network to, from to freq

as.igraph(population, directed=TRUE)

edges <- ToDataFrameNetwork(population, "freq")

vertices <- ToDataFrameTypeCol(population, "freq")

# Then we have to make a 'graph' object using the igraph library:
mygraph <- graph_from_data_frame( edges)

# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
edges=flare$edges
vertices = flare$vertices
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Make the plot
ggraph(mygraph, layout = 'circlepack') + 
  geom_node_circle() +
  theme_void()

# Left: color depends of depth
p=ggraph(mygraph, layout = 'circlepack', weight="freq") + 
  geom_node_circle(aes(fill = depth)) +
  theme_void() + 
  theme(legend.position="FALSE")
p

# Adjust color palette:
p + scale_fill_viridis()
p + scale_fill_distiller(palette = "RdPu") 

p=ggraph(mygraph, layout = 'circlepack') + 
  geom_node_circle(aes(fill = depth)) +
  geom_node_text( aes(label=name, filter=leaf, fill=depth)) +
  theme_void() + 
  theme(legend.position="FALSE") + 
  scale_fill_viridis()
p

# Libraries
library(ggraph)
library(igraph)
library(tidyverse)

# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
edges=flare$edges

# Usually we associate another dataset that give information about each node of the dataset:
vertices = flare$vertices

# Then we have to make a 'graph' object using the igraph library:
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Make the plot
ggraph(mygraph, layout = 'circlepack') + 
  geom_node_circle() +
  theme_void()

## test for comparing sequences (core sequences)

string.diff.ex<-function(a="ATTCGAN",b="attTGTT",exclude=c("n","N","?"),ignore.case=TRUE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case==TRUE)
  {
    a<-toupper(a)
    b<-toupper(b)
  }
  diff.a<-unlist(strsplit(a,split=""))
  diff.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(diff.a,diff.b)
  for(ex.loop in 1:length(exclude))
  {
    diff.d<-diff.d[,!(diff.d[1,]==exclude[ex.loop]|diff.d[2,]==exclude[ex.loop])]
  }
  differences<-sum(diff.d[1,]!=diff.d[2,])
  return(differences)
}

string.diff.ex(a="ACTGACAGG", "ACTTACAGT")

library(seqinr)
alignm <- read.alignment('WGS_SNP/core_SNPs_matrix.fasta', format = 'fasta')
dist_alignm <- dist.alignment(alignm, matrix="identity")
summary(dist_alignm)
plot(hclust(as.dist(dist_alignm), method="average"), main='UPGMA', cex=0.5)

library(ape)
library(phytools)
coreTree <- read.tree('WGS_SNP/core_SNPs_matrix.tre')

coreTree <- midpoint.root(coreTree)
dendrogram <- chronos(coreTree)

cutree(coreTree)

library(tidytree) # install.packages("tidytree")
library(ggtree) # install.packages("ggtree")
library(treeio) # install.packages("treeio")

read.phylip.tree('WGS_SNP/core_SNPs_matrix.tre')




