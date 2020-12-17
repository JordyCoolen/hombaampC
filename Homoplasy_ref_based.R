########                      ########
########  Homoplasy script    ########
########   Reference-based    ########
########  19 September 2019   ########
########                      ########

require(ape)
require(insect)

path = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/"

# set working dir
setwd(path)

core.full.final.aln <- '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/core.full.final.aln'

# read alignment
alignment <- ape::read.dna(core.full.final.aln, format = "fasta")

# extract regions of interest

# read alignment
#alignment_ref <- ape::read.dna("Alignment/Reference.full.aln.fasta", 
#                           format = "fasta")

# mulvey 4470129-4470273 (position may vary if INDELS will be present)
ape::write.dna(alignment[,4470129:4470273], file = "Alignment/mulvey.aln", format = "fasta")

# ampC 4470241-4471374 (position may vary if INDELS will be present)
ape::write.dna(alignment[,4470241:4471374], file = "Alignment/ampC.aln", format = "fasta")

# attenuator 4470198-4470218
ape::write.dna(alignment[,4470198:4470218], file = "Alignment/attenuator.aln", format = "fasta")

# promotor 4470140-4470174
ape::write.dna(alignment[,4470140:4470174], file = "Alignment/promotor.aln", format = "fasta")

# gspL mutation
ape::write.dna(alignment[,833451], file = "gspl_mutation.txt", format="sequential")

library("readxl")

# load overview table
excel = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/20190402_Lijst_AmpC_project_data_promoter_attenuator_correctie_totaal.xlsx"
excel = read_excel(excel)

gspl = read.table("gspl_mutation.txt", sep=" ")
gspl = gspl[2:nrow(gspl),]

combine <- merge(gspl, excel, by.x = "V1", by.y = "unique_ID")
combine <- combine[,c(2,8)]
mel <- melt(table(combine))

# gspL mutation position 833451

#   pampC   not-pampC
# a  33       9+22 (31)      64 total A      
# t  51       10+37 (47)     98 total T

# 162 total strains with this position
# 10 no data on this position

# file paths

promotor <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/promoter/20190430_172Strains_Promoter.fasta"
attenuator <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/attenuater/20190430_172Strains_Attenuator.aln"
chrampC <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/ampC/20190226_All_chrampC_172strains.aln"
coreSNP <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/core.aln"
# coreREF <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/core.full.aln"
coreREF <- "/Users/jordycoolen/Desktop/homoplasy/core.full.final.aln"
coreSNP_tree <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Phylogeny/core.final.tre"

# obtain all chrampC segregating sites and store to disk
ampC_dna <-ape::read.dna(chrampC, format="fasta")
ampC_dna <- ampC_dna[,ape::seg.sites(ampC_dna)]
ape::write.dna(ampC_dna, 
               file = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/ampC/20191003_All_chrampC_172strains_segsites.aln",
               format = "fasta")

########                      ########
########           MSA        ########
########      Visulaization   ########
########                      ########

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("msa")
require(msa)

#TODO
# sort alignment based on order of coreSNP tree

# create sequence logo and msa with msa library
# sequence logo etc
MSA <- readDNAMultipleAlignment(attenuator, format = "fasta")

msaPrettyPrint(MSA, output="tex", 
               paperHeight = 40,
               showNames="right", showLogo="top",
               showLogoScale = "left",
               showNumbering = "none",
               consensusThreshold=50,
               showConsensus="bottom",
               shadingColors = "blues",
               logoColors="rasmol", shadingMode="similar",
               showLegend=TRUE, askForOverwrite=FALSE)

tools::texi2pdf("MSA.tex", clean = TRUE)

###########           ##############
########### TREESPACE ##############
###########           ##############

# Treespace to compare trees
require(treespace)
require(phytools)
require(adegraphics)

if(!file.exists('Phylogeny/core.final.tre')){
  coreSNP.tree <- ape::read.tree("Phylogeny/core.tre")
  coreSNP.tree <- ape::drop.tip(coreSNP.tree, "Reference")
  coreREF.tree <- ape::read.tree("Phylogeny/core.full.tre")
  coreREF.tree <- ape::drop.tip(coreREF.tree, "Reference")
  
  # function to correct the tip names
  correct_name_phy <- function(tree){
    tree$tip.label[tree$tip.label=="ETZ_AWGS170023"] <- "AWGS170023"
    tree$tip.label[tree$tip.label=="ETZ_AWGS170036"] <- "AWGS170036"
    tree$tip.label[tree$tip.label=="B16038396-1"] <- "B16038396-1A"
    tree$tip.label[tree$tip.label=="B16039545-1"] <- "B16039545-1A"
    return(tree)
  }
  
  coreSNP.tree <- correct_name_phy(coreSNP.tree)
  coreREF.tree <- correct_name_phy(coreREF.tree)
  
  # rename tips function
  rename.tips <- function(phy, tiplabels) {
    for (old_names in tiplabels$Old){
      #print(old_names)
      row = match(old_names, tiplabels$Old)
      new_names = tiplabels$ampC_ID[row]
      #print(new_names)
      mpos <- match(old_names, phy$tip.label)
      phy$tip.label[mpos] <- new_names
    }
    return(phy)
  }
  
  # read file to convert tip labels to ampC ID's
  options(stringsAsFactors = FALSE)
  tiplabels = read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/WGS_SNP/Tip_labels.txt',
                         sep="\t", header = 1)
  
  # rename tips to ampC ID's
  coreSNP.tree = rename.tips(coreSNP.tree, tiplabels)
  coreREF.tree = rename.tips(coreREF.tree, tiplabels)
  
  # write new trees
  ape::write.tree(coreSNP.tree, file='Phylogeny/core.final.tre')
  ape::write.tree(coreREF.tree, file='Phylogeny/core.full.final.tre')
} else {
  coreSNP.tree <- ape::read.tree("Phylogeny/core.final.tre")
  coreREF.tree <- ape::read.tree("Phylogeny/core.full.final.tre")
}

# de novo based SNP trees
corekSNP.tree <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/WGS_SNP/core_SNPs_matrix_renamed.tre")
coreRoary.tree <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/Roary_output/core_gene_alignment_renamed.tre")
#denovo.attenuator.tree <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/attenuater/Model_GTR_bootstraps_1000/20190430_172Strains_Attenuator.phy_phyml_tree.txt")
#denovo.promotor.tree <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/promoter/Model_GTR_bootstraps_1000/20190430_172Strains_Promoter.phy_phyml_tree.txt")
#denovo.ampC.tree <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/ampC/Model_GTR_bootstraps_1000/20190226_All_chrampC_172strains.phy_phyml_tree.txt")

denovo.promotor.FT.tre <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/promoter/FastTree_nt_gtr/20190430_172Strains_Promoter.tre")
denovo.attenuator.FT.tre <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/attenuater/FastTree_nt_gtr/20190430_172Strains_Attenuator.tre")
denovo.chrampC.FT.tre <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/ampC/FastTree_nt_gtr/20190226_All_chrampC_172strains.tre")
denovo.chrampC.AA.FT.tre <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/ampC/AA_FastTree_nt_gtr/20190226_All_chrampC_172strains.translated.tre")

#TODO
# rename all tips to ampC coding

Trees <- c(midpoint.root(corekSNP.tree), midpoint.root(coreRoary.tree), 
           #midpoint.root(denovo.attenuator.tree), midpoint.root(denovo.promotor.tree), 
           #midpoint.root(denovo.ampC.tree), 
           midpoint.root(coreSNP.tree), 
           midpoint.root(coreREF.tree),
           midpoint.root(denovo.promotor.FT.tre), midpoint.root(denovo.attenuator.FT.tre),
           midpoint.root(denovo.chrampC.FT.tre), midpoint.root(denovo.chrampC.AA.FT.tre))

names(Trees) <- c("corekSNP","coreRoary",
                  #"attenuator","promotor",
                  #"ampC",
                  "coreSNP",
                  "coreREF",
                  "promotor",
                  "attenuator",
                  "chrampC",
                  "chrampC.AA")

res <- treespace::treespace(Trees, method  = "treeVec", nf = 3)
# table.image
table.image(res$D, nclass=30)
plot(hclust(res$D))

# table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle")

#PCo plot
treespace::plotGroves(res$pco, lab.show=TRUE, lab.cex=1.0)

###########                 ##############
########### HOMOPLASYFINDER ##############
###########                 ##############

# HomoplasyFinder
require(homoplasyFinder)

# create list of alignments
list <- c(chrampC, promotor, attenuator)
names(list) <- c("chrampC", "promotor", "attenuator")

# homoplasy test function on multiple alignment
# aggregates the results in a data.frame
homoplasy_test <- function(list, tree, all_positions=FALSE){
  
  final_table = c()
  
  for (name in names(list)){
    
    print(name)
    print(list[name])
    
    dir = 'out'
    
    dir.create(paste0(getwd(),"/HomoplasyFinder/",dir))
    
    inconsistentPositions <- runHomoplasyFinderInJava(tree, list[name],  
                                                      paste0(getwd(),"/HomoplasyFinder/",dir,"/"),
                                                      includeConsistentSitesInReport=all_positions)
    # read consistentPosition Table
    consis.table <- read.csv2(paste0("HomoplasyFinder/", dir, 
                                     "/consistencyIndexReport_",
                                     format(Sys.time(), "%d-%m-%y"),".txt"), sep="\t")
    
    # add missing positions
    #consis.table <- add_missing_positions(list[name], consis.table)
    
    # change type of column
    consis.table$ConsistencyIndex <- as.numeric(as.character(consis.table$ConsistencyIndex))
    
    # add name column for plots
    consis.table$name <- rep(name, nrow(consis.table))
    
    final_table <- rbind(final_table, consis.table)
    
  }
  
  #final_table$name <- factor(final_table$name, levels = names(list))
  
  return(final_table)
}

markers.path <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/markers_HomoplasyFinder_output.txt"

if(!file.exists(markers.path)){
  markers <- homoplasy_test(list, coreSNP_tree, all_positions = TRUE)
  write.table(markers, file=markers.path,
              sep="\t")
} else {
  markers <- read.table(file=markers.path)
}

markers$MinimumNumberChangesOnTree[markers$MinimumNumberChangesOnTree == "-"] <- 0
markers$MinimumNumberChangesOnTree <- factor(markers$MinimumNumberChangesOnTree)

# transform consistency score to log10 score
markers$ConsistencyIndex <- log10(markers$ConsistencyIndex)

# read in REF_Consistency
REF.CI.path = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/REF_Consistency.txt"

if(!file.exists(REF.CI.path)){
  setwd("/Users/jordycoolen/Desktop/homoplasy/Input")
  
  # function to split alignment into parts and run homoplasyFinder on each part
  split_alignment <- function(alignment, size, coreSNP){
    
    # dividing the sequences in parts to run HomoplasyFinder
    DNAbin <- ape::read.dna(alignment, format ='fasta')
    length <- length(DNAbin[1,])
    
    # divide alignment in parts
    s <- seq(1, length, by = size)
    parts = c()
    
    # create the alignment parts and save to disk
    for (i in seq(1, length(s))){
      start = s[i]
      end = s[i+1] - 1
      if (!is.na(end)){
        print(start)
        print(end)
        dna <- DNAbin[,start:end]
        ape::write.FASTA(dna,
                         file = paste0(i, ".aln"))
      } else {
        print(start)
        print(length)
        dna <- DNAbin[,start:length]
        ape::write.FASTA(dna,
                         file = paste0(i, ".aln"))
      }
    }
    
    # create list of created part of alignment
    file.names <- dir(getwd(), pattern =".aln")
    for(i in 1:length(file.names)){
      print(file.names[i])
    }
    # add names to the values (same as the value)
    names(file.names) <- file.names
    table <- homoplasy_test(file.names, coreSNP, all_positions = TRUE)
    write.table(table, file='Consistency.txt', sep="\t")
    return(table)
  }
  
  # create all parts of alignment and run 
  table <- split_alignment(alignment = coreREF, size = 20000, coreSNP = coreSNP_tree)
  
  ########                      ########
  ########       Homoplasy      ########
  ########                      ########
  
  input = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/REF_Consistency.txt"
  
  setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/")
  
  table <- read.table(input, header = 1)
} else {
  table <- read.table(REF.CI.path, header = 1)
}

# transform consistency score to log10 score
table$ConsistencyIndex <- log10(table$ConsistencyIndex)

# get 1% quantile value (based on Core Roary)
q <- as.data.frame(quantile(table$ConsistencyIndex, probs= c(0.0015)))
colnames(q) <- c('q')

# create plots on ConsistencyIndex scores

require(ggplot2)

# plot density plot with values in regions and quantile values
p <- ggplot(table, aes(x=ConsistencyIndex, mainTitle="ConsistencyIndex")) +
  geom_density(aes(y = ..density..)) +
  geom_point(data=markers, aes(x=ConsistencyIndex, y=name)) +
  geom_vline(data = q, aes(xintercept=q),
             color="blue", linetype="dashed", size=1) +
  scale_x_reverse() +
  theme_minimal()
p

# merge files

table2 <- table
table2$name <- 'genome'
violinplot <- rbind(table2[,1:5], markers)

# Violin plot of ConsistencyIndex per sample
p <- ggplot(violinplot, aes(x=name, y=ConsistencyIndex, addDot=TRUE , dotSize=1 ,
                             mainTitle="ConsistencyIndex")) +
  geom_violin(scale="area", trim=FALSE, fill="gray", color="darkred")+
  #coord_flip() +
  geom_boxplot(width=0.05) +
  #geom_point(alpha=0.3, color="tomato", position = "jitter") +
  theme_classic()
p

# barplot of ConsistencyIndex
ggplot(data=markers[markers$name=="promotor",], aes(x=Position, y=ConsistencyIndex)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_classic()

# barplot of ConsistencyIndex
ggplot(data=markers[markers$name=="chrampC",], aes(x=Position, y=ConsistencyIndex)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_classic()

# barplot of ConsistencyIndex
ggplot(data=markers[markers$name=="attenuator",], aes(x=Position, y=ConsistencyIndex)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_classic()

# not used anymore
markers$MinimumNumberChangesOnTree <- as.numeric(as.character(markers$MinimumNumberChangesOnTree))
#markers$MinimumNumberChangesOnTree <- factor(markers$MinimumNumberChangesOnTree)

# barplot of MinimumNumberChangesOnTree
ggplot(data=markers[markers$name=="promotor",], aes(x=Position, y=MinimumNumberChangesOnTree)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, 40, by = 5)) + 
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 18, by = 2)) +
  expand_limits(y = c(0, 18)) +
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_classic()

# barplot of MinimumNumberChangesOnTree
ggplot(data=markers[markers$name=="attenuator",], aes(x=Position, y=MinimumNumberChangesOnTree)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, 40, by = 5)) + 
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 18, by = 2)) +
  expand_limits(y = c(0, 18)) +
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_classic()

# barplot of MinimumNumberChangesOnTree
ggplot(data=markers[markers$name=="chrampC",], aes(x=Position, y=MinimumNumberChangesOnTree)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  scale_x_continuous(expand = c(0, 0)) + 
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_classic()

########                      ########
########       Homoplasy      ########
########        windows       ########
########                      ########

# avoid scientific notations
options(scipen = 999)

q = 0.125 # 1% quantile
q = 0.05882353 # minimal promotor score (0.07% percentile)
q = 0.07142857 # maximum of two mutations in frame of 11 (0.15% percentile)

# position on promotor with 0.07142857...
q = table[4470150,][2]

# obtain all positions with <= q
CI_q_positions <- table[table$ConsistencyIndex <= as.numeric(q),]

CI_q_positions <- cbind(rep('chr1', nrow(CI_q_positions)), 
      CI_q_positions$Positionnew, 
      CI_q_positions$Positionnew, 
      CI_q_positions$ConsistencyIndex)

write.table(CI_q_positions, file = '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Circos/CI_q_0.07.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# new function test
CI_density <- function(q, table, window, step){
  # count consistency index lower or equal to quartile (q)
  quant_1 <- function(x) {
    return(length(x[x <= as.numeric(q)])/length(x)*100)
  }
  
  # perform sliding window on counts
  res <- evobiR::SlidingWindow(FUN=quant_1,
                               table$ConsistencyIndex,
                               window = window, step = step)
  
  # results to dataframe
  res <- as.data.frame(res)
  res$windownumber <- rownames(res) # window number
  
  # obtain windowrange (start end) of each window
  start <- seq(from = 1, to = (nrow(table) - window), by = step)
  end <- sapply(start, function(x){ x + window})
  windowrange <- cbind(start ,end)
  windowrange <- rbind(windowrange, c(windowrange[dim(windowrange)[1],2], 
                                      nrow(table)))
  
  # calculate counts of last window seperatly
  res <- rbind(res, c(sum(sapply(table$ConsistencyIndex[tail(windowrange, 1)[1]:tail(windowrange, 1)[2]], 
                                 FUN=quant_1)), length(rownames(res))+1))
  res <- cbind(rep('chr1', nrow(res)), windowrange, res)
  colnames(res) <- c('chr', 'start', 'end', 'value', 'options')
  
  # add start end of each window
  
  # l <- cbind(i, dim(res)[1], summary(res$res >= 2)[3])
  # result <- rbind(result, l)
  # print(result)
  return(res)
}

CI_window.table <- CI_density(q, table, 11, 1)

CI_window.table <- CI_window.table[CI_window.table$value >= 10,]

write.table(CI_window.table[,1:4], file = '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Circos/CI_q_0.07_11_1.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# obtain SNP density plot
SNP_density <- function(table, window, step, freq=TRUE){
  # count number of SNPs in window
  SNP.count <- function(x) {
    return(length(x[x != "-"])) # counts if there is a SNP
  }
  
  # count number of SNPs in window
  SNP.count.freq <- function(x) {
    return(length(x[x != "-"])/length(x)*100)
  }
  
  if (freq == TRUE){
    # perform sliding window on counts
    res <- evobiR::SlidingWindow(FUN=SNP.count.freq,
                                 table$MinimumNumberChangesOnTree,
                                 window = window, step = step)
  }
  else{
    # perform sliding window on counts
    res <- evobiR::SlidingWindow(FUN=SNP.count,
                                 table$MinimumNumberChangesOnTree,
                                 window = window, step = step)
  }
  
  # results to dataframe
  res <- as.data.frame(res)
  res$windownumber <- rownames(res) # window number
  
  # obtain windowrange (start end) of each window
  start <- seq(from = 1, to = (nrow(table) - window), by = step)
  end <- sapply(start, function(x){ x + window})
  windowrange <- cbind(start ,end)
  windowrange <- rbind(windowrange, c(windowrange[dim(windowrange)[1],2], 
                                      nrow(table)))
  
  # calculate counts of last window seperatly
  res <- rbind(res, c(sum(sapply(table$MinimumNumberChangesOnTree[tail(windowrange, 1)[1]:tail(windowrange, 1)[2]], 
                                 FUN=SNP.count)), length(rownames(res))+1))
  res <- cbind(rep('chr1', nrow(res)), windowrange, res)
  colnames(res) <- c('chr', 'start', 'end', 'value', 'options')

  return(res)
}

# create Consistency Index count per window
SNPdensity.table <- SNP_density(table, 20, 20)

write.table(SNPdensity.table[,1:4], file = '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Circos/SNPdensity_20_20_freq.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# obtain CI average over window
CI_average <- function(table, window, step){
  
  # perform sliding window on counts
  res <- evobiR::SlidingWindow(FUN=mean,
                               table$ConsistencyIndex,
                               window = window, step = step)
  
  # results to dataframe
  res <- as.data.frame(res)
  res$windownumber <- rownames(res) # window number
  
  # obtain windowrange (start end) of each window
  start <- seq(from = 1, to = (nrow(table) - window), by = step)
  end <- sapply(start, function(x){ x + window})
  windowrange <- cbind(start ,end)
  windowrange <- rbind(windowrange, c(windowrange[dim(windowrange)[1],2], 
                                      nrow(table)))
  
  # calculate counts of last window seperatly
  res <- rbind(res, c(sapply(table$ConsistencyIndex[tail(windowrange, 1)[1]:tail(windowrange, 1)[2]], 
                                 FUN=mean), length(rownames(res))+1))
  res <- cbind(rep('chr1', nrow(res)), windowrange, res)
  colnames(res) <- c('chr', 'start', 'end', 'value', 'options')
  
  return(res)
}

# create Consistency Index count per window
CI_average.table <- CI_average(table, 1000, 1000)

write.table(CI_average.table[,1:4], file = '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Circos/CI_average_1000.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)


## 2020 24 January report Consistency Index on pampC harbouring strains
setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/pampC/")

# HomoplasyFinder
require(homoplasyFinder)

# input is coreSNP tree
coreSNP.tree = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Phylogeny/core.final.tre"
# input is mock alignment sequence file of A non-pampC vs. T pampC (.aln)
pampC = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/pampC/pampC_mock.aln"

inconsistentPositions <- runHomoplasyFinderInJava(coreSNP.tree, pampC,  
                                                  path=paste0(getwd(),"/"),
                                                  includeConsistentSitesInReport=TRUE)

## phangorn method
require(phangorn)
pampC_phyDat = as.phyDat(read.dna(pampC, format = "fasta"))
CI(read.tree(coreSNP.tree), pampC_phyDat, cost = NULL, sitewise = FALSE)



###########                                         ###########
########### 20200311 SNP to methylation correlation ###########
###########                                         ###########

require('evobiR')

SNP.input = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/REF_Consistency.txt"
SNP.table <- read.table(SNP.input, header = 1)

# load SNPdensity function (see above in script)

# create Consistency Index count per window
SNPdensity.table <- SNP_density(SNP.table, 100, 100, freq=FALSE)

# load methylation data

methyl.input = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/methylation/ampC_0069_modifications.gff"
methyl.table <- read.table(methyl.input, skip = 4) # some positions are called multiple times
colnames(methyl.table)<- c("seqid","source","type","start","end","score","strand","phase","attributes")

# TODO
methyl.table <- read.table(text = gsub(";", "\t", readLines(methyl.input)), fill = TRUE) # some positions are called multiple times
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("coverage=", "", x)}))
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("context=", "", x)}))
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("IPDRatio=", "", x)}))
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("frac=", "", x)}))
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("fracLow=", "", x)}))
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("fracUp=", "", x)}))
methyl.table <- data.frame(lapply(methyl.table, function(x) {gsub("identificationQv=", "", x)}))
colnames(methyl.table)<- c("seqid","source","type","start","end","score","strand","phase",
                           "coverage","context","IPDRatio","frac","fracLow","fracUp","identificationQv")

# transform columns to numeric (for calculations)
#methyl.table <- transform(methyl.table, coverage = as.numeric(coverage))
#methyl.table <- transform(methyl.table, context = as.numeric(context))
#methyl.table <- transform(methyl.table, IPDRatio = as.numeric(IPDRatio))
#methyl.table <- transform(methyl.table, frac = as.numeric(frac))
#methyl.table <- transform(methyl.table, fracLow = as.numeric(fracLow))
#methyl.table <- transform(methyl.table, fracUp = as.numeric(fracUp))
#methyl.table <- transform(methyl.table, identificationQv = as.numeric(identificationQv))

summary(methyl.table)

# remove context column
drop <- c("context")
methyl.table = methyl.table[,!(names(methyl.table) %in% drop)]

# TODO set numeric
summary(methyl.table[,9:14])

boxplot(methyl.table[,9:14])

# all genomic positions
positions <- seq(1,5056572)

# filter on type, take position
m6A.table <- methyl.table[methyl.table$type == 'm6A',]
m6A.table$start = as.numeric(as.character(m6A.table$start))
m6A.table$end = as.numeric(as.character(m6A.table$end))
m6A.table$identificationQv = as.numeric(as.character(m6A.table$identificationQv))
#m6A.table <- m6A.table[m6A.table$identificationQv >= 85,]
# detect full and hemi methylated spots
m6A.table$mtype <- sapply(m6A.table$start, function(x) {any(m6A.table$start==x+1 | m6A.table$start==x-1)})
m6A.table$mtype[m6A.table$mtype == TRUE] <- 'full'
m6A.table$mtype[m6A.table$mtype == FALSE] <- 'hemi'
# filter on only full
m6A.table <- m6A.table[m6A.table$mtype == 'hemi',]

# remove mtype column
drop <- c("mtype")
m6A.table = m6A.table[,!(names(m6A.table) %in% drop)]

# add all missing postions and add no methylation
# create list of genomic positions not in methyl data
missing <- setdiff(positions, m6A.table$start)

# create empty matrix
missing.table <- matrix('', nrow = length(missing), ncol = 14)
missing.table <- as.data.frame(missing.table)

# filling missing.table with mock data
colnames(missing.table)<- c("seqid","source","type","start","end","score","strand","phase",
                           "coverage","IPDRatio","frac","fracLow","fracUp","identificationQv")
missing.table$seqid <- 1
missing.table$source <- 'R'
missing.table$type <- '-'
missing.table$start <- missing
missing.table$end <- missing
missing.table$score <- '-'
missing.table$strand <- '-'
missing.table$phase <- '-'
missing.table$coverage <- '-'
missing.table$IPDRatio <- '-'
missing.table$frac <- '-'
missing.table$fracLow <- '-'
missing.table$fracUp <- '-'
missing.table$identificationQv <- '-'

missing.table$start = as.numeric(as.character(missing.table$start))
missing.table$end = as.numeric(as.character(missing.table$end))

missing.table <- transform(missing.table, type = as.factor(type))

# combine missing and methylation position
allmethyl.table <- rbind(missing.table, m6A.table)
rownames(allmethyl.table) <- allmethyl.table$start
allmethyl.table <- allmethyl.table[order(allmethyl.table$start),]

# fix column classes
allmethyl.table$type <- factor(allmethyl.table$type)
allmethyl.table$frac = as.numeric(as.character(allmethyl.table$frac))

# bin methylation data (density data)
methyl_density <- function(table, window, step, freq=TRUE){
  # count number of SNPs in window
  methyl.count <- function(x) {
    return(length(x[x != "-"])) # counts if there is a SNP
  }
  
  # count number of SNPs in window
  methyl.count.freq <- function(x) {
    return(length(x[x != "-"])/length(x)*100)
  }
  
  if (freq == TRUE){
    # perform sliding window on counts
    res <- evobiR::SlidingWindow(FUN=methyl.count.freq,
                                 table$type,
                                 window = window, step = step)
  }
  else{
    # perform sliding window on counts
    res <- evobiR::SlidingWindow(FUN=methyl.count,
                                 table$type,
                                 window = window, step = step)
  }

  # results to dataframe
  res <- as.data.frame(res)
  res$windownumber <- rownames(res) # window number
  
  # obtain windowrange (start end) of each window
  start <- seq(from = 1, to = (nrow(table) - window), by = step)
  end <- sapply(start, function(x){ x + window})
  windowrange <- cbind(start ,end)
  windowrange <- rbind(windowrange, c(windowrange[dim(windowrange)[1],2], 
                                      nrow(table)))
  
  # calculate counts of last window seperatly
  res <- rbind(res, c(sum(sapply(table$type[tail(windowrange, 1)[1]:tail(windowrange, 1)[2]], 
                                 FUN=methyl.count)), length(rownames(res))+1))
  res <- cbind(rep('chr1', nrow(res)), windowrange, res)
  colnames(res) <- c('chr', 'start', 'end', 'value', 'options')
  
  return(res)
}

# create count number of methylation per window
methyldensity.table <- methyl_density(allmethyl.table, 100, 100, freq=FALSE)

# plot SNPdensity vs methyl density
par(mfrow=c(3,1))

require("ggpubr")
require("grid")
require("gridExtra")

t <- textGrob("Results m6A and SNP on 1000 window")

p <- ggscatter(SNPdensity.table, x = "start", y = "value",
          xlab = "genomic position", ylab = "SNP count per 1000")

r <- ggscatter(methyldensity.table, x = "start", y = "value",
          xlab = "genomic position", ylab = "methyl count per 1000")

merged <- cbind(SNPdensity.table$value, methyldensity.table$value)
merged <- as.data.frame(merged)
colnames(merged) <- c("SNP", "methyl")

z <- ggscatter(merged, x = "SNP", y = "methyl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "SNP window 1000", ylab = "m6A window 1000")

grid.arrange(p, t, r, z, ncol=2)

# test extreem SNP locations

extreemSNPs.table <- SNPdensity.table[SNPdensity.table$value >= 80,]
correspondingmethyl.table <- methyldensity.table[rownames(extreemSNPs.table),]

plot(extreemSNPs.table$value, correspondingmethyl.table$value)
cor(extreemSNPs.table$value, correspondingmethyl.table$value)

# plot on genomic position
SNP.table$SNP[SNP.table$MinimumNumberChangesOnTree != '-'] <- 1
SNP.table$full_methyl[m6A.table$mtype != '-'] <- 1

methyl.table$start = as.numeric(as.character(methyl.table$start))
methyl.table$end = as.numeric(as.character(methyl.table$end))
methyl.table$identificationQv = as.numeric(as.character(methyl.table$identificationQv))
methyl.table$mtype <- sapply(methyl.table$start, function(x) {any(methyl.table$start==x+1 | methyl.table$start==x-1)})
methyl.table$mtype[methyl.table$mtype == TRUE] <- 'full'
methyl.table$mtype[methyl.table$mtype == FALSE] <- 'hemi'

total <- merge(SNP.table, methyl.table, by.x="Positionnew", by.y="start", all = TRUE)


# match all kmers of GATC
library(stringr)


# load genome sequence
require(seqinr)
genome <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Reference/assembly_chromosome.fasta"

genome <- read.fasta(file = genome, as.string = TRUE, forceDNAtolower = FALSE)


# detect all gatc positions in forward strand and create bed file
gatc_bed <- str_locate_all(pattern ='GATC', genome$`1`[1])
gatc_bed <- as.data.frame(gatc_bed)
gatc_bed$name <- rownames(gatc_bed)
gatc_bed$seq_id <- 1
gatc_bed$start <- gatc_bed$start-1
gatc_bed <- gatc_bed[, c(4, 1, 2, 3)]

# write to .bed file
write.table(gatc_bed, file='gatc_positions.bed',sep="\t", row.names = F, col.names = F, quote = F)
# modify bed to indicate all A positions on forward and reverse strand
A_positions <- gatc_bed # forward
A_positions$start <- gatc_bed$start + 1
A_positions$end <- A_positions$start
A_rev_positions <- gatc_bed # reverse
A_rev_positions$start <- gatc_bed$start + 2
A_rev_positions$end <- A_rev_positions$start

A_pos_bed <- rbind(A_positions, A_rev_positions)
write.table(A_pos_bed, file='A_positions.bed',sep="\t", row.names = F, col.names = F, quote = F)

# find a way to detect non methylated GATC positions
# find m6A not on GATC positions
# couple does positions to mutation locations

## April 23 2020
## loading gubbins gff file
## parse for use in circos

path = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins"

# set working dir
setwd(path)

gubbins <- read.table(file="ampC.recombination_predictions_min100.gff", sep="\t")
gubbins2 <- read.table(file="ampC.recombination_predictions_min100.gff", sep=";")

library(stringr)

gubbins$count <- sapply(gubbins2$V3, function(x){ str_count(x, pattern="ampC")})
gubbins$chr <- rep('chr1', nrow(gubbins))

gubbins <- gubbins[c(11, 4,5,10)]

write.table(gubbins, file="gubbins_circos.txt", 
            col.names=FALSE, row.names=FALSE, quote=FALSE,
            sep="\t")

plot(gubbins$V4 + gubbins$V5)

newdata <- gubbins[order(gubbins$V4, gubbins$V5),]

# HomoplasyFinder on correct recombination tree from gubbins
path = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins"

# set working dir
setwd(path)
require(homoplasyFinder)
gubbins.tree <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins/ampC.final_tree.tre"
gubbins.tree <- ape::read.tree("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins/ampC.final_tree.tre")

core.tree <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Phylogeny/core.full.final_midpoint.tre"

runHomoplasyFinderInJava(gubbins.tree, promotor, path=getwd(), includeConsistentSitesInReport = TRUE)

####                  ####
#### Rebuttal         ####
#### 29 October 2020  ####
####                  ####

# create venndiagrom of pyseer and our method

# load pyseer data
# position pvalue

# load pyseer data removed recombination with gubbins
# position pvalue

# load our fisher exact method
# position pvalue

# Load library
library(VennDiagram)

# Dataset
# our method
set1 <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_FDR.txt")

list1 <- as.numeric(set1$V2)
#set1 <- set1[1:100]

set1 <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_final_positions_circos.txt")
list1 <- as.numeric(set1$V2)

#pyseer
require(reshape2)
set2 <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/pyseer.assoc.p1",
                   header = 1, fill = TRUE)
split <- set2$variant
split <- colsplit(split, "_", c("reference", "position","REF","ALT"))
set2 <- cbind(split , set2)
set2 <- set2[set2$lrt.pvalue<=0.05,]
#set2 <- set2[set2$REF == 'A' | set2$REF == 'T' | set2$REF == 'C' | set2$REF == 'G',]
#set2 <- set2[set2$ALT == 'A' | set2$ALT == 'T' | set2$ALT == 'C' | set2$ALT == 'G',]
#set2 <- paste(set2$position)
list2 <- as.numeric(set2$position)
list2 <- unique(list2)

#pyseer gubbins
set3 <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/pyseer.assoc.gubbins.corrected.p1",
                   header = 1, fill = TRUE)
split <- set3$variant
split <- colsplit(split, "_", c("reference", "position","REF","ALT"))
set3 <- cbind(split , set3)
set3 <- set3[set3$lrt.pvalue<=0.05,]
#set3 <- set3[set3$REF == 'A' | set3$REF == 'T' | set3$REF == 'C' | set3$REF == 'G',]
#set3 <- set3[set3$ALT == 'A' | set3$ALT == 'T' | set3$ALT == 'C' | set3$ALT == 'G',]
#set3 <- paste(set3$position)
list3 <- as.numeric(set3$position)
list3 <- unique(list3)

setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(list1, list2, list3),
  category.names = c("Our" , "pyseer" , "pyseer-rec"),
  filename = '#14_venn_diagramm_lrtpvalue.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

intersect(list1, list2)
intersect(list1, list3)

setdiff(set2, set1)

#####                                                     ##### 
##### For pyseer tre names need to be the original names  #####  
##### renaming tip labels back to original labels         ##### 
#####                                                     #####

# rename tips function
rename.tips <- function(phy, tiplabels) {
  for (old_names in tiplabels$Old){
    #print(old_names)
    row = match(old_names, tiplabels$Old)
    new_names = tiplabels$New[row]
    #print(new_names)
    mpos <- match(old_names, phy$tip.label)
    phy$tip.label[mpos] <- new_names
  }
  return(phy)
}

# full tree (with ampC names)
full.tree <- ape::read.tree('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/core.full.final.tre')

# read file to convert tip labels to ampC ID's
options(stringsAsFactors = FALSE)
tiplabels = read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/Tip_labels_reverse.txt',
                       sep="\t", header = 1)

# rename tips to original ID's
full.tree = rename.tips(full.tree, tiplabels)

# write new trees
ape::write.tree(full.tree, file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/core.full.final.oldnames.tre')

##### also rename gubbins corrected tree to original ID's
gub.tree <- ape::read.tree('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins/ampC.final_tree.tre')

# rename tips to original ID's
gub.tree = rename.tips(gub.tree, tiplabels)

# write new trees
ape::write.tree(gub.tree, file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins/ampC.final_tree.oldnames.tre')

#####                                       ##### 
##### Matching the 24 found positions       ##### 
##### and compare them to the pyseer output ##### 
#####                                       ##### 

require(reshape2)
# pyseer gubbins corrected results
pyseer <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/pyseer.assoc.gubbins.corrected.p1",
                   header = 1, fill = TRUE)
# modify the matrix
variant <- pyseer$variant
variant <- colsplit(variant, "_", c("reference", "position","REF","ALT"))
pyseer <- cbind(variant, pyseer)

# filter on ltr-pvalue <= 0.05
filter <- pyseer[pyseer$filter.pvalue<=0.05,]
# count unique positions
length(unique(filter$position))

filter <- filter[filter$lrt.pvalue<=0.05,]
# count unique positions
length(unique(filter$position))

# 24 positions
set1 <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_final_positions_circos.txt")
list1 <- as.numeric(set1$V2)

# create new matrix with only the SNPs (removing complex en INDELS)
# pyseer gubbins corrected results
pyseer <- read.table("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/pyseer/pyseer.assoc.gubbins.corrected.p1",
                     header = 1, fill = TRUE)
# modify the matrix
variant <- pyseer$variant
variant <- colsplit(variant, "_", c("reference", "position","REF","ALT"))
pyseer <- cbind(variant, pyseer)
new <- pyseer[pyseer$REF == 'A' | pyseer$REF == 'T' | pyseer$REF == 'C' | pyseer$REF == 'G',]

new <- test[row.names(unique(test$position)),]

# create empty data.frame
output <- pyseer[pyseer$position == 0,1:8]

# get pyseer results of 24 positions and save in output variable
for (i in list1){
  output <- rbind(output, pyseer[pyseer$position == i,1:8])
}

# position 843887
pyseer[pyseer$position>1862994 & pyseer$position<1863014, 1:8]

pyseer[pyseer$position>1862980 & pyseer$position<1863018, 1:8]

#2057518
pyseer[pyseer$position>2057408 & pyseer$position<2057538, 1:8]

2057518

#810791
pyseer[pyseer$position>810781 & pyseer$position<810801, 1:8]

#824522
pyseer[pyseer$position>824512 & pyseer$position<824532, 1:8]

#830684
pyseer[pyseer$position>830674 & pyseer$position<830694, 1:8]

#830708
pyseer[pyseer$position>830698 & pyseer$position<830718, 1:8]

#832152
pyseer[pyseer$position>832142 & pyseer$position<832162, 1:8]

#832287
pyseer[pyseer$position>832247 & pyseer$position<832317, 1:8]

#833451
pyseer[pyseer$position>833401 & pyseer$position<833521, 1:8]

#1952745
pyseer[pyseer$position>1952735 & pyseer$position<1952755, 1:8]

#2051214
pyseer[pyseer$position>2051204 & pyseer$position<2051224, 1:8]

#####                           ##### 
##### HomoplasyFinder           ##### 
##### On INDEL positions        ##### 
#####  promotor and attenuator  ##### 
#####                           ##### 

# HomoplasyFinder
require(homoplasyFinder)

# inputs
coreSNP_tree <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Phylogeny/core.final.tre"
INDELS <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/INDEL_promotor_attenuator/presence-absence_INDELS.csv"

setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/")

dir.create(paste0(getwd(),"/homoplasy_INDELS"))

inconsistentPositions <- runHomoplasyFinderInJava(coreSNP_tree, presenceAbsenceFile = INDELS,  
                                                  path = paste0(getwd(),"/homoplasy_INDELS","/"),
                                                  includeConsistentSitesInReport = TRUE)