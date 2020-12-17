######    2019 23 October     ######
######    JORDY COOLEN        ######
######    Mutation to trait   ######

# functions

# function to split alignment into parts and run homoplasyFinder on each part
split_alignment <- function(alignment, size){
  
  # dividing the core_Roary sequences to run HomoplasyFinder
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
      if (i < 10){
        n = paste0('00',i)
      } else if (i >= 10 & i < 100){
        n = paste0('0',i)
      }
      else{
        n = i
      }
      ape::write.FASTA(dna,
                       file = paste0(n, ".aln"))
    } else {
      print(start)
      print(length)
      dna <- DNAbin[,start:length]
      if (i < 10){
        n = paste0('00',i)
      } else if (i >= 10 & i < 100){
        n = paste0('0',i)
      } else{
        n = i
      }
      ape::write.FASTA(dna,
                       file = paste0(n, ".aln"))
    }
  }
  
  # create list of created part of alignment
  file.names <- dir(getwd(), pattern =".aln", full.names = TRUE)
  for(i in 1:length(file.names)){
    print(file.names[i])
  }
  # add names to the values (same as the value)
  names(file.names) <- base::basename(file.names)
  return(file.names)
}

setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/segments")

# read alignment file
alignment <- ape::read.dna('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/promotor.aln',
                           format = 'fasta', as.character = TRUE)

alignment <- ape::read.dna('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/ampC.aln',
                           format = 'fasta', as.character = TRUE)

alignment <- ape::read.dna('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/core.aln',
                           format = 'fasta', as.character = TRUE)

# perform on the full alignment (core.full.final.aln)
# this will split the full alignment in smaller segments
segments = split_alignment('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/core.full.final.aln',
                size = 10000)

# alternative if segments are already made
segments <- dir(getwd(), pattern =".aln", full.names = TRUE)

# test on a small set of segments
#segments_short <- segments[1:4]

phenpath <- '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191028_phenotype_CTX_corrected.txt'

# this function will loop over all segments in file list and perform
# fisher exact test on traits and FDR
fisher_on_genome <- function(phenpath, segments){
  
  # read phenotype data
  # read CTX R and S for non-pampC
  phenotable <- read.table(phenpath, 
                     sep = "\t", header = F)
  
  for (segment in segments){
    alignment <- ape::read.dna(segment, format = 'fasta', as.character = TRUE)
    segment <- base::basename(segment)
    
    # remove samples of which no phenotypic data is provided
    alignment <- alignment[rownames(alignment) %in% phenotable$V1,]
    
    # sort data so the alignment matches with the phenotype data
    alignment <- alignment[match(phenotable$V1, rownames(alignment)),]
    
    # check if order is correct before proceding
    rownames(alignment) == phenotable$V1
    
    # obtain only R and S
    phen <- phenotable$V2
    
    # results
    tbl <- list()
    res <- list()
    
    # function to count nucleotide
    # detect_nucleotide <- function(x, nucleotide){
    #   if (x == nucleotide){
    #     return(1)
    #   } else{
    #     return(0)
    #   }
    # }
    
    # create A matrix vs non
    
    #snps <- apply(alignment, MARGIN=c(1,2), function(x) detect_nucleotide(x, 'a'))
    snps <- alignment
    snps[snps == 'a'] <- 1
    snps[snps != 1] <- 0
    colnames(snps) <- seq(1, ncol(snps))
    
    tbl$A_snps <- snps
    
    # create T matrix vs non
    
    #snps <- apply(alignment, MARGIN=c(1,2), function(x) detect_nucleotide(x, 't'))
    
    snps <- alignment
    snps[snps == 't'] <- 1
    snps[snps != 1] <- 0
    colnames(snps) <- seq(1, ncol(snps))
    
    tbl$T_snps <- snps
    
    # create C matrix vs non
    
    #snps <- apply(alignment, MARGIN=c(1,2), function(x) detect_nucleotide(x, 'c'))
    
    snps <- alignment
    snps[snps == 'c'] <- 1
    snps[snps != 1] <- 0
    colnames(snps) <- seq(1, ncol(snps))
    
    tbl$C_snps <- snps
    
    # create G matrix vs non
    
    #snps <- apply(alignment, MARGIN=c(1,2), function(x) detect_nucleotide(x, 'g'))
    
    snps <- alignment
    snps[snps == 'g'] <- 1
    snps[snps != 1] <- 0
    colnames(snps) <- seq(1, ncol(snps))
    
    tbl$G_snps <- snps
    
    # calculate Fisher
    
    # example
    # stats::fisher.test(table(factor(tbl$T_snps[,1], levels=c(0,1)), phen))
    
    # fisher for A
    pval <- apply(tbl$A_snps, 2, function(e)
      stats::fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
    
    res$A_fisher <- pval
    
    #pval.corrected.fdr <- p.adjust(pval, method="fdr")
    
    #res$A_fdr <- pval.corrected.fdr
    
    # fisher for T
    pval <- apply(tbl$T_snps, 2, function(e)
      stats::fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
    
    res$T_fisher <- pval
    
    #pval.corrected.fdr <- p.adjust(pval, method="fdr")
    
    #res$T_fdr <- pval.corrected.fdr
    
    # fisher for C
    pval <- apply(tbl$C_snps, 2, function(e)
      stats::fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
    
    res$C_fisher <- pval
    
    #pval.corrected.fdr <- p.adjust(pval, method="fdr")
    
    #res$C_fdr <- pval.corrected.fdr
    
    # fisher for G
    pval <- apply(tbl$G_snps, 2, function(e)
      stats::fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
    
    res$G_fisher <- pval
    
    #pval.corrected.fdr <- p.adjust(pval, method="fdr")
    
    #res$G_fdr <- pval.corrected.fdr
    
    # write result table to disk
    write.table(res$A_fisher, 
                file = paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/results/',
                              segment, '.A_fisher.txt'), col.names = F,
                sep ="\t")
    write.table(res$T_fisher, 
                file = paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/results/',
                              segment, '.T_fisher.txt'), col.names = F,
                sep ="\t")
    write.table(res$C_fisher, 
                file = paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/results/',
                              segment, '.C_fisher.txt'), col.names = F,
                sep ="\t")
    write.table(res$G_fisher, 
                file = paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/results/',
                              segment, '.G_fisher.txt'), col.names = F,
                sep ="\t")
    print(segment)
  }
  print("Completed")
}

fisher_on_genome(phenpath, segments)

# merge fisher results

nucl = 'G'

path <- '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/results'
results <- dir(path, paste0('*',nucl,'_fisher.txt'), full.names = TRUE)

res <- plyr::ldply(results, function(x) read.table(x, sep = '\t'))
colnames(res) = c('Position','Fisher')
res$Position <- rownames(res)

res$FDR <- p.adjust(res$Fisher, method="fdr")

write.table(res, 
            file=paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/',nucl,'_all.txt') ,
            row.names = F, col.names = F, quote = F,
            sep="\t")

###             ###
### CIRCOS plot ###
###             ###

# create CIRCOS output
options(scipen = 999)
options(digits=10)

res$chr1 <- rep('chr1', nrow(res))

res <- res[res$Fisher <= 0.05,]
res$Fisher <- round(res$Fisher, 6)

write.table(cbind(res$chr1, res$Position, res$Position, res$Fisher), 
            file=paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/',nucl,'_fisher.txt') ,
            row.names = F, col.names = F, quote = F,
            sep="\t")

res <- res[res$FDR <= 0.05,]
res$FDR <- round(res$FDR, 6)
write.table(cbind(res$chr1, res$Position, res$Position, res$FDR), 
            file=paste0('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/',nucl,'_FDR.txt') ,
            row.names = F, col.names = F, quote = F,
            sep="\t")

# plot
par(mfrow=c(2,2))

plot(res$Fisher, res$FDR, main="Fisher vs FDR")

heatmap(as.matrix(cbind(res$A_fisher, res$T_fisher, res$C_fisher, res$G_fisher)))

# write result table (needs work)
#write.table(res, 
#            file = '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/results.txt',
#            sep ="\t")

# correlate FDR to CI

CI <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/REF_Consistency.txt',
           header =T)

A_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/A_all.txt')
colnames(A_all) <- c('postion','fisher', 'FDR')
A_all$nucl <- "A"

T_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/T_all.txt')
colnames(T_all) <- c('postion','fisher', 'FDR')
T_all$nucl <- "T"

C_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/C_all.txt')
colnames(C_all) <- c('postion','fisher', 'FDR')
C_all$nucl <- "C"

G_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/G_all.txt')
colnames(G_all) <- c('postion','fisher', 'FDR')
G_all$nucl <- "G"

###             ###
### CIRCOS plot ###
###             ###

#total FDR
all <- rbind(A_all, C_all, G_all, T_all)
all$FDR_ALL <- p.adjust(all$fisher, method="fdr")

A_all <- all[all$nucl == "A",]
T_all <- all[all$nucl == "T",]
C_all <- all[all$nucl == "C",]
G_all <- all[all$nucl == "G",]

# positie  in promotor
#4470140 position 1
#4470150 postion 11

# positie on ampC
# 4470906 position 666 on ampC gene

CI_score <- CI[4470140,][2]
CI[4470150,]
CI[4470906,]

# new table
total <- cbind(CI$Positionnew, CI$ConsistencyIndex, A_all$FDR_ALL, T_all$FDR_ALL, C_all$FDR_ALL, G_all$FDR_ALL)
colnames(total) <- c('position','CI','A_FDR','T_FDR','C_FDR','G_FDR')
total <- as.data.frame(total)

best_positions <- total[total$CI <= as.numeric(CI_score),]

final_positions <- best_positions[best_positions$A_FDR <= 0.05,]
final_positions <- rbind(final_positions , best_positions[best_positions$T_FDR <= 0.05,])
final_positions <- rbind(final_positions , best_positions[best_positions$C_FDR <= 0.05,])
final_positions <- rbind(final_positions , best_positions[best_positions$G_FDR <= 0.05,])

# count other positions with same scores as position 1 in promotor

final_positions <- unique(final_positions)

test <- cbind(final_positions, CI[match(final_positions$position, CI$Positionnew),])

final <- data.frame(test$position, test$CI, test$A_FDR, test$T_FDR, test$C_FDR, test$G_FDR,
           test$CountsACGT, test$MinimumNumberChangesOnTree)

colnames(final) <- c('position', 'CI', 'A_FDR','T_FDR', 'C_FDR', 'G_FDR', 
                     'CountsACGT', 'MinimumNumberChangesonTree')

final <- final[order(final$position),]

write.table(final, file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_final_positions.txt', 
row.names=F, col.names=T, sep='\t')

circos <- data.frame(rep('chr1', nrow(final)), final$position, final$position, final$CI)

write.table(circos, file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_final_positions_circos.txt', 
            row.names=F, col.names=F, sep='\t', quote = F)

### find lowest scoring nucleotide per position ## (Martijn Huijnen 20191109)
total <- data.frame(A_all$postion, A_all$FDR_ALL, C_all$FDR_ALL, G_all$FDR_ALL, T_all$FDR_ALL)

# add minimal value index to column best
total$best <- apply(total, 1, FUN=which.min)

# add minimal value to column min
total$min <- apply(total, 1, FUN=min)

# filter on min FDR <= 0.05
total_filter <- total[total$min <= 0.05,]

total$best[total$best == 2] <- 'A'
total$best[total$best == 3] <- 'C'
total$best[total$best == 4] <- 'G'
total$best[total$best == 5] <- 'T'

circos <- data.frame(rep('chr1', nrow(total)), total$A_all.postion, total$A_all.postion, total$best)

write.table(circos, file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_FDR.txt', 
            row.names=F, col.names=F, sep='\t', quote = F)


circos <- read.table(file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191119_FDR.txt')

###                                   ###
### rebuttal 8 November 2020          ###
### INDEL on promotor and attenuator  ###
### Fisher exact                      ###
###                                   ###

# set working location
#setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/segments")

phenpath <- '/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/20191028_phenotype_CTX_corrected.txt'

# read phenotype data
# read CTX R and S for non-pampC
phenotable <- read.table(phenpath, 
                         sep = "\t", header = F)

# obtain only R and S
phen <- phenotable$V2

segment <- segments[1]

alignment <- ape::read.dna(segment, format = 'fasta', as.character = TRUE)

# load absence/presence INDEL table
INDELS = "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/INDEL_promotor_attenuator/presence-absence_INDELS.csv"
INDELS <- read.table(INDELS, 
                         sep = ",", header = T)
INDELS <- INDELS[3:ncol(INDELS)]
INDELS <- t(INDELS)

# remove samples of which no phenotypic data is provided
INDELS <- INDELS[rownames(INDELS) %in% phenotable$V1,]

# sort data so the alignment matches with the phenotype data
INDELS <- INDELS[match(phenotable$V1, rownames(INDELS)),]

# check if order is correct before proceding
rownames(alignment) == phenotable$V1

tbl <- list()
tbl$INDELS <- INDELS

# fisher for A
pval <- apply(tbl$INDELS, 2, function(e)
  stats::fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)

INDELS_table <- apply(tbl$INDELS, 2, function(e)
  stats::fisher.test(table(factor(e, levels=c(0,1)), phen)))

###                                       ###
### get contingency table of 24 positions ###
###               Fisher exact            ###
###                                       ###

# get alignment
DNAbin <- ape::read.dna('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Alignment/core.full.final.aln',
                        format ='fasta', as.character = TRUE)

# remove samples of which no phenotypic data is provided
DNAbin <- DNAbin[rownames(DNAbin) %in% phenotable$V1,]

# sort data so the alignment matches with the phenotype data
DNAbin <- DNAbin[match(phenotable$V1, rownames(DNAbin)),]

# extract the 24 positions from alignment
twentyfour_positions <- c(810581,810791,815680,824522,828695,830684,830708,830732,
               831564,832152,832287,833451,843887,1863004,1911551,
               1946016,1946067,1946072,1952745,2051214,2051220,
               2057518,2068593,4470140)

# calculate fisher and obtain contingency table
test <- DNAbin[,2057518]

snps <- test
snps[snps == 'g'] <- 1
snps[snps != 1] <- 0

contigency <- table(factor(snps, levels=c(0,1)), phen)
contigency
formatC(stats::fisher.test(contigency)$p.value, format = "e", digits = 3)

###                         ###
### GUBBINS VCF file parse  ###
###                         ###

setwd('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins/')
vcf = read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/gubbins/ampC.summary_of_snp_distribution.vcf',
                 header = F)

# set format for filter
vcf$V6 <- 'PASS'

# set format to GT
vcf$V9 <- 'GT'

###                   ###
### Table 1 piecharts ###
###                   ###

# Apply blank theme
library(scales)
library(ggplot2)
library(RColorBrewer)

# pie plot function for table 1
piefunc <- function(data, legend){
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"),
      text = element_text(size=50),
      legend.position = legend
    )
  
  pie <- ggplot(data, aes(x="", y=percentage, fill=phylogroup)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)
  
  # Use brewer palette
  pie <- pie + scale_fill_brewer(palette = "Accent") + blank_theme +
    theme(axis.text.x=element_blank())
  
  pie
  return(pie)
}

# pampC
pampC <- data.frame(percentage = c(11.1, 13.3, 27.8, 2.2, 31.1, 3.3, 11.1, 0),
                  phylogroup  = c("A", "B1", "B2", "C", "D", "E", "F", "clade IV"))
piefunc(pampC, "none")

# Putative hyperproducers
phyper <- data.frame(percentage = c(1.64, 31.15, 50.82, 8.2, 4.92, 0, 1.64, 1.64),
                     phylogroup =c("A","B1","B2","C","D","E","F","clade IV"))
piefunc(phyper, "none")

# Putative low-level AmpC producers
plow <- data.frame(percentage = c(4.76, 9.52, 38.1, 9.52, 33.3, 4.76, 0, 0),
                   phylogroup =c("A","B1","B2","C","D","E", "F", "clade IV"))
piefunc(plow, "none")

# Total
total <- data.frame(percentage = c(7, 19.2, 37.2, 5.2, 22.1, 2.3, 6.4, 0.6),
                   phylogroup =c("A","B1","B2","C","D","E","F","clade IV"))
piefunc(total, "bottom")


########### make one big file of all positions consistency index + fisher + FDR ##########


CI <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/HomoplasyFinder/REF_Consistency.txt',
                 header =T)

A_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/A_all.txt')
colnames(A_all) <- c('postion','fisher', 'FDR')
A_all$nucl <- "A"

T_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/T_all.txt')
colnames(T_all) <- c('postion','fisher', 'FDR')
T_all$nucl <- "T"

C_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/C_all.txt')
colnames(C_all) <- c('postion','fisher', 'FDR')
C_all$nucl <- "C"

G_all <- read.table('/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/G_all.txt')
colnames(G_all) <- c('postion','fisher', 'FDR')
G_all$nucl <- "G"

#total FDR
all <- rbind(A_all, C_all, G_all, T_all)
all$FDR_ALL <- p.adjust(all$fisher, method="fdr")

A_all <- all[all$nucl == "A",]
T_all <- all[all$nucl == "T",]
C_all <- all[all$nucl == "C",]
G_all <- all[all$nucl == "G",]

# new table
total <- cbind(CI$Positionnew, CI$ConsistencyIndex, CI$CountsACGT, CI$MinimumNumberChangesOnTree, 
               A_all$fisher, A_all$FDR_ALL, 
               C_all$fisher, C_all$FDR_ALL, 
               G_all$fisher, G_all$FDR_ALL,
               T_all$fisher, T_all$FDR_ALL) 

colnames(total) <- c('position','CI', 'CountsACGT', 'MinimumNumberChangesOnTree',
                     'A_pvalue', 'A_FDR',
                     'C_pvalue', 'C_FDR',
                     'G_pvalue', 'G_FDR',
                     'T_pvalue', 'T_FDR')

total <- as.data.frame(total)

write.table(total, file='/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/ReferenceBased/Fisher/total_all.txt', 
            row.names = F, col.names = T, quote = F, sep="\t")

