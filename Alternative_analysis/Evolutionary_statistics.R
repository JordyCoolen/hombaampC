#############################
####                     ####
####   JPM Coolen        ####
####   02 Aug 2019       ####
####                     ####
####  (p)AmpC genotype   ####
####  Evolutionary stats ####
#############################


# installing packages
#install.packages("pegas")
#install.packages("spider")
#install.packages("ggseqlogo")
require("pegas")
require("spider")
require("ggplot2")
setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/")

# core SNP data
core_SNP = "WGS_SNP/core_SNPs_matrix_renamed.fa"
core_Roary = "Roary_output/core_gene_alignment_renamed.aln"
promotor = "promoter/20190430_172Strains_Promoter.aln"
ampC = "ampC/20190226_All_chrampC_172strains.aln"
attenuator = "attenuater/20190430_172Strains_Attenuator.fasta"

# calculate tajimaD on fragment
TD <- function(filepath, name){
  # read alignment files core
  dna <- ape::read.dna(filepath, format="fasta")
  
  # delete gaps in alignment (experimental)
  dna <- del.colgapsonly(dna, threshold = 0.00001, freq.only = FALSE)
  
  res <- tajima.test(dna)
  
  res <- c(res, name)
  
  return(res)
}

# function to test Tajima's D on window
window <- function(filepath, window, name){
  
  # read alignment files core
  dna <- ape::read.dna(filepath, format="fasta")
  
  # delete gaps in alignment (experimental)
  dna <- del.colgapsonly(dna, threshold = 0.00001, freq.only = FALSE)
  
  # calculate tajima's D statistic in window
  tj_test_window = lapply(slidingWindow(dna, window, interval = window),
                          tajima.test)
  
  # create loop to obtain all window tajima D test statistics
  tajima = vector()
  position = vector()
  pos = 0
  
  # loop to extract stat vs windowsize
  for (i in tj_test_window){
    tajima <- c(tajima, i$D)
    position <- c(position, pos)
    df <- data.frame(tajima, position, name)
    pos <- pos + window
  }
  
  return(df)
}

# create table of Tajima's D values
if(!file.exists('Routput/TajimasD/Tajima_final_table.txt')){
  # calculate Tajima'D on fragment
  t <- TD(filepath=core_SNP, 'core_SNP')
  t2 <- TD(filepath=core_Roary, 'core_Roary')
  t3 <- TD(filepath=promotor, 'promotor')
  t4 <- TD(filepath=ampC, 'ampC')
  t5 <- TD(filepath=attenuator, 'attenuator')
  
  t_final <- rbind(t, t2)
  t_final <- rbind(t_final, t3)
  t_final <- rbind(t_final, t4)
  t_final <- rbind(t_final, t5)
  
  colnames(t_final) <- c("D","Pval.normal","Pval.beta","name")
  t_final <- as.matrix(t_final)
  
  write.table(t_final, file = "Routput/TajimasD/Tajima_final_table.txt", sep="\t")
} else {
  t_final <- read.table(file= "Routput/TajimasD/Tajima_final_table.txt")
}

# plot the table
ggplot(t_final, aes(x=name, y=D,
                     mainTitle="Tajima's D score over window")) +
  geom_point() +
  theme_classic()

# create table of Tajima's D using sliding window
if(!file.exists('Routput/TajimasD/Tajima_final_table_slidingwindow.txt')){
  # calculate window
  df <- window(filepath=core_SNP, 1000, "core_SNP_1000")
  df2 <- window(filepath=core_Roary, 5000, "core_Roary_5000")
  df3 <- window(filepath=promotor, 10, "promotor_10")
  df4 <- window(filepath=ampC, 50, "ampC_50")
  df5 <- window(filepath=attenuator, 5, "attenuater_5")
  
  df_final <- rbind(df, df2)
  df_final <- rbind(df_final, df3)
  df_final <- rbind(df_final, df4)
  df_final <- rbind(df_final, df5)
  
  df_final <- as.matrix(df_final)
  write.table(df_final, file = "Routput/TajimasD/Tajima_final_table_slidingwindow.txt", sep="\t")
} else {
  df_final <- read.table(file= "Routput/TajimasD/Tajima_final_table_slidingwindow.txt")
}

# TODO
ggplot(data=df_final, aes(x=position, y=tajima)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_minimal()

# Tajima's D boxplot plot sliding window
ggplot(df_final, aes(x=name, y=tajima, addDot=TRUE , dotSize=1 ,
                     mainTitle="Tajima's D score over window")) +
  geom_boxplot(width=0.5) + 
  geom_point(alpha=0.3, color="gray", position = "jitter") +
  theme_classic()

# Tajima's D violin plot sliding window
ggplot(df_final, aes(x=name, y=tajima, addDot=TRUE , dotSize=1 ,
                     mainTitle="Tajima's D score over window")) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + 
  geom_point(alpha=0.3, color="gray", position = "jitter") +
  theme_classic()

###### HomoplasyFinder #######

# Install devtools package - necessary for installing homoplasyFinder
#install.packages("devtools")

# Load devtools
#require("devtools")

# Install homoplasyFinder
#install_github("JosephCrispell/homoplasyFinder")

require(homoplasyFinder)
setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/")
# load original tree
#core_SNP.name <- "WGS_SNP/core_SNPs_matrix.tre"
#core_SNP.tree <- ape::read.tree(core_SNP.name)

# rename tips to make compatible with alignments
#core_SNP.tree = rename.tips(core_SNP.tree, tiplabels)

# store new newick tree with renamed tips to disk
#write.tree(core_SNP.tree, file = 'WGS_SNP/core_SNPs_matrix_renamed.tre')

# rename tips
#core_Roary.name = "Roary_output/core_gene_alignment.tre"
#core_Roary.tree <- ape::read.tree(core_Roary.name)
#core_Roary.tree = rename.tips(core_Roary.tree, tiplabels)
#write.tree(core_Roary.tree, file = 'Roary_output/core_gene_alignment_renamed.tre')

# reroot tree (test to see if results change)
#core_Roary.tree <- ape::read.tree(core_Roary.name)
#core_Roary.tree <- root(core_Roary.tree, "ampC_0128")
#write.tree(core_Roary.tree, file = 'Roary_output/core_gene_alignment_renamed_root_ampC_0128.tre')

#phylotools::rename.fasta(core_SNP, tiplabels)

# location to renamed newick tree
core_SNP.name <- "WGS_SNP/core_SNPs_matrix_renamed.tre"
core_Roary.name = "Roary_output/core_gene_alignment_renamed.tre"
#core_Roary.name = "Roary_output/core_gene_alignment_renamed_root_ampC_0128.tre"

setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/")

# create list of alignments
list <- c(core_SNP, ampC, promotor, attenuator)
names(list) <- c("core_SNP", "ampC", "promotor", "attenuator")

# homoplasy test function on multiple alignment
# aggregates the results in a data.frame
homoplasy_test <- function(list, tree, all_positions=FALSE){
  
  final_table = c()
  
  for (name in names(list)){
    
    print(name)
    print(list[name])
    
    dir.create(paste0(getwd(),"/HomoplasyFinder/",name))
    
    inconsistentPositions <- runHomoplasyFinderInJava(tree, list[name],  
                                                      paste0(getwd(),"/HomoplasyFinder/",name,"/"),
                                                      includeConsistentSitesInReport=all_positions)
    # read consistentPosition Table
    consis.table <- read.csv2(paste0("HomoplasyFinder/", name, 
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

if(!file.exists("HomoplasyFinder/markers_HomoplasyFinder_output.txt")){
  markers <- homoplasy_test(list, core_SNP.name, all_positions = TRUE)
  write.table(markers, file="HomoplasyFinder/markers_HomoplasyFinder_output.txt", sep="\t")
} else {
  markers <- read.table(file="HomoplasyFinder/markers_HomoplasyFinder_output.txt")
}

# homoplasy index 1- consistency index

markers_filtered <- markers[markers$ConsistencyIndex != 1,]

# Violin plot of ConsistencyIndex per sample
ggplot(markers_filtered, aes(x=name, y=ConsistencyIndex, addDot=TRUE , dotSize=1 ,
                            mainTitle="ConsistencyIndex")) +
  geom_violin(scale="width", trim=FALSE, fill="gray", color="darkred")+
  #coord_flip() +
  geom_boxplot(width=0.1) +
  geom_point(alpha=0.3, color="tomato", position = "jitter") +
  theme_classic()

# get 1% quantile value
q <- as.data.frame(quantile(markers[markers$name == "core_SNP",]$ConsistencyIndex, probs= c(0.01)))
colnames(q) <- c('q')

# plot density plot with values in regions and quantile values
p <- ggplot(markers[markers$name == "core_SNP",], 
            aes(x=ConsistencyIndex,
            mainTitle="ConsistencyIndex")) +
  geom_density() +
  geom_point(data=markers, aes(x=ConsistencyIndex, y=name)) +
  geom_vline(data = q, aes(xintercept=q),
             color="blue", linetype="dashed", size=1) +
  scale_x_reverse() +
  theme_minimal()
p

# select only 0.50 - 0.00 (to make visual better)
markers_0.50_0.00 <- markers[markers$ConsistencyIndex<=0.5,]

library(plotly)
# plot density plot with values in regions and quantile values
p <- ggplot(markers_0.50_0.00[markers_0.50_0.00$name == "core_SNP",], 
            aes(x=ConsistencyIndex,
            mainTitle="ConsistencyIndex")) +
  geom_density() +
  geom_point(data=markers_0.50_0.00, aes(x=ConsistencyIndex, y=name)) +
  geom_vline(data = q, aes(xintercept=q),
             color="blue", linetype="dashed", size=1) +
  scale_x_reverse() +
  theme_minimal()
p
ggplotly(p)


# count number of < 1% quantile

count_q

count_q <- cbind('ampC', 
                 summary(markers[markers$name == "ampC",]$ConsistencyIndex < as.numeric(q))[3])
count_q <- rbind(count_q, cbind('core_SNP' , 
                                summary(markers[markers$name == "core_SNP",]$ConsistencyIndex < as.numeric(q))[3]))
count_q <- rbind(count_q, cbind('attenuator', 
                                summary(markers[markers$name == "attenuator",]$ConsistencyIndex < as.numeric(q))[3]))
count_q <- rbind(count_q, cbind('promotor', 
                                summary(markers[markers$name == "promotor",]$ConsistencyIndex < as.numeric(q))[3]))
# counts lower then quantile
count_q


########                      ########
########       Core Roary     ########
########                      ########


# if file exists for dividing core Roary results
setwd("/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/TEST")
if(!file.exists("core_Roary_HomoplasyFinder_output.txt")){
  # dividing the core_Roary sequences to run HomoplasyFinder
  core_Roary_DNAbin <- ape::read.dna(core_Roary, format ='fasta')
  core_Roary.length <- length(core_Roary_DNAbin[1,])
  
  # divide core Roary alignment in parts
  s <- seq(1, core_Roary.length, by = 50000)
  parts = c()
  
  # create the core Roary alignment parts and save to disk
  for (i in seq(1, length(s))){
    start = s[i]
    end = s[i+1] - 1
    if (!is.na(end)){
      print(start)
      print(end)
      dna <- core_Roary_DNAbin[,start:end]
      ape::write.FASTA(dna, 
                       file = paste0(start,"_", end, "_core_Roary.aln"))
    } else {
      print(start)
      print(core_Roary.length)
      dna <- core_Roary_DNAbin[,start:core_Roary.length]
      ape::write.FASTA(dna, 
                       file = paste0(start,"_", core_Roary.length, "_core_Roary.aln"))
    }
  }
  
  # create list of created part of Core Roary
  file.names <- dir(getwd(), pattern =".aln")
  for(i in 1:length(file.names)){
    print(file.names[i])
  }
  
  # add names to the values (same as the value)
  names(file.names) <- file.names
  core_SNP.full <- "/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/WGS_SNP/core_SNPs_matrix_renamed.tre"
  roary_table <- homoplasy_test(file.names, core_SNP.full, all_positions = TRUE)
  write.table(roary_table, file='core_Roary_HomoplasyFinder_output.txt', sep="\t")
} else {
  roary_table <- read.table(file="core_Roary_HomoplasyFinder_output.txt")
}

# manually rename the name files to only numbers
# regex _(.*)\w+.aln
# reorder rows
roary_table <- roary_table[order(roary_table$name),]
# renumber rownames to correct positions
rownames(roary_table) <- seq(length=nrow(roary_table))

# convert rownames to Position
roary_table$Position <- rownames(roary_table)
roary_table$name <- as.character(roary_table$name)

# get rid of all consistency values of 1 positions
roary_table_filtered <- roary_table[roary_table$ConsistencyIndex != 1,]

# get 1% quantile value (based on Core Roary)
#q <- as.data.frame(quantile(roary_table$ConsistencyIndex, probs= c(0.01)))
#colnames(q) <- c('q')

# plot density plot with values in regions and quantile values
p <- ggplot(roary_table, aes(x=ConsistencyIndex, mainTitle="ConsistencyIndex")) +
  geom_density() +
  geom_point(data=markers, aes(x=ConsistencyIndex, y=name)) +
  geom_vline(data = q, aes(xintercept=q),
             color="blue", linetype="dashed", size=1) +
  scale_x_reverse() +
  theme_minimal()
p

########                      ########
########  Combine everything  ########
########                      ########

roary_table$name <- 'core_Roary'

final_table <- rbind(markers, roary_table)

final_table_unique <- unique(final_table[c("ConsistencyIndex", "name")])

final_table_filtered <- final_table[final_table$ConsistencyIndex != 1,]

# plot density plot with values in regions and quantile values
p <- ggplot(final_table[final_table$name == "core_Roary",], 
            aes(x=ConsistencyIndex, mainTitle="ConsistencyIndex")) +
  geom_density() +
  geom_point(data=final_table_unique, 
             aes(x=ConsistencyIndex, y=name)) +
  geom_vline(data = q, aes(xintercept=q),
             color="blue", linetype="dashed", size=1) +
  scale_x_reverse() +
  theme_minimal()
p


# select only 0.25 - 0.00 (to make visual better)
final_table_0.25_0.00 <- final_table[final_table$ConsistencyIndex<=0.25,]

p <- ggplot(final_table_0.25_0.00[final_table_0.25_0.00$name == "core_Roary",], 
            aes(x=ConsistencyIndex, mainTitle="ConsistencyIndex")) +
  geom_density() +
  geom_point(data=final_table_unique, 
             aes(x=ConsistencyIndex, y=name)) +
  geom_vline(data = q, aes(xintercept=q),
             color="blue", linetype="dashed", size=1) +
  scale_x_reverse() +
  theme_minimal()
p

# count < 1% quantile in core Roary
cbind('core_Roary', 
      summary(final_table[final_table$name == "core_Roary",]$ConsistencyIndex < as.numeric(q))[3])

###########                                    ###########
########### sliding window of consistencyindex ###########
###########                                    ###########

# function to perform sliding window
# https://www.rdocumentation.org/packages/evobiR/versions/1.1/topics/SlidingWindow
require(evobiR)

windows = c(2,3,5,11,22,39,50,100,250,500,1000,2000,5000)

windows = c(11)

count_q_in_window <- function(list, q, table){
  result <- c('windowsize', 'total','count')
  for (i in list){
    print(i)
    # count consistency index lower or equal to 1% quartile
    quant_1 <- function(x) {
      return(length(x[x < as.numeric(q)])) # <=0.06250000 min consistency index of promotor
      # < 0.08333333 1% quantile core SNP
    }

    # perform sliding window on counts
    res <- evobiR::SlidingWindow(FUN=quant_1,
                                 table$ConsistencyIndex,
                                 window = i, step = 1)

    # results to dataframe
    res <- as.data.frame(res)
    res$window <- rownames(res)

    l <- cbind(i, dim(res)[1], summary(res$res >= 2)[3])
    print(l)
    result <- rbind(result, l)
  }
  return(res)
}

out <- count_q_in_window(windows, q, roary_table)


ggplot(data=res, aes(x=res)) +
  geom_violin()

ggplot(data=out, aes(x=window, y=res)) +
  geom_point()

# count < 1% quartile scores

# plot counts

# problem is because of arbritary order of CDS sequences in Roary
# some of the junctions between genes could give a unreal consistency value

# renaming of the trees will not work because
tree <- readAnnotatedTree(paste0(getwd(),"/HomoplasyFinder/promotor/"))

# get only results of one alignment
consis.promotor <- final_table_all[final_table_all$name=="promotor",]

# plot Annotatedtree
plotAnnotatedTree(tree, consis.promotor$Position, promotor)

# barplot of ConsistencyIndex for the promotor
ggplot(data=consis.promotor, aes(x=Position, y=ConsistencyIndex)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  #geom_text(aes(label=round(ConsistencyIndex,2)), vjust=1.6, color="black", size=3) +
  theme_minimal()


library(ggthemes)
scale_colour_Publication()
theme_Publication()


########### TREESPACE ##############

# Treespace to compare trees
require(treespace)
require(phytools)
require(adegraphics)
core_SNP.tree <- ape::read.tree("WGS_SNP/core_SNPs_matrix_renamed.tre")
core_Roary.tree <- ape::read.tree("Roary_output/core_gene_alignment_renamed.tre")
ampC.tree <- ape::read.tree("ampC/Model_GTR_bootstraps_1000/20190226_All_chrampC_172strains.phy_phyml_tree.txt")
promotor.tree <- ape::read.tree("promoter/Model_GTR_bootstraps_1000/20190430_172Strains_Promoter.phy_phyml_tree.txt")
attenuator.tree <- ape::read.tree("attenuater/Model_GTR_bootstraps_1000/20190430_172Strains_Attenuator.phy_phyml_tree.txt")

Trees <- c(midpoint.root(core_SNP.tree), midpoint.root(core_Roary.tree), 
           midpoint.root(ampC.tree), midpoint.root(promotor.tree), midpoint.root(attenuator.tree))
names(Trees) <- c("core_SNP","core_Roary","ampC","promotor","attenuator")
res <- treespace(Trees, nf = 3)
# table.image
table.image(res$D, nclass=30)

# table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(5))

#PCo plot
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)
