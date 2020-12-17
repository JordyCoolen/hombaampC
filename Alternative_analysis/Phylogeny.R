## 2019 23 April
## Evert and Jordy
## phylogeny ampC genotypes

# install phangorn
#install.packages("phangorn", dependencies=TRUE)
# install via anaconda navigator

library(phangorn)    # load the phangorn library
library(magrittr)
library(phytools)
library(treespace)
library(rlist)
library(phylolm)
#library(adegene)
library(adephylo)
library(adegraphics)
library(rgl)
library(reshape2)

setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/")

promoter_tree_name = "promoter/Model_GTR_bootstraps_1000/20190430_172Strains_Promoter.phy_phyml_tree.txt"
attenuator_tree_name = "attenuater/Model_GTR_bootstraps_1000/20190430_172Strains_Attenuator.phy_phyml_tree.txt"
mulvey_tree_name = "mulvey/Model_GTR_bootstraps_1000/20190328_172Strains_Mulvey.phy_phyml_tree.txt"
chrampc_tree_name = "ampC/Model_GTR_bootstraps_1000/20190226_All_chrampC_172strains.phy_phyml_tree.txt"
coreSNP_tree_name = "WGS_SNP/Model_GTR_bootstraps_5/core_SNPs_matrix.phy_phyml_tree.txt"
coreSNP_fasttree_name = "WGS_SNP/core_SNPs_matrix.tre"
SNPs_fasttree_name = "WGS_SNP/SNPs_all_matrix.tre"

# read newick file
promoter.tree = midpoint.root(read.tree(promoter_tree_name))
attenuator.tree = midpoint.root(read.tree(attenuator_tree_name))
mulvey.tree = midpoint.root(read.tree(mulvey_tree_name))
chrampc.tree = midpoint.root(read.tree(chrampc_tree_name))
chrampc2.tree = midpoint.root(read.tree(chrampc2_tree_name))
coreSNP.tree = midpoint.root(read.tree(coreSNP_tree_name))
coreSNP.fasttree = midpoint.root(read.tree(coreSNP_fasttree_name))
SNPs.fasttree = midpoint.root(read.tree(SNPs_fasttree_name))

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
tiplabels = read.table('WGS_SNP/Tip_labels.txt', sep="\t", header = 1)

# rename tips to ampC ID's
coreSNP.tree = rename.tips(coreSNP.tree, tiplabels)
coreSNP.fasttree = rename.tips(coreSNP.fasttree, tiplabels)
SNPs.fasttree = rename.tips(SNPs.fasttree, tiplabels)

# create function to randomize tiplabels
randomize.tips <- function(n, tree, name){
  trees = list()
  names = list()
  tips <- tree$tip.label
  for (i in 1:n){
    # randomize order of tip labels
    tree$tip.label <- sample(tips)
    names[[i]] <- paste0(name,"_",i)
    trees[[i]] <- tree
  }
  # add names to list
  names(trees) <- names
  return(trees)
}

# generate single tip random tree
n = 10
random.trees <- randomize.tips(n, mulvey.tree, "m")
random.trees2 <- randomize.tips(n, chrampc.tree, "c")
random.trees3 <- rmtree(n, 172, tip.label = chrampc.tree$tip.label)
random.trees3 <- lapply(random.trees3, midpoint.root)
names(random.trees3) <- paste0("r", 1:n, sep = "") 
trees <- c(random.trees, random.trees2, random.trees3)
trees <- list.append(trees,
                     "coreSNP" = coreSNP.tree, 
                     "coreSNP.fasttree" = coreSNP.fasttree,
                     "attenuater" = attenuator.tree,
                     "promoter" = promoter.tree,
                     "mulvey" = mulvey.tree,
                     "wgSNP.fasttree" = SNPs.fasttree,
                     "chrampc" = chrampc.tree,
                     "chrampc2" = chrampc.tree)
                     
# set class
class(trees) <- "multiPhylo"

# plot trees to show randomization
par(mfrow=c(3,2))
plot(trees$mulvey, main="mulvey", cex=0.5)
plot(trees$m_1, main="tip_random_mulvey", cex=0.5)
plot(trees$chrampc, main="chrampC", cex=0.5)
plot(trees$c_1, main="tip_random_chrampC", cex=0.5)
plot(trees$r1, main="random_1", cex=0.5)
plot(trees$r10, main="random_10", cex=0.5)

other_dist_plot <- function(trees, list){
  par(mfrow=c(1,3))
  for (i in list){
    print(i)
    res <- treespace(trees, nf=3, processors=4, method=i)
    res
    plot(hclust(res$D, method = "complete"), main=i)
  }
  par(mfrow=c(1,1))
}

other_dist_plot(trees, c("treeVec", "KF", "RF"))

# calculate treespace
res <- treespace(trees, nf=3, processors=4, lambda=0.5)

# get distance matrix
dist <- res$D
#dist <- cophenetic(trees[[1]])

label <- c()
for (tree in trees){label <- c(label, tree$name)}
par(mfrow=c(1,3))
plot(hclust(dist, method = "average"), labels=label)
plot(hclust(dist, method = "complete"))

# plot 2 dimensions
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5, 
           xax=1, yax=2, scree.size=0.2)
screeplot(res$pco)
round(res$pco$eig/sum(res$pco$eig)*100, 1)

# extract Distancematrix
test <- as.matrix(res$D)

# chrampc distance to random mulvey
one = test[rownames(test) == "chrampc", 1:n]
# mulvey distance to random chrampc
two = test[rownames(test) == "mulvey", (n+1):(n*2)]
# chrampc distance to random trees
three = test[rownames(test) == "chrampc", (n*2+1):(n*3)]
# mulvey distance to random trees
four = test[rownames(test) == "mulvey", (n*2+1):(n*3)]
# distance between chrampc and mulvey
five = test[rownames(test) == "chrampc", (n*3)+1]
# distance between mulvey and chrampc (same as previous)
six = test[rownames(test) == "mulvey", (n*3)+2]

# boxplot of differences between random and other distance
boxplot(one, two, three, four, five, six,
        names = c("chrampC_vs_rmulvey","mulvey_vs_rchrampC",
                  "chrampC_vs_random","mulvey_vs_random",
                  "chrampC_vs_mulvey","mulvey_vs_chrampC"),
        las = 2)

boxdata <- melt(cbind(one, two, three, four, five, six))

library(ggplot2)
# Basic violin plot
p <- ggplot(data, aes(x=Var2, y=value, fill=Var2)) + geom_violin()
p + geom_boxplot(width=0.1) + theme_minimal()

########## PAco ############
library(paco)
library(ape)
library(vegan)
library(gtools)

setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/Core_genes/")
setwd("/Volumes/Macintosh HD/Users/jordycoolen/surfdrive/PhD_RUMC/Ecoli_pampC_genotypes/Data/Core_genes/test")


# PACo analysis
# read trees

set.seed(5)

# edit 2019 22 Mei (Tom Ederveen) core vs core fasttree
FastTreeCore <- read.tree(coreSNP_fasttree_name)
FastTreeCore = rename.tips(FastTreeCore, tiplabels)
p = Pipeline_PACo(TreeCore, FastTreeCore, 172, "Core", "Core_FastTree", 100)

# get all tre files and combinations to compare
file.names <- dir(getwd(), pattern =".tre")
#file.names <- file.names[1:4]
# check if tre files contain all samples
core.names <- c()
for (label in file.names){
  t <- read.tree(label)
  if (length(t$tip.label) != 172){
    file.remove(label)
  }
  else{
    core.names <- c(core.names, label)
  }
}

core.names <- file.names
core.names <- core.names[200:277]
# file.names <- file.names[1:4] # take subset
combis <- combinations(length(core.names), 2, core.names)
dim(combis)

#  correct tiplabels
options(stringsAsFactors = FALSE)
tiplabels = read.table('Tip_labels.txt', sep="\t", header = 1)

# get all m2 results
result <- data.frame()
for (i in 1:nrow(combis)){
  row <- combis[i,]
  TreeH <- read.tree(row[1])
  TreeH <- roary_to_R(TreeH)
  TreeH <- rename.tips(TreeH, tiplabels)
  TreeH$name <- row[1]
  TreeP <- read.tree(row[2])
  TreeP <- roary_to_R(TreeP)
  TreeP <- rename.tips(TreeP, tiplabels)
  TreeP$name <- row[2]
  p = Pipeline_PACo(TreeH, TreeP, 172, TreeH$name, TreeP$name, 1, "pdf", F)
  result <- rbind(result, p$pwc)
}

# change headers
colnames(result) <- c("X1","X2","m2")

# get node names
nodes_names <- unique(c(result$X1,result$X2))

# set up weight matrix
m <- matrix(0, nrow=length(nodes_names), ncol=length(nodes_names))
rownames(m) <- colnames(m) <- nodes_names

# populate edge weights in appropriate matrix elements
for(dd in 1:nrow(result)) {
  row_id <- result[dd, "X1"]
  col_id <- result[dd, "X2"]
  m[row_id, col_id] <- as.numeric(result[dd, "m2"])
  m[col_id, row_id] <- as.numeric(result[dd, "m2"])
}

# creating aggregation results
plot(hclust(as.dist(m)))
plot(hclust(as.dist(m), method="average"))
biplot(pcoa(m))
heatmap(m, cexRow = 0.8, cexCol = 0.8)

# alternative
# for combinations
# https://stat.ethz.ch/R-manual/R-devel/library/utils/html/combn.html

#https://stackoverflow.com/questions/4493287/generating-a-very-large-matrix-of-string-combinations-using-combn-and-bigmemor
".combinadic" <- function(n, r, i) {
  
  # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
  # http://en.wikipedia.org/wiki/Combinadic
  
  if(i < 1 | i > choose(n,r)) stop("'i' must be 0 < i <= n!/(n-r)!")
  
  largestV <- function(n, r, i) {
    #v <- n-1
    v <- n                                  # Adjusted for one-based indexing
    #while(choose(v,r) > i) v <- v-1
    while(choose(v,r) >= i) v <- v-1        # Adjusted for one-based indexing
    return(v)
  }
  
  res <- rep(NA,r)
  for(j in 1:r) {
    res[j] <- largestV(n,r,i)
    i <- i-choose(res[j],r)
    n <- res[j]
    r <- r-1
  }
  res <- res + 1
  return(res)
}

# treespace alternative

library(rlist)
library(phytools)
trees = list()
for (i in file.names){
  tree <- rename.tips(roary_to_R(midpoint.root(read.tree(i))), tiplabels)
  tree$name <- i
  trees <- list.append(trees, i = tree)
}

class(trees) <- "multiPhylo"

# calculate treespace
res <- treespace(trees, nf=3, processors=4, lambda=0.5, method = "KF")

# get distance matrix
dist <- res$D
#dist <- cophenetic(trees[[1]])

label <- c()
for (tree in trees){label <- c(label, tree$name)}
par(mfrow=c(1,3))
plot(hclust(dist, method = "average"), labels=label)
plot(hclust(dist, method = "complete"), labels=label)

# test core genes 2019 23 Mei
adkname = "adk.tre"
Treeadk = read.tree(adkname)
Treeadk$name <- adkname

ubiEname = "ubiE.tre"
TreeubiE = read.tree(ubiEname)

thiLname = "thiL.tre"
TreethiL = read.tree(thiLname)

soxRname = "soxR.tre"
TreesoxR = read.tree(soxRname)

mutLname = "mutL.tre"
TreemutL = read.tree(mutLname)

mutSname = "mutS.tre"
TreemutS = read.tree(mutSname)

# replace "_genenumber" to "" inorder to match genename
roary_to_R <- function(Tree){
  l = c()
  for (n in Tree$tip.label){
    n <- gsub("_(.+)","",n, perl=TRUE)
    #print(n)
    l = c(l, n)
  }
  Tree$tip.label <- l
  return(Tree)
}

Treeadk <- roary_to_R(Treeadk)
Treeadk <- rename.tips(Treeadk, tiplabels)
p = Pipeline_PACo(TreeCore, Treeadk, 172, "Core", "adk", 100)

TreeubiE <- roary_to_R(TreeubiE)
TreeubiE <- rename.tips(TreeubiE, tiplabels)
p = Pipeline_PACo(TreeCore, TreeubiE, 172, "Core", "ubiE", 100)

TreethiL <- roary_to_R(TreethiL)
TreethiL <- rename.tips(TreethiL, tiplabels)
p = Pipeline_PACo(TreeCore, TreethiL, 172, "Core", "thiL", 100)

TreesoxR <- roary_to_R(TreesoxR)
TreesoxR <- rename.tips(TreesoxR, tiplabels)
p = Pipeline_PACo(TreeCore, TreesoxR, 172, "Core", "soxR", 100)

TreemutL <- roary_to_R(TreemutL)
TreemutL <- rename.tips(TreemutL, tiplabels)
p = Pipeline_PACo(TreeCore, TreemutL, 172, "Core", "mutL", 100)

TreemutS <- roary_to_R(TreemutS)
TreemutS <- rename.tips(TreemutS, tiplabels)
p = Pipeline_PACo(TreeCore, TreemutS, 172, "Core", "mutS", 100)

Trees <- c(Treeadk, TreeubiE, TreethiL, TreesoxR, TreemutL, TreemutS)
nTrees <- seq(1, length(Trees)) 
require(gtools)
combis <- combinations(n=length(nTrees),r=2,v=nTrees)




test <- data.frame(rbind(p$pwc,c("Core","mutS", 0), c("mutS","mutL",5)))
test <- reshape(test, direction="wide", idvar="X2", timevar="X1")



df <- read.table(textConnection("
A1  A1  0.90
                                A1  B1  0.85
                                A1  C1  0.45
                                A1  D1  0.96
                                B1  B1  0.90
                                B1  C1  0.85
                                B1  D1  0.56
                                C1  C1  0.55
                                C1  D1  0.45
                                D1  D1  0.90"))


require(adephylo)
par(mfrow=c(1,2))
Hdist <- distTips(TreeCore, tips = "all", method = c("patristic"), useC = TRUE)
Hdist <- sqrt(distH)
hc_H <- hclust(as.dist(distH))
plot(hc_H, hang = -1, cex = 0.2, labels= FALSE, main = paste0("Distance_Tree ", 'nameH'))
Pdist <- distTips(TreesoxR, tips = "all", method = c("patristic"), useC = TRUE)
Pdist <- sqrt(distP)
hc_P <- hclust(as.dist(distP))
plot(hc_P, hang = -1, cex = 0.2, labels= FALSE, main = paste0("Distance_Tree ", 'nameP'))
D <- prepare_paco_data(Hdist, Pdist, HP)
D <- add_pcoord(D, correction ="cailliez")


TreechrampC <- read.tree(chrampc_tree_name)
TreeRandom <- rmtree(1, 172, tip.label = TreechrampC$tip.label)[[1]]
Treemulvey <- read.tree(mulvey_tree_name)
Treeattenuator <- read.tree(attenuator_tree_name)
Treepromoter <- read.tree(promoter_tree_name)

# coreSNP tree
TreeCore <- read.tree(coreSNP_tree_name)
TreeCore = rename.tips(TreeCore, tiplabels)

Trees = c(TreechrampC, Treepromoter, Treeattenuator, Treemulvey, TreeCore, TreeRandom)

results = data.frame()

# chrampC as TreeH
p = Pipeline_PACo(TreechrampC, Treepromoter, 172, "chrampC", "promoter", 10)
results <- cbind(p)
p = Pipeline_PACo(TreechrampC, Treeattenuator, 172, "chrampC", "attenuator", 100)
results <- cbind(results, p)
p = Pipeline_PACo(TreechrampC, Treemulvey, 172, "chrampC", "mulvey", 100)
results <- cbind(results, p)
p = Pipeline_PACo(TreechrampC, TreeCore, 172, "chrampC", "Core", 100)
results <- cbind(results, p)
p = Pipeline_PACo(TreechrampC, TreeRandom, 172, "chrampC", "random", 100)
results <- cbind(results, p)

# attenuator as TreeH
p = Pipeline_PACo(Treeattenuator, Treepromoter, 172, "attenuator", "promoter", 100)
results <- cbind(results, p)
p = Pipeline_PACo(Treeattenuator, Treemulvey, 172, "attenuator", "mulvey", 100)
results <- cbind(results, p)
p = Pipeline_PACo(Treeattenuator, TreeCore, 172, "attenuator", "Core", 100)
results <- cbind(results, p)
p = Pipeline_PACo(Treeattenuator, TreeRandom, 172, "attenuator", "random", 100)
results <- cbind(results, p)

# promoter as TreeH
p = Pipeline_PACo(Treepromoter, Treemulvey, 172, "promoter", "mulvey", 100)
results <- cbind(results, p)
p = Pipeline_PACo(Treepromoter, TreeCore, 172, "promoter", "coreSNP", 100)
results <- cbind(results, p)
p = Pipeline_PACo(Treepromoter, TreeRandom, 172, "promoter", "random", 100)
results <- cbind(results, p)

# mulvey as TreeH
p = Pipeline_PACo(Treemulvey, TreeRandom, 172, "mulvey", "random", 100)
results <- cbind(results, p)

# add padjust values to table
results <- rbind(results, padjust)

####### multiple testing correction #######
# FDR correction on results
padjust <- p.adjust(results[4,], method ="BH")

# write results
require(openxlsx)
write.xlsx(results, file= 'PACo_results.xlsx', row.names=TRUE)

################################################################################

# method in development to see residuals and annotation
rs <- data.frame(residuals(p$proc))

# subset data to see which samples have large procrustes error
above_mean <- subset(rs, residuals.p.proc. > mean(rs$residuals.p.proc.))
above_mean <- cbind(above_mean, rownames(above_mean))
colnames(above_mean) <- c('residuals', 'names')

require(openxlsx)
metadata <- "20190409_Lijst_AmpC_project_data_promoter_attenuator_correctie_mulvey_typen.xlsx"
metadata <- read.xlsx(metadata)

subset(metadata, unique_ID == rownames(above_mean))

merged <- merge(above_mean, metadata, by.x = "names", by.y = "unique_ID")

rownames(above_mean)==metadata$unique_ID

# protest analysis as alternative for the PACo permutation
ptest <- protest(p$H_PCo, p$P_PCo, permutations = 10000, symmetric = T)






############# ############# ############# ############# 

# TreeDiff plot
plotTreeDiff(midpoint.root(TreeH), midpoint.root(TreeP), use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", cex=0.5)

HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)
colnames (SQres.jackn) <- paste(rownames(HostX),rownames(ParY), sep="-")
t.critical = qt(0.975,NLinks-1)
for(i in c(1:NLinks))
{HP.ind <- HP
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(host.D, para.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind}
SQres.jackn <- SQres.jackn**2
SQres <- (residuals (HP.proc)**2)
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres))
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks)
pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.25, col="white", xlab="Host-parasite link", 
                   ylab= "Squared residuals", ylim=c(0, max(phi.UCI)), cex.lab=1.2)
text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels = colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.6)
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2)

# test on random tree to show relevance of the study
# in a more understandable way


########## TEST ############

bla <- rmtree(3,172)
bla <- lapply(bla, midpoint.root)
bla.randomtrees <- randomize.tips(100, bla[[1]], "x")
bla <- c(bla, bla.randomtrees)

class(bla) <- "multiPhylo"
bla.res <- treespace(bla, nf=5, processors=4)

# plot 2 dimensions
plotGroves(bla.res$pco, lab.show=TRUE, lab.cex=1.5)

########## END ############

# create random trees
random.trees = rmtree(39, 172, tip.label=tiplabels$ampC_ID)

# todo random label assignment

# plot tree
plot(mulvey.tree, no.margin=TRUE, edge.width=2)
plot(SNPs.tree, no.margin=TRUE, edge.width=2)
plot(coreSNP.tree, no.margin=TRUE, edge.width=2)

# plot random tree
plot(random.trees[[1]], no.margin=TRUE, edge.width=2)

# create supertree
rf_st <- superTree(random.trees, method = "RF")

# add tree to random
random.trees[[30]] <- midpoint.root(random.trees[[30]])
random.trees[[40]] <- midpoint.root(coreSNP.tree)
random.trees[[41]] <- midpoint.root(SNPs.tree)
random.trees[[42]] <- midpoint.root(mulvey.tree)
random.trees[[43]] <- midpoint.root(rf_st)

# calculate distance
dist = RF.dist(random.trees) # pairwise distance

# heatmap test
heatmap(as.matrix(dist))

# minimal spanning tree
MST <- mst(dist)
plot(MST)

# use of treespace
#install.packages("treespace")
#install.packages("phylolm")


# test trees
test.trees <- rmtree(10, 20)

# use treespace
res <- treespace(random.trees, nf=3)
names(res)

# table.value with some customization
table.value(res$D, nclass=5, method="color", symbol="circle", col=redpal(5))
table.image(res$D)

# heatmap on treespace matrix
heatmap(as.matrix(res$D))

# plot 2 dimensions
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)

# Compare median trees from clusters 1 and 2:
plotTreeDiff(random.trees[[40]],random.trees[[41]], use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)

# phylolm

# read data
require(openxlsx)
traits = read.xlsx('20190409_Lijst_AmpC_project_data_promoter_attenuator_correctie_mulvey_typen.xlsx')
rownames(traits) <- traits$unique_ID

# rename function used to rename all pampC variant types to pampC
# needed to run the models
rename_pampC <- function(data) {
  data$ampC.type[data$ampC.type != "negative"] <- "ampC"
  data$ampC.type <- factor(data$ampC.type)
  return(data)
}

traits <- rename_pampC(traits)

# model fit
OUshifts(traits$CAZ, SNPs.tree, method = c("mbic"), 10, check.pruningwise = TRUE)


phyloglm(promoter.mutation~CAZ, traits, SNPs.tree, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot = 0, full.matrix = TRUE)

rTrait(n=1, SNPs.tree, model=c("BM"),
       parameters = NULL, plot.tree=TRUE)

# other package
library(phylopath)
models <- define_model_set(
  one   = c(ampC.type ~ CAZ),
  two   = c(ampC.type ~ CTX),
  three = c(ampC.type ~ FOX)
)

plot_model_set(models)

result <- phylo_path(models, data = traits, tree = SNPs.tree, model = 'lambda')
(s <- summary(result))
plot(s)

# phylosignal
install.packages("phylosignal")
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)

trts <- cbind(traits$CAZ, traits$CTX, traits$FOX)
rownames(trts) <- rownames(traits)

p4d <- phylo4d(mulvey.tree, trts)

barplot(p4d)
dotplot(p4d)

