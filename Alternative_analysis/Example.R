library(paco)
library(ape)

# create two random trees
# perform paco
TreeH = rmtree(1, 172)[[1]]
TreeP = rmtree(2, 172)[[1]]

# perform pcoa with ape

# create association matrix
HP <- matrix( rep( 0, len=29584), nrow = 172,
              dimnames = list(TreeP$tip.label, 
                              TreeP$tip.label))
diag(HP) <- 1

# calculate distance
Hdist <- cophenetic(TreeH)
Pdist <- cophenetic(TreeP)

# prepare paco
D <- prepare_paco_data(Hdist, Pdist, HP)

# calculate using paco
D <- paco::add_pcoord(D, correction ="cailliez")

# calculate using ape::pcoa
test <- ape::pcoa(D$H, "none")

# create the trees based on 
par(mfrow=c(2,3))

# plot MDS plots
plot(D$H_PCo, main="MDS TreeH paco cailliez")
plot(test$vectors, main="MDS TreeP ape::pcoa none")
plot(D$H_PCo[,1], test$vectors[,1], main="comparison of first axis")

# calculate using ape::pcoa now with correction
test <- ape::pcoa(D$H, "cailliez")

# plot MDS plots
plot(D$H_PCo, main="MDS TreeH paco cailliez")
plot(test$vectors.cor, main="MDS TreeP ape::pcoa cailliez")
plot(D$H_PCo[,1], test$vectors.cor[,1], main="comparison of first axis")

