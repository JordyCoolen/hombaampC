Pipeline_PACo <- function(TreeH, TreeP, n, nameH, nameP, perm, out ,plot = T){
  ############# alternative method PACo #############
  # example of running
  # Pipeline_PACo(TreeH, TreeP, 172, "chrampC", "promoter", 2)
  
  # check if labels are equal
  check <- setequal(TreeH$tip.label, TreeP$tip.label)
  if(check != T){
    stop(paste0('not the same tiplabels ',nameH, ' ',nameP))
  }
  
  # create association matrix
  HP <- matrix( rep( 0, len=n*n), nrow = n,
                dimnames = list(TreeP$tip.label, 
                                TreeP$tip.label))
  diag(HP) <- 1
  dim(HP)
  
  # check HP
  #heatmap(HP)
  
  # host and parasite trees are then transformed 
  # for further analysis into respective matrices of
  # patristic distances (host.D and para.D):
  Hdist <- cophenetic(TreeH)
  Pdist <- cophenetic(TreeP)
  
  hc_H <- hclust(as.dist(Hdist))
  hc_P <- hclust(as.dist(Pdist))
  
  # create class with distance
  D <- prepare_paco_data(Hdist, Pdist, HP)
  
  # # adjusted by adding eigenvalues to function output
  add_pcoord <- function (D, correction = "none") 
  {
    HP_bin <- which(D$HP > 0, arr.ind = TRUE)
    if (correction == "none") {
      H_PCo <- coordpcoa(D$H, correction)
      P_PCo <- coordpcoa(D$P, correction)
      H_PCoVec <- H_PCo$vectors
      P_PCoVec <- P_PCo$vectors
      D$H_PCo <- H_PCoVec[HP_bin[, 1], ]
      D$P_PCo <- P_PCoVec[HP_bin[, 2], ]
    }
    else {
      H_PCo <- coordpcoa(D$H, correction)
      P_PCo <- coordpcoa(D$P, correction)
      if (H_PCo$note == "There were no negative eigenvalues. No correction was applied") {
        H_PCoVec <- H_PCo$vectors
        D$H_PCo <- H_PCoVec[HP_bin[, 1], ]
      }
      else {
        H_PCoVec <- H_PCo$vectors.cor
        D$H_PCo <- H_PCoVec[HP_bin[, 1], ]
      }
      if (P_PCo$note == "There were no negative eigenvalues. No correction was applied") {
        P_PCoVec <- P_PCo$vectors
        D$P_PCo <- P_PCoVec[HP_bin[, 2], ]
      }
      else {
        P_PCoVec <- P_PCo$vectors.cor
        D$P_PCo <- P_PCoVec[HP_bin[, 2], ]
      }
    }
    D$H_eg <- H_PCo$values$Rel_corr_eig
    D$P_eg <- P_PCo$values$Rel_corr_eig
    D$correction <- correction
    D$note <- list(H_note = H_PCo$note, P_note = P_PCo$note)
    return(D)
  }
  
  # transform distance to MDS/PCoA and add to class
  D <- add_pcoord(D, correction ="cailliez")
  
  # perform procrustes on the dataset
  D <- PACo(D, nperm=perm, seed=5, method="backtrack", symmetric=TRUE, shuffled=TRUE)
  #D <- PACo_test(D, nperm=perm, seed=5, method="backtrack", symmetric=TRUE, shuffled=TRUE)
  
  get_results <- function(D){
    # results
    cat("The observed m2 is:", D$gof$ss,
        "\nMean m2 permutation:", mean(D$shuffled),
        "\nP-value =",D$gof$p, "based on", D$gof$n,"permutations.")
  }
  # obtain results
  get_results(D)
  
  # create result table
  table = matrix(c(D$gof$ss, mean(D$shuffled), perm, D$gof$p),
                 nrow=4, ncol=1, dimnames = list(c('m2','mean m2.perm','perms','p-value'), 
                 c(paste0(nameH,"_",nameP))))
  
  # add pairwise m2 score
  D$pwc <- c(nameH, nameP, D$proc$ss)
  
  # plot results
  if (plot == T){
    # open a pdf to fill with results
    pdf(paste0(out, "/",nameH,"_",nameP,"_",perm,".pdf"), paper="a4")
    
    # create the trees based on 
    par(mfrow=c(4,2))
    
    # Put the labels at the same height: hang = -1
    plot(hc_H, hang = -1, cex = 0.2, labels= FALSE, main = paste0("Distance_Tree ", nameH))
    plot(hc_P, hang = -1, cex = 0.2, labels= FALSE, main = paste0("Distance_Tree ", nameP))
    
    # plot MDS plots
    plot(D$H_PCo, main=paste0("MDS ", nameH), xlab=round(D$H_eg[1], 2), ylab=round(D$H_eg[2], 3))
    plot(D$P_PCo, main=paste0("MDS ", nameP), xlab=round(D$P_eg[1], 2), ylab=round(D$P_eg[2], 3))
    
    hist(D$shuffled, xlim = c(0.2,1), main = 'Histogram of m2')
    abline(v=D$gof$ss,col="blue")
    abline(v=mean(D$shuffled), col="red")
    legend(x = "topleft", # location of legend within plot area
           c("m2","mean m2.perm"),
           col = c("blue","red"),
           lwd = c(2, 2, 2))
    
    #D <- paco_links(D)
    #plot(D$proc)
    #plot(D$proc, kind=2)
    
    library(plotrix)
    plot.new()
    addtable2plot(0,0,table,
                  display.rownames = TRUE)
    dev.off()
  }
  
  return(c(table, D))
}
