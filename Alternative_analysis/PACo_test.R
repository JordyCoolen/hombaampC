PACo_test <- function (D, nperm = 1000, seed = NA, method = "r0", symmetric = FALSE, 
          proc.warnings = TRUE, shuffled = FALSE) 
{
  correction <- D$correction
  method <- match.arg(method, c("r0", "r1", "r2", "r00", "c0", 
                                "swap", "tswap", "backtrack", "quasiswap"))
  if (!("H_PCo" %in% names(D))) 
    D <- add_pcoord(D, correction = correction)
  if (proc.warnings == TRUE) {
    proc <- vegan::procrustes(X = D$H_PCo, Y = D$P_PCo, symmetric = symmetric)
  }
  else {
    proc <- suppressWarnings(vegan::procrustes(X = D$H_PCo, 
                                               Y = D$P_PCo, symmetric = symmetric))
  }
  Nlinks <- sum(D$HP)
  m2ss <- proc$ss
  pvalue <- 0
  if (!is.na(seed)) 
    set.seed(seed)
  null_model <- vegan::nullmodel(D$HP, method)
  randomised_matrices <- stats::simulate(null_model, nsim = nperm)
  if (shuffled == TRUE) {
    rands <- c()
    for (n in c(1:nperm)) {
      permuted_HP <- randomised_matrices[, , n]
      permuted_HP <- permuted_HP[rownames(D$HP), colnames(D$HP)]
      pdf(paste0('test_', n,'.pdf'))
      heatmap(permuted_HP)
      dev.off()
      
      perm_D <- list(H = D$H, P = D$P, HP = permuted_HP)
      perm_paco <- add_pcoord(perm_D, correction = correction)
      if (proc.warnings == TRUE) {
        perm_proc_ss <- vegan::procrustes(X = perm_paco$H_PCo, 
                                          Y = perm_paco$P_PCo, symmetric = symmetric)$ss
      }
      else {
        perm_proc_ss <- suppressWarnings(vegan::procrustes(X = perm_paco$H_PCo, 
                                                           Y = perm_paco$P_PCo, symmetric = symmetric)$ss)
      }
      if (perm_proc_ss <= m2ss) 
        pvalue <- pvalue + 1
      rands <- c(rands, perm_proc_ss)
    }
  }
  else {
    for (n in c(1:nperm)) {
      permuted_HP <- randomised_matrices[, , n]
      permuted_HP <- permuted_HP[rownames(D$HP), colnames(D$HP)]
      perm_D <- list(H = D$H, P = D$P, HP = permuted_HP)
      perm_paco <- add_pcoord(perm_D, correction = correction)
      if (proc.warnings == TRUE) {
        perm_proc_ss <- vegan::procrustes(X = perm_paco$H_PCo, 
                                          Y = perm_paco$P_PCo, symmetric = symmetric)$ss
      }
      else {
        perm_proc_ss <- suppressWarnings(vegan::procrustes(X = perm_paco$H_PCo, 
                                                           Y = perm_paco$P_PCo, symmetric = symmetric)$ss)
      }
      if (perm_proc_ss <= m2ss) 
        pvalue <- pvalue + 1
    }
  }
  pvalue <- pvalue/nperm
  D$proc <- proc
  D$gof <- list(p = pvalue, ss = m2ss, n = nperm)
  D$method <- method
  D$symmetric <- symmetric
  D$correction <- correction
  if (exists("rands") == TRUE) 
    D$shuffled <- rands
  return(D)
}
