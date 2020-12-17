## install devtools, if necessary:
install.packages("devtools", dep=TRUE)
library(devtools)

library(gdal)
install.packages("sf", configure.args = "--with-proj-lib=/usr/local/lib/")

## install treeWAS from github:
install_github("caitiecollins/treeWAS", build_vignettes = TRUE)
library(treeWAS)

#### homoplasyfinder (INDEL)

devtools::install_github("JosephCrispell/homoplasyFinder")
library(homoplasyFinder)


runHomoplasyFinderInJava