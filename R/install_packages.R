# install CRAN packages and load libraries
packagesCRAN <- c("broom", "devtools", "edgeR", "ggpubr", "tidyverse", "openxlsx", "plotly")

for (package in packagesCRAN) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Custom package for RNAseq analysis
if (!requireNamespace("RNAseqswissknife", quietly = TRUE)) {
  devtools::install_gitlab("arnaud.polizzi/rnaseqswissknife", host = "https://forgemia.inra.fr")
}
library(RNAseqSwissknife)

rm(package, packagesCRAN)

