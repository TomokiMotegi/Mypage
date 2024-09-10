rm(list = ls(all.names = T))
getwd()

#install_packages
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install(c("affy", "hgu133plus2cdf", 
                       "ggplot2", "limma", "gplots", "TCC"), force = T, update = F)
