### 今いるフォルダのチェックと移動 ###
getwd()
rm(list = ls(all.names = T))
## 動く先のPATHを指定 ##
## 動く先がわからない人はTABで指定 ##
## 動かない硬派な人はそのままでOK ##
setwd("C:/Users/240909_handson/")

### 使うパッケージのインストール ###
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("genefilter","factoextra",
                       "circlize", "TCC", "ggplot2", "msigdbr",
                       "ConsensusClusterPlus", "ComplexHeatmap", "biomaRt",
                       "clusterProfiler", "enrichplot", "ggnewscale",
                       "GSEABase","GSVA"),force = TRUE)

## 入れたパッケージが動くか確認 ##
library(ggplot2)
library(TCC)
library(genefilter)
library(factoextra)
library(circlize)
library(ComplexHeatmap)
library(biomaRt)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(GSVA)
library(GSEABase)

### 使うデータのダウンロード ###
## MsigDB ##
MsigDB = "https://drive.usercontent.google.com/download?id=1qSr6Jmj0ROCCGQgFtdHguiWAVL_hPFVi&confirm=xxx"
o1 = "msigdb.v2024.1.Hs.symbols.gmt"
download.file(url = MsigDB, destfile = o1, method = "curl")
## Master ##
Master_Data = "https://drive.usercontent.google.com/download?id=1IQEqdkmW9z-v-wXuozHcghmlS6imagFV&confirm=xxx"
o2 = "TCGA_PRAD_100_Data_master.txt"
download.file(url = Master_Data, destfile = o2, method = "curl")
## Extracted Data ##
EX_PRAD_Data = "https://drive.usercontent.google.com/download?id=1xfzk1To23Laq0hku3ZkOD8PEnKzvPEhI&confirm=xxx"
o3 = "TCGA_PRAD_100_Data_master.txt"
download.file(url = Master_Data, destfile = o3, method = "curl")
## Sample Info ##
Sample_data = "https://drive.usercontent.google.com/download?id=1n20_KxvuXA-dd1BkDNSW7tXjlQGvDtR9&confirm=xxx"
o4 = "Sample_Info.txt"
download.file(url = Sample_data, destfile = o4, method = "curl")

