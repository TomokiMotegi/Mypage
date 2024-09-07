### 今いるフォルダのチェックと移動 ###
getwd()
rm(list = ls(all.names = T))
## 動く先のPATHを指定 ##
## 動く先がわからない人はTABで指定 ##
## 動かない硬派な人はそのままでOK ##
setwd("C:/Users/240909_handson/")

### 必要なファイルがあるか確認 ###
list.files("./", all.files=T)

### データのセット ###
f1 = "TCGA_PRAD_Ex_Data.txt"
f2 = "Sample_Info.txt"
o1 = "hc_dist_hclust.png"
o2 = "hc_asdist_hclust_Ave.png"
o3 = "hc_asdist_hclust_WD2.png"
o4 = "hc_TCC_Ave.png"
o5 = "hc_TCC_WD2.png"
o6 = "heatmap_full.png"
o7 = "heatmap_limited_range.png"
o8 = "heatmap_final.png"
o9 = "GSVA.txt"


### Hierarchical Clustering ###
master = read.table(file = f1, header = TRUE, sep = "\t", row.names = 1)
end_data = ncol(master)-7
data = lapply(master[,1:end_data], as.numeric)
data = master[,1:end_data]
data[1,]

## 古典的手法によるHierarchical Clustering ##
rd = dist(t(data), method="euclidean")
rc = hclust(d = rd, method="ward.D2")
png(filename = o1, width = 1400, height = 400)#出力ファイルの各種パラメータを指定
par(mar=c(0, 4, 1, 0)) #下、左、上、右の順で余白（行）を指定
plot(rc, sub = "", xlab = "", cex.lab = 1, 
     cex = 1, main="", ylab="Height")
dev.off()

## RNA-seqの0カウントによる影響を排除して平均距離を計算 ##
rd = dist(1 - cor(data, method="spearman"))
rc = hclust(d = rd, method="average")
png(filename = o2, width = 1400, height = 400)#出力ファイルの各種パラメータを指定
par(mar=c(0, 4, 1, 0)) #下、左、上、右の順で余白（行）を指定
plot(rc, sub = "", xlab = "", cex.lab = 1, 
     cex = 1, main="", ylab="Height")
dev.off()

## RNA-seqの0カウントによる影響を排除してWard法で計算 ##
rd = dist(1 - cor(data, method="spearman"))
rc = hclust(d = rd, method="ward.D2")
png(filename = o3, width = 1400, height = 400)#出力ファイルの各種パラメータを指定
par(mar=c(0, 4, 1, 0)) #下、左、上、右の順で余白（行）を指定
plot(rc, sub = "", xlab = "", cex.lab = 1, 
     cex = 1, main="", ylab="Height")
dev.off()


## TCCパッケージでHCを実施して平均距離で計算 ##
library(TCC)
out = clusterSample(data, dist.method="spearman",
                    hclust.method="average", unique.pattern=TRUE)
png(filename = o4, width = 1400, height = 400)#出力ファイルの各種パラメータを指定
par(mar=c(0, 4, 1, 0)) #下、左、上、右の順で余白（行）を指定
plot(out, sub = "", xlab = "", cex.lab = 1, 
     cex = 1, main="", ylab="Height")
dev.off()

## TCCパッケージでHCを実施してWard法で距離計算 ##
out = clusterSample(data, dist.method="spearman",
                    hclust.method="ward.D2", unique.pattern=TRUE)
png(filename = o5, width = 1400, height = 400)#出力ファイルの各種パラメータを指定
par(mar=c(0, 4, 1, 0)) #下、左、上、右の順で余白（行）を指定
plot(out, sub = "", xlab = "", cex.lab = 1, 
     cex = 1, main="", ylab="Height")
dev.off()



### k-mean Clustering ###
## ConsensusClsterPlus to samples ##
library(ConsensusClusterPlus)
rcc = ConsensusClusterPlus(as.matrix(data), maxK = 10, reps = 100, 
                           pItem = 0.8, pFeature = 1, seed = 2024,
                           distance = "euclidean", clusterAlg = "km",
                           innerLinkage="ward.D2", finalLinkage="ward.D2",
                           plot = "png", title = "CCP_kmeans")
nc = 4
clust.col = rcc[[nc]]$consensusTree


### ヒートマップの作製 ##
## Zスコア化 ##
library(genefilter)
exccp_z = genescale(data, axis = 1, method = "Z")
library(ggplot2)
library(ComplexHeatmap)
set.seed(2024)
gg = ComplexHeatmap::Heatmap(matrix = exccp_z,
                             column_title = "TCGA_PRAD_DEGs",
                             cluster_columns = clust.col,
                             column_split = nc,
                             cluster_row_slices = TRUE,
                             show_row_names = FALSE,
                             clustering_distance_rows = "euclidean",
                             clustering_method_rows = "ward.D2",
                             column_names_side = "bottom",
                             row_names_gp = gpar(fontsize=7),
                             column_names_gp = gpar(fontsize=7),
                             heatmap_legend_param = list(
                               title = "Z-Socred Central Coverage",
                               title_position = "leftcenter-rot",
                               title_gp = gpar(fontsize = 8, 
                                               fontface = "bold"),
                               labels_gp = gpar(fontsize = 8)
                             ),
                             row_dend_reorder = FALSE)
png(file = o6, width = 10, height = 4, units="in", res=600)
set.seed(2024)
draw(gg)
dev.off()


## Zスコアのレンジ制限 ##
exccp_z[exccp_z < -2] = -2
exccp_z[exccp_z > 2] = 2
gg = ComplexHeatmap::Heatmap(matrix = exccp_z,
                             column_title = "TCGA_PRAD_DEGs",
                             cluster_columns = clust.col,
                             column_split = nc,
                             show_row_names = FALSE,
                             clustering_distance_rows = "euclidean",
                             clustering_method_rows = "ward.D2",
                             column_names_side = "bottom",
                             row_names_gp = gpar(fontsize=7),
                             column_names_gp = gpar(fontsize=7),
                             heatmap_legend_param = list(
                               title = "Z-Socred Central Coverage",
                               title_position = "leftcenter-rot",
                               title_gp = gpar(fontsize = 8, 
                                               fontface = "bold"),
                               labels_gp = gpar(fontsize = 8)
                             ),
                             row_dend_reorder = FALSE)
png(file = o7, width = 10, height = 4, units="in", res=600)
set.seed(2024)
draw(gg)
dev.off()


## 多次元データの一括描画 ##
tissue_info = read.table(file = f2, header = TRUE, sep = "\t", row.names = 1)
sample_name = rownames(tissue_info)
cluster_info = as.data.frame(rcc[[nc]]$consensusClass)
library(dplyr)
sample_info = merge(cluster_info, tissue_info, by = "row.names", 
                    all = TRUE, sort = FALSE)
colnames(sample_info) = c("Sample_ID", "Cluster", "Sample_Type")
rownames(sample_info) = sample_info$Sample_ID

# 各パラメータの色を定義 #
library(circlize)
ann_colors = list(
  Sample_Type = c(Primary_Tumor = "#f6aa00", Solid_Tissue_Normal = "#666666"),
  Cluster = c("1" = "#FFFF99", "2" = "#ff8082", "3" = "#ff4b00", "4" = "#804000")
)

# 各パラメータの情報を設定 #
library(ggplot2)
library(ComplexHeatmap)
Annotation = HeatmapAnnotation(Cluster = sample_info$Cluster,
                               Sample_Type = sample_info$Sample_Type,
                               show_annotation_name = TRUE,
                               annotation_name_side = "right",
                               annotation_name_gp = gpar(fontsize=8),
                               annotation_legend_param = list(
                                 title_position= "topleft",
                                 title_gp = gpar(fontsize = 8, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 8)),
                               simple_anno_size = unit(2,"mm"),
                               col = ann_colors)

# ヒートマップの描画 ##
set.seed(2024)
gg = ComplexHeatmap::Heatmap(matrix = exccp_z,
                             column_title = "TCGA_PRAD_DEGs",
                             cluster_columns = clust.col,
                             column_split = nc,
                             cluster_row_slices = TRUE,
                             row_split = 4,
                             show_row_names = FALSE,
                             clustering_distance_rows = "euclidean",
                             clustering_method_rows = "ward.D2",
                             column_names_side = "bottom",
                             row_names_gp = gpar(fontsize=7),
                             column_names_gp = gpar(fontsize=7),
                             col =  circlize::colorRamp2(c(-2, 0, 2), 
                                                         c("blue", "white", "red")),
                             heatmap_legend_param = list(
                               at = c(-2, 2),
                               labels = c("Low(-2)","High(2)"),
                               title = "Z-Socred TPM",
                               title_position = "leftcenter-rot",
                               title_gp = gpar(fontsize = 8, 
                                               fontface = "bold"),
                               labels_gp = gpar(fontsize = 8)
                             ),
                             top_annotation = Annotation,
                             row_dend_reorder = TRUE,
                             column_dend_reorder = FALSE)
png(file = o8, width = 10, height = 6, units="in", res=600)
set.seed(2024)
draw(gg)
dev.off()



### エンリッチメント解析 ###
## prepare module data ##
set.seed(2024)
row_hc = kmeans(exccp_z, centers = 4)
module_number = 3
Module_gene_list = names(row_hc$cluster[row_hc$cluster == module_number])

## Ensembl ID to Universal gene symbols ##
library(biomaRt)
mart = useMart(biomart="ensembl",
               dataset = "hsapiens_gene_ensembl",
               host = "https://useast.ensembl.org")
a = listAttributes(mart)
gene_name_data = getBM(
  mart = mart,
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filter = "ensembl_gene_id",
  values = Module_gene_list,
  uniqueRows = TRUE)
input_gene_list = unique(gene_name_data$external_gene_name)
any(duplicated(input_gene_list)) #遺伝子名の重複を確認


## EnricheR ##
library(msigdbr)
m_df = msigdbr(species = "Homo sapiens")

msig_category = "C5"
m_t2g = msigdbr(species = "Homo sapiens", category = msig_category) %>% 
  dplyr::select(gs_name, gene_symbol)

library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(ggnewscale)
enrich_res = enricher(input_gene_list, TERM2GENE=m_t2g, 
                      minGSSize = 5, # Minimum gene set size
                      maxGSSize = 500, # Maximum gene set set
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")


if(!dir.exists(paste0("Module", module_number, "_", msig_category))){
  dir.create(paste0("Module", module_number, "_", msig_category))
}

write.table(enrich_res@result, 
            file = paste0("Module", module_number, "_", msig_category,"/result.txt"), 
            quote = FALSE, append = FALSE, sep = "\t")

ggsave(filename = paste0("Module", module_number, "_", msig_category,"/dotplot.png"), 
       plot = dotplot(enrich_res), 
       dpi = 300, width = 6, height = 8)

ggsave(filename = paste0("Module", module_number, "_", msig_category,"/qvalute_dotplot.png"), 
       dotplot(enrich_res, x="p.adjust"), 
       dpi = 300, width = 6, height = 8)





### GSEA ###
### prepare module data ###
set.seed(2024)
row_hc = kmeans(exccp_z, centers = 4)
module_number = 1
Module_gene_list = names(row_hc$cluster[row_hc$cluster == module_number])
Module_data = master[Module_gene_list,]

library(biomaRt)
mart = useMart(biomart="ensembl",
               dataset = "hsapiens_gene_ensembl",
               host = "https://useast.ensembl.org")
a = listAttributes(mart)
gene_name_data = getBM(
  mart = mart,
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filter = "ensembl_gene_id",
  values = Module_gene_list,
  uniqueRows = TRUE)
rownames(gene_name_data) = gene_name_data$ensembl_gene_id
any(duplicated(gene_name_data$external_gene_name))

Ex_data = merge(gene_name_data, Module_data, by = "row.names")
mapped_df = subset(x = Ex_data, subset = external_gene_name != "")
dup_gene_symbols <- mapped_df %>%
  dplyr::filter(duplicated(external_gene_name)) %>%
  dplyr::pull(external_gene_name)

mapped_df %>%
  dplyr::filter(external_gene_name %in% dup_gene_symbols) %>%
  dplyr::arrange(external_gene_name)

filtered_mapped_df <- mapped_df %>%
  dplyr::arrange(dplyr::desc(abs(m.value))) %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

# recheck_duplication
any(duplicated(filtered_mapped_df$external_gene_name))

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- filtered_mapped_df$m.value
names(lfc_vector) <- filtered_mapped_df$external_gene_name

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
# Look at first entries of the ranked log2 fold change vector
head(lfc_vector)

library(msigdbr)
m_df = msigdbr(species = "Homo sapiens")

msig_category = "C5"
m_t2g = msigdbr(species = "Homo sapiens", category = msig_category) %>% 
  dplyr::select(gs_name, gene_symbol)

library(clusterProfiler)
# GESA_Set_seed
gsea_results = GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = 2024, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = m_t2g
)

if(!dir.exists(paste0("GESA_Module", module_number, "_", msig_category))){
  dir.create(paste0("GESA_Module", module_number, "_", msig_category))
}


library(enrichplot)
ggsave(filename = paste0("GESA_Module", module_number, "_", msig_category,"/dotplot.png"), 
       plot = dotplot(gsea_results), 
       dpi = 300, width = 6, height = 8)

ggsave(filename = paste0("GESA_Module", module_number, "_", msig_category,"/qvalute_dotplot.png"), 
       dotplot(gsea_results, x="p.adjust"), 
       dpi = 300, width = 6, height = 8)

# emap_data = pairwise_termsim(gsea_results)
# ggsave(filename = paste0("GESA_Module", module_number, "_", msig_category,"/treeplot.png"),
#        treeplot(gsea_results),
#        dpi = 300, width = 16, height = 9)


# Visualizaion #
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

# install.packages("readr")
library(readr)
readr::write_tsv(
  gsea_result_df, file = paste0("GESA_Module", module_number, "_", msig_category,"/", 
                                msig_category, "_analsis.txt")
)

for (name in gsea_result_df$ID) {
  most_positive_nes_plot <- enrichplot::gseaplot(
    gsea_results,
    geneSetID = name,
    title = name,
    color.line = "#0d76ff"
  )
  png(filename = paste0("GESA_Module", module_number, "_", msig_category,"/", 
                        msig_category, "_", name,".png"), 
      width = 600, height = 800)
  print(most_positive_nes_plot)
  dev.off()
}


### GSVA ###
library(GSVA)
library(GSEABase)

gene_name_data = getBM(
  mart = mart,
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filter = "ensembl_gene_id",
  values = rownames(master),
  uniqueRows = TRUE)
rownames(gene_name_data) = gene_name_data$ensembl_gene_id
any(duplicated(gene_name_data$external_gene_name))

Ex_data = merge(gene_name_data, master, by = "row.names")
mapped_df = subset(x = Ex_data, subset = external_gene_name != "")
dup_gene_symbols <- mapped_df %>%
  dplyr::filter(duplicated(external_gene_name)) %>%
  dplyr::pull(external_gene_name)

mapped_df %>%
  dplyr::filter(external_gene_name %in% dup_gene_symbols) %>%
  dplyr::arrange(external_gene_name)

filtered_mapped_df <- mapped_df %>%
  dplyr::arrange(dplyr::desc(abs(m.value))) %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)
rownames(filtered_mapped_df) = filtered_mapped_df$external_gene_name


sub1 = subset(x = tissue_info, subset = tissue_info$Sample_Type == "Solid_Tissue_Normal")
normal_name = rownames(sub1)
sub2 = subset(x = tissue_info, subset = tissue_info$Sample_Type == "Primary_Tumor")
tumor_name = rownames(sub2)
sample_list = c(normal_name, tumor_name)
Tumor_data = round(filtered_mapped_df[,sample_list],0)

## GSVA required matrix data
## if show ERROR: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
## update R versions
in_f2 = "msigdb.v2024.1.Hs.symbols.gmt"
geneset = getGmt(in_f2, geneIdType=SymbolIdentifier(),#in_f2で指定したファイルの読み込み
                 collectionType=BroadCollection(category="c5"))
library(GSEABase)
library(GSVAdata)

data(c2BroadSets)

# TPMのような連続データではなくカウントデータの場合はkcdf="Poisson"
outparam = gsvaParam(as.matrix(Tumor_data), geneset,
                     minSize=5, maxSize=500, kcdf="Gaussian")
out = gsva(outparam)
dim(out)

#後処理(wilcox.testを実行)
pvalue <- NULL   
for(i in 1:nrow(out)){                 #遺伝子セット数に相当するnrow(out)回だけループを回す
  test_name = rownames(out)[i]
  wc_data_1 = out[i, 1:length(normal_name)]
  wc_data_2 = out[i, (length(normal_name)+1):length(sample_list)]
  test = wilcox.test(wc_data_1, wc_data_2)
  pvalue <- append(pvalue,test$p.value) #p値を計算した結果をpvalueに格納
}

#ファイルに保存
tmp <- cbind(rownames(out), pvalue, out)#保存したい情報をtmpに格納
tmp <- tmp[order(pvalue),]             #pvalue順にソートした結果をtmpに格納
write.table(x = tmp, file = o9 , sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存


## ssGSEA ##
outparam = ssgseaParam(as.matrix(Tumor_data), geneset,
                       minSize=5, maxSize=500)
out = gsva(outparam)
dim(out)


