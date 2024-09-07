### 今いるフォルダのチェックと移動 ###
getwd()
rm(list = ls(all.names = T))
## 動く先のPATHを指定 ##
## 動く先がわからない人はTABで指定 ##
## 動かない硬派な人はそのままでOK ##
setwd("C:/Users/240909_handson/")

### 必要なファイルがあるか確認 ###
list.files("./", all.files=T)

### DEGの抽出（お試し） ###
f0 = "TCGA_PRAD_100_Data_master.txt"
f1 = "TCGA_PRAD_Ex_Data.txt"
f2 = "Sample_Info.txt"
o0 = "TCGA_PRAD_DEG.txt"

## データの読み込み ##
master = read.table(file = f0, header = TRUE, sep = "\t", row.names = 1)

## サンプル情報の確認 ##
hoge = read.table(file = f2, header = TRUE, sep = "\t")
sub1 = subset(x = hoge, subset = hoge$Sample_Type == "Solid_Tissue_Normal")
normal_name = sub1$Sample_ID
sub2 = subset(x = hoge, subset = hoge$Sample_Type == "Primary_Tumor")
tumor_name = sub2$Sample_ID
sample_list = c(normal_name, tumor_name)
data = master[,sample_list]

## Normalizing ##
## このデータはNormalize済みなのでやってはいけない行為 ##
library(TCC)
param_G1 = length(normal_name)
param_G2 = length(tumor_name)
param_FDR = 0.05
data.cl = c(rep(1, param_G1), rep(2, param_G2))
tcc = new("TCC", data, data.cl)
tcc = calcNormFactors(tcc, norm.method="tmm", test.method="edger",#正規化を実行した結果をtccに格納
                      iteration=3, FDR=0.1, floorPDEG=0.05)
normalized = getNormalizedData(tcc)   #正規化後のデータを取り出してnormalizedに格納

## DEGの抽出 ##
tcc = estimateDE(tcc, test.method="edger", FDR=param_FDR)#DEG検出を実行した結果をtccに格納
result = getResult(tcc, sort=FALSE)   #p値などの計算結果をresultに格納

## results ##
tmp = cbind(rownames(tcc$count), normalized, result)#正規化後のデータの右側にDEG検出結果を結合したものをtmpに格納
tmp = tmp[order(tmp$rank),]           #発現変動順にソートした結果をtmpに格納
write.table(tmp, o0, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存

Ex_data = subset(x = tmp, subset = tmp$estimatedDEG == 1)
write.table(Ex_data, f1, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存


