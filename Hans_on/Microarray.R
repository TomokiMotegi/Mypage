#move_working_directory
getwd()
rm(list = ls(all.names = T))


#clustering_rma
data = read.table("hoge_rma.txt", header=TRUE, row.names=1, sep="\t", quote="")
dim(data)
data1 = t(data)
data.dist = dist(data1) #距離はデフォルトがeuclidean
png("hc_rma_data1.png", width = 1100, height = 500)
out = hclust(data.dist, method = "ward.D2")
plot(out, cex = 1.0, font = 1)
dev.off()


#clustering_mas5
data = read.table("hoge_mas5.txt", header=TRUE, row.names=1, sep="\t", quote="")
dim(data)
data1 = t(data)
data.dist = dist(data1) #距離はデフォルトがeuclidean
png("hc_mas5_data1.png", width = 1100, height = 500)
out = hclust(data.dist, method = "ward.D2")
plot(out, cex = 1.0, font = 1)
dev.off()


#発現変動解析_rma
in_f = "hoge_rma.txt"           #入力ファイル名を指定してin_fに格納
out_f1 = "Limma_ex_rma.txt"                  #出力ファイル名を指定してout_f1に格納
out_f2 = "Limma_MAplot_rma.png"                  #出力ファイル名を指定してout_f2に格納
param_G1 = 38                         #G1群のサンプル数を指定
param_G2 = 92                         #G2群のサンプル数を指定
param_FDR = 0.1                      #false discovery rate (FDR)閾値を指定
param_fig = c(600, 800)               #ファイル出力時の横幅と縦幅を指定(単位はピクセル)

#必要なパッケージをロード
library(limma)                         #パッケージの読み込み

#入力ファイルの読み込みとラベル情報の作成
data = read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")#in_fで指定したファイルの読み込み
data.cl = c(rep(1, param_G1), rep(2, param_G2))#G1群を1、G2群を2としたベクトルdata.clを作成

#本番
design = model.matrix(~data.cl)       #デザイン行列を作成した結果をdesignに格納
fit = lmFit(data, design)             #モデル構築(ばらつきの程度を見積もっている)
out = eBayes(fit)                     #検定(経験ベイズ)
hoge = topTable(out,coef=colnames(design)[ncol(design)], adjust="BH",number= nrow(data))#発現変動順にソートした結果をhogeに格納
sum(hoge$adj.P.Val < 0.1)             #FDR < 0.1を満たす遺伝子数を表示

#ファイルに保存(テキストファイル)
tmp = cbind(rownames(hoge), data[rownames(hoge),], hoge)#入力データの右側にDEG検出結果を結合したものをtmpに格納
write.table(tmp, out_f1, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存

#ファイルに保存(M-A plot)
M = hoge$logFC                        #M-A plotのM値(y軸の値)に相当するものをMに格納
A = hoge$AveExpr                      #M-A plotのA値(x軸の値)に相当するものをAに格納
png(out_f2, pointsize=20, width=param_fig[1], height=param_fig[2])#出力ファイルの各種パラメータを指定
plot(A, M, xlab="AveExpr", ylab="logFC", cex=0.5, pch=20)#M-A plotを描画
grid(col="gray", lty="dotted")         #指定したパラメータでグリッドを表示
obj = as.logical(hoge$adj.P.Val < param_FDR) #条件を満たすかどうかを判定した結果をobjに格納
points(A[obj], M[obj], col="magenta", cex=0.1, pch=20)#objがTRUEとなる要素のみ指定した色で描画
legend("topright", c(paste("DEG(FDR<", param_FDR, ")", sep=""), "non-DEG"),#凡例を作成している
       col=c("magenta", "black"), pch=20)#凡例を作成している
dev.off()  


#発現変動解析_mas5
in_f = "hoge_mas5.txt"           #入力ファイル名を指定してin_fに格納
out_f1 = "Limma_ex_mas5.txt"                  #出力ファイル名を指定してout_f1に格納
out_f2 = "Limma_MAplot_mas5.png"                  #出力ファイル名を指定してout_f2に格納
param_G1 = 38                         #G1群のサンプル数を指定
param_G2 = 92                         #G2群のサンプル数を指定
param_FDR = 0.1                      #false discovery rate (FDR)閾値を指定
param_fig = c(600, 800)               #ファイル出力時の横幅と縦幅を指定(単位はピクセル)

#必要なパッケージをロード
library(limma)                         #パッケージの読み込み

#入力ファイルの読み込みとラベル情報の作成
data = read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")#in_fで指定したファイルの読み込み
data.cl = c(rep(1, param_G1), rep(2, param_G2))#G1群を1、G2群を2としたベクトルdata.clを作成

#本番
design = model.matrix(~data.cl)       #デザイン行列を作成した結果をdesignに格納
fit = lmFit(data, design)             #モデル構築(ばらつきの程度を見積もっている)
out = eBayes(fit)                     #検定(経験ベイズ)
hoge = topTable(out,coef=colnames(design)[ncol(design)], adjust="BH",number= nrow(data))#発現変動順にソートした結果をhogeに格納
sum(hoge$adj.P.Val < 0.1)             #FDR < 0.1を満たす遺伝子数を表示

#ファイルに保存(テキストファイル)
tmp = cbind(rownames(hoge), data[rownames(hoge),], hoge)#入力データの右側にDEG検出結果を結合したものをtmpに格納
write.table(tmp, out_f1, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存

#ファイルに保存(M-A plot)
M = hoge$logFC                        #M-A plotのM値(y軸の値)に相当するものをMに格納
A = hoge$AveExpr                      #M-A plotのA値(x軸の値)に相当するものをAに格納
png(out_f2, pointsize=20, width=param_fig[1], height=param_fig[2])#出力ファイルの各種パラメータを指定
plot(A, M, xlab="AveExpr", ylab="logFC", cex=0.5, pch=20)#M-A plotを描画
grid(col="gray", lty="dotted")         #指定したパラメータでグリッドを表示
obj = as.logical(hoge$adj.P.Val < param_FDR) #条件を満たすかどうかを判定した結果をobjに格納
points(A[obj], M[obj], col="magenta", cex=0.1, pch=20)#objがTRUEとなる要素のみ指定した色で描画
legend("topright", c(paste("DEG(FDR<", param_FDR, ")", sep=""), "non-DEG"),#凡例を作成している
       col=c("magenta", "black"), pch=20)#凡例を作成している
dev.off()         


# 発現変動が大きい遺伝子のheatmap_rma
in_f = "Limma_ex_rma.txt"      #入力ファイル名を指定してin_fに格納
in_f2 = "GPL570-55999_modify.txt"
out_f = "EX_genename_rma.txt"                   #出力ファイル名を指定してout_fに格納
param_G1 = 38                         #G1群のサンプル数を指定
param_G2 = 92 
param = mean

data = read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")#in_fで指定したファイルの読み込み
sum(data$adj.P.Val < 0.1) 
659*0.05
sum(data$logFC > 1) 
sum(data$logFC < -1) 
data_ex1 = subset(data, adj.P.Val < 0.1 & logFC > 1 | adj.P.Val < 0.1 & logFC < -1)
dim(data_ex1)
data_ex2 = subset(data, adj.P.Val < 0.1 & logFC > 0.5 | adj.P.Val < 0.1 & logFC < -0.5)
dim(data_ex2)
data_ex2 = data_ex2[,1:130]

sym = read.table(in_f2, header=TRUE, row.names=1, sep="\t", quote="")#in_f2で指定したファイルの読み込み
EX_probe_name = c(rownames(data_ex2))
EX_probe_name
IDs = sym[EX_probe_name,]
IDs = as.vector(IDs[,2])              #Gene symbol情報をベクトルに変換し、IDsに格納

names(IDs) = EX_probe_name            #IDsを行名で対応づけられるようにしている
uniqID = unique(IDs)                  #non-redundant ID情報を抽出し、uniqIDに格納
uniqID = uniqID[uniqID != ""]         #uniqIDの中から指定したIDがないものを除く
uniqID = uniqID[!is.na(uniqID)]       #uniqIDの中から指定したIDが"NA"のものを除く
uniqID = uniqID[!is.nan(uniqID)]      #uniqIDの中から指定したIDが"NaN"のものを除く

hoge = t(apply(as.matrix(uniqID), 1, function(i, d = data_ex2, s = IDs, p = param) {#uniqIDを一つずつ処理する
  apply(d[which(s == i), ], 2, p, na.rm = TRUE)#dataの中から現在のgene symbol(i)と同じprobesを全て抽出し、その平均値(mean)を返す
}, data_ex2, IDs, param))          #apply関数でdata,IDs,paramを使えるように代入
rownames(hoge) = uniqID               #non-redundant IDをhogeの行の名前として利用

#ファイルに保存
tmp = cbind(rownames(hoge), hoge)     #指定したIDの列を行列hogeの左端に挿入し、結果をtmpに格納
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存

#heatmapを描く
library(gplots)
data = read.table(out_f, header=TRUE, row.names=1, sep="\t", quote="")#in_fで指定したファイルの読み込み
data = na.omit(data)
hm_data = data.matrix(data) #出てくる結果はd = Euclidean, method = "ward.D2"
pdf("heatmap2_rma_q0.1_F2.pdf", height = 7)
heatmap.2(hm_data, margins = c(6,5), col=redgreen(256), trace = "none", lhei = c(4,8), cexRow = 0.4, cexCol = 0.3)
dev.off()


# 発現変動が大きい遺伝子のheatmap_mas5
in_f = "Limma_ex_mas5.txt"      #入力ファイル名を指定してin_fに格納
in_f2 = "GPL570-55999_modify.txt"
out_f = "EX_genename_mas5.txt"                   #出力ファイル名を指定してout_fに格納
param_G1 = 38                         #G1群のサンプル数を指定
param_G2 = 92 
param = mean

data = read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")#in_fで指定したファイルの読み込み
sum(data$adj.P.Val < 0.1) 
96*0.05
sum(data$logFC > 1) 
sum(data$logFC < -1) 
data_ex2 = subset(data, adj.P.Val < 0.1 & logFC > 1 | adj.P.Val < 0.1 & logFC < -1)
data_ex2 = data_ex2[,1:130]

sym = read.table(in_f2, header=TRUE, row.names=1, sep="\t", quote="")#in_f2で指定したファイルの読み込み
EX_probe_name = c(rownames(data_ex2))
EX_probe_name
IDs = sym[EX_probe_name,]
IDs = as.vector(IDs[,2])              #Gene symbol情報をベクトルに変換し、IDsに格納

names(IDs) = EX_probe_name            #IDsを行名で対応づけられるようにしている
uniqID = unique(IDs)                  #non-redundant ID情報を抽出し、uniqIDに格納
uniqID = uniqID[uniqID != ""]         #uniqIDの中から指定したIDがないものを除く
uniqID = uniqID[!is.na(uniqID)]       #uniqIDの中から指定したIDが"NA"のものを除く
uniqID = uniqID[!is.nan(uniqID)]      #uniqIDの中から指定したIDが"NaN"のものを除く

hoge = t(apply(as.matrix(uniqID), 1, function(i, d = data_ex2, s = IDs, p = param) {#uniqIDを一つずつ処理する
  apply(d[which(s == i), ], 2, p, na.rm = TRUE)#dataの中から現在のgene symbol(i)と同じprobesを全て抽出し、その平均値(mean)を返す
}, data_ex2, IDs, param))          #apply関数でdata,IDs,paramを使えるように代入
rownames(hoge) = uniqID               #non-redundant IDをhogeの行の名前として利用

#ファイルに保存
tmp = cbind(rownames(hoge), hoge)     #指定したIDの列を行列hogeの左端に挿入し、結果をtmpに格納
write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存

#heatmapを描く
library(gplots)
data = read.table(out_f, header=TRUE, row.names=1, sep="\t", quote="")#in_fで指定したファイルの読み込み
data = na.omit(data)
hm_data = data.matrix(data) #出てくる結果はd = Euclidean, method = "ward.D2"
pdf("heatmap2_mas5_q0.1_F2.pdf", height = 7)
heatmap.2(hm_data, margins = c(6,5), col=redgreen(256), trace = "none", lhei = c(4,8), cexRow = 0.6, cexCol = 0.3)
dev.off()
