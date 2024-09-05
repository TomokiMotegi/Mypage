#move_working_directory
getwd()
rm(list = ls(all.names = T))
setwd("./RNA_seq_Training/01_Microarray/NEJM_CEL/")


#rma
out_f = "hoge_rma.txt"  #出力ファイル名を指定してout_fに格納
library(affy)         #パッケージの読み込み
library(hgu133plus2cdf)         #パッケージの読み込み
hoge = ReadAffy()    #*.CELファイルの読み込み
eset = rma(hoge)     #RMAを実行し、結果をesetに保存
write.exprs(eset, file=out_f)#結果を指定したファイル名で保存
rm(list = "eset")


#mas5
out_f = "hoge_mas5.txt"  #出力ファイル名を指定してout_fに格納
library(affy)         #パッケージの読み込み
library(hgu133plus2cdf)
hoge = ReadAffy()    #*.CELファイルの読み込み
eset = mas5(hoge)     #mas5.0を実行し、結果をesetに保存
exprs(eset)[exprs(eset) < 1] = 1
exprs(eset) = log(exprs(eset), 2)
write.exprs(eset, file=out_f)          #結果を指定したファイル名で保存
rm(list = ls(all.names = T))


