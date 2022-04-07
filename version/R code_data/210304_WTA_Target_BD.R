#background subtration of BD Rhapsody WTA data based on bimodal distribution of lognormal raw read counts 
#of each gene by using mclust package
#Written by Shigeyuki Shichino, ver0.2 20190603
#read raw count table (gene x cell BC table, number of cell BC is 2x of inflection threshold count)

#ライブラリの読み込み
source('library_source_Seurat.R') #first load reticulate python extension settings
library(rDBEC)
library(Seurat)
library(stringr)
library(viridis)
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
#BD Rhapsody WTA analysis with Seurat package
set.seed(seed = 42)

#データ名の指定
input_name <- "Aoki_Seurat.rda"

#クラスタリングの解像度の指定
resol <- 0.6

#出力名の指定

dir.name=("Seurat_plots/")
sample.name=("BDtarget_Seurat_Aoki")

#=========================================================
#元データの読み込み
load(input_name)

#################クラスタリングの実行#######################
#クラスタリングに使うPCの数を指定 (dims.use; 遺伝子セットが濃縮されているPC)
#pvalueの閾値を満たすPCを自動で抽出
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
dims= c(1:(min(tmp1)-1))

#クラスターをどこまで分けるかはマニュアル (resolution)
mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution =resol, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

#####################################################################################
#tSNE結果の書き出し (スーラ)
#pt.size; 細胞あたりのポイントの大きさ (細胞数に応じて調整)
#group.by; メタデータのどのパラメタで色分けするか
#head(mBC@meta.data)でメタデータにアクセス可能 (res.X; クラスタリングのresolutionの値)
#デフォルトでやると最後のクラスタリング結果によって色分けされる
p1 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 0.3, no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 
p2 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = FALSE, 
             label.size = 10,group.by = "TagIDs",
             do.return = TRUE, vector.friendly = TRUE, pt.size = 0.5, no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

#cowplotでlegendだけ抜き出し, save_plotで並べ直す (legendの文字数が長い時用)
legend1 <- cowplot::get_legend(p1)
legend2 <- cowplot::get_legend(p2)
p1 = p1 + theme(legend.position = 'none')
p2 = p2 + theme(legend.position = 'none')
file.name=paste(dir.name, sample.name, "FItSNE_reso", resol, ".png", sep='')
save_plot(file = file.name, plot_grid(p1, legend1, p2, legend2, ncol=2, nrow=2), device="png", 
          units="in", dpi = 600, base_width = 10, base_height = 10, limitsize=FALSE)

#####################################################################################
#スーラオブジェクトの出力
file.name=paste(sample.name, "_res", resol, ".rda", sep='')
save(mBC, file=file.name)

#####################################################################################
#tSNE結果に各細胞の由来をハイライト (orig.ident)
tmp = sort(as.vector(unique(mBC@meta.data$TagIDs)))
file.name=paste("HighlightSampleOrigin", sep='')
dir.create(file.name)

for (i in c(1:length(tmp))) {
  file.name = sprintf("./HighlightSampleOrigin/%s.png", tmp[i])
  cols = rep("red", length(tmp))
  cols[-i] = color=rgb(0.3,0.3,0.3,alpha=0.1)
  p = DimPlot(object = mBC, group.by = "TagIDs", reduction.use = "FItSNE", 
              do.return = TRUE, pt.size = 1.0, cols.use=cols)
  legend <- cowplot::get_legend(p)
  p = p + theme(legend.position = 'none')
  save_plot(file = file.name, plot_grid(p, legend, ncol=2, nrow=1), device="png", 
            units="in", dpi = 300, base_width = 10, base_height = 5, limitsize=FALSE)
}

#################################################################################################
#マーカー遺伝子の抽出・出力(七野さん自作。並列処理)
#min.pct; あるクラスター内で発現している細胞の割合
#pvalue threshold = 0.05
#nthreads; データの大きさとスペックに依存
#cell_ident=as.factor(mBC@meta.data$res.2)
#names(cell_ident)=rownames(mBC@meta.data)
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 12,
                                               adj.p.val.threshold=0.05)
#前処理
mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)
file.name=paste(sample.name, "ALLmarkers_minpct0.1_Adj_p0.05.txt", sep='')
fwrite(mBC.markers, file.name, row.names=F, col.names=T, sep="\t", quote=F)
file.name=paste(sample.name, "ALLmarkers_minpct0.1_Adj_p0.05.rda", sep='')
save(mBC.markers, file=file.name)

#出力
top10 = mBC.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10 = as.data.frame(top10)
top10  = top10 [!duplicated(top10$gene),]
top10 = top10 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top10 = as.data.frame(top10)
file.name=paste(sample.name, "marker_res", resol, ".png", sep='')
p = DoHeatmap2(object = mBC, genes.use = top10$gene, genes.ident = top10$cluster,
               slim.col.label = TRUE, remove.key = FALSE, cex.row=10, disp.min = -2.5, disp.max = 2.5)
ggsave(file = file.name, plot = p, device="png", units="in", dpi = 300,
       width = 20, height = 25, limitsize=FALSE)

#####################################################################################
#Lineageの確認
tmp_spe = c("Trbc2", "Cd3e", "Cd8a", "Cd4", "Cd14", "Lyz1")
sub_name <- "population"
file.name=paste("MarkerGenePlots_", sub_name, sep='')
dir.create(file.name)

for (i in c(1:length(tmp_spe))) {
  hoge = sprintf("/%s.png", tmp_spe[i])
  file.name = paste("MarkerGenePlots_", sub_name, hoge, sep='')
  png(file.name, width = 512, height = 512)
  FeaturePlot(object = mBC, features.plot = tmp_spe[i], cols.use = c("grey", "red"),
              reduction.use = "FItSNE", no.legend = FALSE,  pt.size = 0.3)
  dev.off()
}

sub_name <- "population"
file.name=paste("VlnPlots_population", "_res", resol, sep='')
dir.create(file.name)

for (i in 1:length(tmp_spe)) {
  hoge = sprintf("/%s.png", tmp_spe[i])
  file.name = paste("VlnPlots_", sub_name, "_res", resol, hoge, sep='')
  png(file.name, width = 768, height = 512)
  p <- VlnPlot(object = mBC, features.plot = tmp_spe[i], point.size.use = 0.1)
  plot(p)
  dev.off()
}

png("population_marker.png", width = 1536, height = 512)
p <- DoHeatmap(object = mBC, genes.use = tmp_spe, 
               slim.col.label = TRUE, group.label.rot = TRUE)
plot(p)
dev.off()

#####################################Cell BCの出力##################################
#T cell関連遺伝子の転写量も併せて出力

#raw.dataの取り出し
raw.data <- mBC@raw.data
raw.data.mat <- t(as.matrix(raw.data))

#遺伝子発現量をmetadataに付記
mBC = AddMetaData(object = mBC, metadata = log2(1+raw.data.mat[,"Trbc2"]), col.name = "Trbc2")
mBC = AddMetaData(object = mBC, metadata = log2(1+raw.data.mat[,"Trbc1"]), col.name = "Trbc1")
mBC = AddMetaData(object = mBC, metadata = log2(1+raw.data.mat[,"Trac"]), col.name = "Trac")
mBC = AddMetaData(object = mBC, metadata = log2(1+raw.data.mat[,"Cd3e"]), col.name = "Cd3e")
mBC = AddMetaData(object = mBC, metadata = log2(1+raw.data.mat[,"Cd4"]), col.name = "Cd4")
mBC = AddMetaData(object = mBC, metadata = log2(1+raw.data.mat[,"Cd8a"]), col.name = "Cd8a")

meta.data <- as.data.frame(mBC@meta.data)
write.csv(meta.data, "meta.data.csv", row.names = TRUE)

