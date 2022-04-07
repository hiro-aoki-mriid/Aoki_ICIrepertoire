#background subtration of BD Rhapsody WTA data based on bimodal distribution of lognormal raw read counts 
#of each gene by using mclust package
#Written by Shigeyuki Shichino, ver0.2 20190603
#read raw count table (gene x cell BC table, number of cell BC is 2x of inflection threshold count)

#Input; Scored_Sub1_BDtarget_Seurat_Aoki_res0.6scTCRmerged.OL.rda

#方針
#1. T cell以外のコンタミと思われるクラスターを除去
#2. Pmel由来の細胞を除去

#ライブラリの読み込み
source('library_source_Seurat.R') #first load reticulate python extension settings
library(rDBEC)
library(Seurat)
library(stringr)
library(viridis)
library(extrafont)
loadfonts("win")
library(patchwork)
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
#BD Rhapsody WTA analysis with Seurat package
set.seed(seed = 42)

#データ名の指定
input_name <- "Scored_Sub1_BDtarget_Seurat_Aoki_res0.6scTCRmerged.OL.rda"

#出力名の指定

dir.name=("Seurat_plots/")
sample.name=("Sub1_BDtarget_Seurat_Aoki")
dir.create(dir.name)

#=======================================================================
#データの読み込み
load(input_name)

#####################################################################################
#tSNE結果の書き出し (スーラ)
#pt.size; 細胞あたりのポイントの大きさ (細胞数に応じて調整)
#group.by; メタデータのどのパラメタで色分けするか
#head(mBC@meta.data)でメタデータにアクセス可能 (res.X; クラスタリングのresolutionの値)
#デフォルトでやると最後のクラスタリング結果によって色分けされる
p1 = DimPlot(object = mBC, reduction.use = "FltSNE", do.label = FALSE, label.size = 8,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 0.75, no.legend=TRUE,
             group.by = "res.0.6") +
  theme(axis.title.x = element_text(size=8, family = "Arial"), 
        axis.title.y = element_text(size=8, family = "Arial"), 
        axis.text.x = element_text(size=6, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 6, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 0.75)) 

#cowplotでlegendだけ抜き出し, save_plotで並べ直す (legendの文字数が長い時用)
file.name=paste("Fig3B", "FltSNE", "tiff", sep='.')
save_plot(file = file.name, plot(p1), device="tiff", 
          units="in", dpi = 600, base_width = 2.4, base_height = 2.4, limitsize=FALSE)

#################################################################################################
#Markergene, Heatmapの作成

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features.plot = feature, point.size.use = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(family = "Arial", size = 6, angle = 0, colour = "black", face = "plain"), 
          axis.text.y = element_blank(), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(family = "Arial", size=8, face = "plain"), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1, heights = 1)
  return(p)
}

genes <- c("Mki67", "Isg15",
           "Pdcd1", "Havcr2", 
           "Ifng", "Gzmb", "Tnfrsf9", 
           "Tcf7", "Il7r", "Sell")

file.name2 <- "Fig3C.Vlnplot.Cluster.tiff"
ppi <- 600
tiff(file.name2, width=2*ppi, height=2.4*ppi, res=ppi)
p <- StackedVlnPlot(obj = mBC, features = genes)
plot(p)
dev.off()

###############################################################################
#T cell clusterの解析
dir.name2=paste(dir.name, "FigS2A", sep='')
dir.create(dir.name2)
tmp_spe = c("Pdcd1", "Havcr2", "Tigit", "Tnfrsf9", 
            "Ifng", "Gzmb", "Klrg1", "Cd69",
            "Ccr7", "Tcf7", "Il7r", "Sell", 
            "Cd38", "Mki67", "Mx1", "Isg15")

for (i in c(1:length(tmp_spe))) {
  hoge = sprintf("/%s.tiff", tmp_spe[i])
  file.name = paste(dir.name2, hoge, sep='')
  ppi <- 600
  plist <- FeaturePlot(object = mBC, features.plot = tmp_spe[i],
              cols.use = c("grey", "red"),
              reduction.use = "FltSNE", no.legend = TRUE, no.axes = TRUE,
              pt.size = 2, do.return = TRUE, vector.friendly = TRUE)
  p1 <- plist[[1]] + 
    theme(plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=4, colour = 1, family = "Arial"),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
          panel.border = element_blank()) +
    theme(plot.background = element_rect(fill = "transparent",color = NA))
  save_plot(file = file.name, plot(p1), device="tiff", 
            units="in", dpi = 600, base_width = 0.7, base_height = 0.7, limitsize=FALSE) 
}


############################################################################
#Tex term, Tex prog, cytotoxic and proliferation singnatures 
#Figure 3D

files  <- c("Sig.Tumor.Prog.txt", "Sig.Tumor.Term.txt")

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features.plot = feature, point.size.use = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(family = "Arial",size = 8), 
          axis.ticks.x = element_line(size=0.2), 
          axis.title.y = element_blank(), 
          axis.text.y = element_text(family = "Arial",size = 6), 
          plot.margin = plot.margin ) 
  return(p)
}


#Gene scoreの付加

for(i in files){
  sig_name <- str_remove(i, "Sig.")
  sig_name <- str_remove(sig_name, ".txt")
  feature <- paste(sig_name, "1", sep = "")
  
  file.name2 <- paste(sig_name, "tiff", sep = ".")
  ppi <- 600
  tiff(file.name2, width=1.5*ppi, height=1*ppi, res=ppi)
  p <- modify_vlnplot(obj = mBC, feature = feature, pt.size = 0)
  plot(p)
  dev.off()
}

################################################################
#Figure 3E
#Oligoclonal / Polyclonal / Non-overlapping cloneのプロット

#meta.dataの呼び出し
meta.data <- as.data.frame(mBC@meta.data)

tmp_raw <- as.vector(row.names(meta.data))
clones_array <- c("nonOL", "Oligo", "Poly")
for(i in 1:length(clones_array)){
  clone_class <- clones_array[i]
  data_sub <- dplyr::filter(meta.data, class == clone_class)
  clone.id <- as.vector(row.names(data_sub))
  tmp2 = tmp_raw %in% clone.id
  p1 <- DimPlot(object = mBC, group.by = "res.0.6", do.label = FALSE,
                cells.use = row.names(meta.data)[tmp2],
                reduction.use = "FltSNE", do.return = TRUE, vector.friendly = TRUE,
                pt.size = 2, no.legend = TRUE, no.axes = TRUE) +
    theme(axis.title = element_blank(), 
          axis.text = element_blank()) +
    theme(panel.border = element_rect(fill = NA, size = 0.1)) 
  
  file.name <- paste("clone_class", clone_class, "tiff", sep=".")
  file.name = paste(file.name, sep='/')
  save_plot(file = file.name, plot(p1), device="tiff", 
            units="in", dpi = 600, base_width = 1.5, base_height = 1.5, limitsize=FALSE)
}

#######################################################################
#Figure 3G
#Gene score同士 / gene scoreとclone sizeの間で相関解析を行う
#Input; scTCR_analysis/clone_within_cluster_calculated.txt

# preparing the libraries
library(ggplot2)
library(stringr)
library(dplyr)

file.name  <- "clone_within_cluster.OLclass.csv"
d <- read.csv(file.name, header = TRUE)
d_sub <- dplyr::select(d, c("freq_tumor", "Cytotoxicity1", "Proliferation1", "Tumor.Prog1", "Tumor.Term1", "class"))
d_sub <- dplyr::filter(d_sub, freq_tumor > 0.001)
d_sub <- dplyr::filter(d_sub, class %in% c("Poly", "Oligo"))

axis_array <- c("freq_tumor", "Tumor.Prog1", "Tumor.Term1")

#Tumor中の頻度を%に変換
d_sub$freq_tumor <- d_sub$freq_tumor*100

dir.name <- "Scatter"
dir.create(dir.name)

for(i in 1:(length(axis_array)-1)){
  for(j in (i+1):length(axis_array)){
    #軸名の指定
    x.name <- axis_array[i]
    y.name <- axis_array[j]
    
    #データの抽出
    d_sub2 <- dplyr::select(d_sub, c("class", x.name, y.name))
    names(d_sub2) <- c("class", "x_axis", "y_axis")
    
    #出力名の指定
    ppi <- 600
    name_out <- str_c("x", x.name, "y", y.name, "tiff", sep = ".")
    name_out <- str_c(dir.name, name_out, sep = "/")
    tiff(name_out, width=1.5*ppi, height=1.5*ppi, res=ppi)
    p <- ggplot(d_sub2, aes(x=x_axis, y=y_axis, colour=class)) +
      scale_colour_manual(values=c(Poly="blue", Oligo="red"))+
      geom_point(size=0.8) +
      scale_x_log10() +
      theme_bw(base_size = 6) +
      theme(
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(size=8, family="Arial"),
        axis.text.y = element_text(size=8, family="Arial")) +
      guides(colour=FALSE)
    plot(p)
    dev.off()
  }
}
