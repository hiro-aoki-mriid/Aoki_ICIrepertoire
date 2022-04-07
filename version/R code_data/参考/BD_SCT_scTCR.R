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
input_name <- "Scored_Sub1_BDtarget_Seurat_Aoki_res0.6scTCRmerged.OL.rda"
sample.name <- "res0.6scTCRmerged.OL"
resol <- 0.6

#データの読み込み
load(input_name)
dir.name <- "."

##########################################################
#meta.dataの呼び出し
meta.data <- as.data.frame(mBC@meta.data)

#########################################################################################
#Oligoclonal / Polyclonal / Non-overlapping cloneのプロット

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
  file.name = paste(dir.name, file.name, sep='/')
  save_plot(file = file.name, plot(p1), device="tiff", 
            units="in", dpi = 600, base_width = 1.8, base_height = 1.8, limitsize=FALSE)
}

data$freq_tumor <- data$temp_count / sum(data$temp_count)

mBC@meta.data <- meta.data
write.csv(meta.data, "meta.data.afterscTCR.csv", row.names = TRUE)
write.csv(data, "clone_within_cluster.OLclass.csv", row.names = FALSE)

#スーラオブジェクトの出力
file.name=paste(sample.name, "scTCRmerged.OL", ".rda", sep='')
save(mBC, file=file.name)

###################################################################################
#Polyclonal fractionとOligoclonal Fractionの間のDEGの抽出

ident <- mBC@ident 
names.ident <- names(ident)
ident <- as.factor(mBC@meta.data$class)
names(ident) <- names.ident
mBC@ident <- ident

DEG.subset.1 = FindMarkers(mBC, ident.1 = "Oligo", ident.2 = "Poly",
                           test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                           features.use = NULL, nthreads = 12,
                           adj.p.val.threshold=0.05)
DEG.subset.2 = FindMarkers(mBC, ident.1 = "Poly", ident.2 = "Oligo",
                           test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                           features.use = NULL, nthreads = 12,
                           adj.p.val.threshold=0.05)

write.csv(DEG.subset.1, "Oligo.Up.csv", row.names = TRUE)
write.csv(DEG.subset.2, "Poly.Up.csv", row.names = TRUE)

#DotPlotによる評価
mBC2 <- SubsetData(mBC, ident.use = c("Poly", "Oligo"))
genes <- c("Tnfrsf4", "Fas", "Il7r", "Bcl6", "Xcl1", "Ccr7", "Tcf7",
           "Entpd1", "Havcr2", "Gzmb")

p <- DotPlot(mBC2, genes, cols.use = c("blue", "red"),
        dot.min = 0, dot.scale = 4,
        scale.by = "size", scale.min = NA, scale.max = NA, group.by=ident,
        plot.legend = TRUE, do.return = TRUE) +
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(size=8, colour = 1, family = "Arial"),
        axis.text.x = element_text(size=8, colour = 1, family = "Arial", angle = 90)) +
  theme(panel.border = element_rect(fill = NA, size = 0.1),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size=6, colour = 1, family = "Arial"),
        legend.title = element_text(size=6, colour = 1, family = "Arial"),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
file.name <- paste("DEG.OligovsPoly", "tiff", sep=".")
save_plot(file = file.name, plot(p), device="tiff", 
          units="in", dpi = 600, base_width = 3.5, base_height = 1.6, limitsize=FALSE)

