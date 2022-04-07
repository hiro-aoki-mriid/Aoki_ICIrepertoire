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
input_name <- "Scored_Sub1_BDtarget_Seurat_Aoki_res0.6.rda"

#出力名の指定

dir.name <- "scTCR_analysis"
dir.create(dir.name)
sample.name="Scored_Sub1_BDtarget_Seurat_Aoki_res0.6"


#=========================================================
#元データの読み込み
load(input_name)

##########################################################
#meta.dataへのClone.idの統合

#scTCR, クローン情報の読み込み
clone_id_table <- read.csv("CellTable.paired.TCR.csv", header = TRUE)

#OL対象先, クローン情報の読み込み
target_table <- read.table("BDscTCR_dLN.txt", header = TRUE)
target_table$B_strict <- paste(target_table$cdr3nt, target_table$v, target_table$j, sep="_")

#OL対象先にクローンが存在するかを検索
B.clone.id.list <- unique(as.vector(clone_id_table$B_strict))

name_table <- names(clone_id_table)
clone_id_table_OL <- data.frame()

for(i in 1:length(B.clone.id.list)){
  cloneid <- B.clone.id.list[i]
  table_sub <- clone_id_table %>% dplyr::filter(B_strict == cloneid)
  target_sub <- target_table %>% dplyr::filter(B_strict == cloneid)
  if(nrow(target_sub) > 0){
    table_sub$freq.OL <- target_sub$freq
  } else{
    table_sub$freq.OL <- 0
  }
  clone_id_table_OL <- rbind(clone_id_table_OL, table_sub)
}

###Meta dataへの統合
#Clone.idのSeurat object, metadataへの統合
clone.id.list <- clone_id_table_OL$`clone.id`
names(clone.id.list)=clone_id_table_OL$X
mBC = AddMetaData(object = mBC, metadata = clone.id.list, col.name = "clone.id")

#OL先の頻度情報のSeurat object, metadataへの統合
OL.stat.list <- clone_id_table_OL$`freq.OL`
names(OL.stat.list)=clone_id_table_OL$X
mBC = AddMetaData(object = mBC, metadata = OL.stat.list, col.name = "freq.OL")

#meta.dataの出力
meta.data <- as.data.frame(mBC@meta.data)

############################################################################################
#スーラオブジェクトの出力
file.name=paste(sample.name, "scTCRmerged", ".rda", sep='')
save(mBC, file=file.name)

#############################################################################################
#クローンレベルでの集計
clone_id_table_OL$clone.id.freq <- str_c(clone_id_table_OL$clone.id, clone_id_table_OL$freq, sep = "__")
tmp = clone_id_table_OL %>% group_by(`clone.id.freq`) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::arrange(desc(count))
tmp = as.data.frame(tmp)
tmp[,3]=tmp[,2]/sum(tmp[,2])
colnames(tmp)=c("ntSeq_TRA_TRB_freq", "CloneCount", "CloneFreq")
tmp <- tmp[order(tmp$ntSeq_TRA_TRB_freq),]

tmp = cbind(tmp, orig.ident=rep("CD8T", nrow(tmp)))

#dLN中の頻度情報の抜き出し
tmp1 = as.vector(tmp$ntSeq_TRA_TRB_freq)
tmp1_table <- str_split(tmp1, pattern = "__", simplify = TRUE)
tmp_array <- as.vector(tmp1_table[,1])
freq_array <- as.vector(tmp1_table[,2])
temp_count <- as.vector(tmp$CloneCount)

#TCRクローンごとにGene scoreを集計
meta.data.calc <- dplyr::select(meta.data, c("clone.id", "Cytotoxicity1", "Proliferation1", "Tumor.Prog1", "Tumor.Term1"))
out_score <- aggregate(x=meta.data.calc[c("Cytotoxicity1", "Proliferation1", "Tumor.Prog1", "Tumor.Term1")],
                 by=list(meta.data.calc$clone.id), FUN=mean)

#TCRクローンごとに頻度を集計
tmp_out = table(mBC@ident, mBC@meta.data$clone.id)
tmp_out = t(tmp_out)
file.name=paste(dir.name, "clone_within_cluster.txt", sep='/')
write.table(tmp_out, file.name, row.names=T, col.names=T, sep="\t", quote=F)
tmp_out2 <- read.table(file.name, header = TRUE)
tmp_out2$names <- row.names(tmp_out2)
tmp_out2 <- tmp_out2[order(tmp_out2$names),]

tmp_out3 <- cbind(tmp_out2, freq_array, temp_count, out_score)
write.table(tmp_out3, file.name, row.names=T, col.names=T, sep="\t", quote=F)

#########################################################################################
#Oligoclonal / Polyclonal / Non-overlapping cloneのプロット

#Oligoclonal, Polyclonal, NonOLの分類
#OLクローンにランクを付加
data <- read.table(file.name, header = TRUE)

data_OL <- subset(data, freq_array > 0)
data_OL <- data_OL[order(data_OL$temp_count, decreasing=T),]
data_OL$rank <- c(1:nrow(data_OL))

data_nonOL <- subset(data, freq_array == 0)
data_nonOL$rank <- 0

data <- rbind(data_OL, data_nonOL)

#OL, oligo/polyclonalityに基づく分類
data$class <- "nonOL"
data$class[which(data$rank > 0 & data$rank <= 10)] <- "Oligo"
data$class[which(data$rank > 10)] <- "Poly"
tmp_raw <- as.vector(meta.data$clone.id)

meta.data$class <- "nd"
dir.name <- "Clone_plot"
dir.create(dir.name)
clones_array <- c("nonOL", "Oligo", "Poly")
for(i in 1:length(clones_array)){
  clone_class <- clones_array[i]
  data_sub <- dplyr::filter(data, class == clone_class)
  clone.id <- as.vector(data_sub$Group.1)
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
            units="in", dpi = 600, base_width = 1.1, base_height = 1.1, limitsize=FALSE)
  
  #meta.dataへの反映
  meta.data$class[tmp2] <- clone_class
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

#####################################################################################

#10細胞以上のクローンをtSNE上にプロットする
tmp_raw <- as.vector(meta.data$clone.id)
dir.name <- "Clone_plot"
dir.create(dir.name)
for (i in c(1:length(temp_count))) {
#for (i in c(1:10)) {
  count <- temp_count[i]
  if(count > 10){
    clone.id <- tmp_array[i]
    freq <- sprintf("%2.3f", (as.numeric(freq_array[i])*100))
    hoge = paste(freq, "png", sep = ".")
    file.name <- paste("count", count, hoge, sep=".")
    file.name = paste(dir.name, file.name, sep='/')
    png(file.name, width = 512, height = 512)
    tmp2 = tmp_raw %in% clone.id
    cols=rep(rgb(0.3,0.3,0.3,alpha=0.2), length(tmp2))
    cols[tmp2]="red"
    DimPlot(object = mBC, group.by = "clone.id", 
            cells.highlight = row.names(meta.data)[tmp2],
            sizes.highlight = 1.8,
            reduction.use = "FltSNE", do.return = FALSE,
            pt.size = 1, cols.use=cols, no.legend = TRUE)
    dev.off()
  } else{
  }
}

