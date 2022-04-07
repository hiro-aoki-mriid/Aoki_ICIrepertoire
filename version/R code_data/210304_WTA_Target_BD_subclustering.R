#background subtration of BD Rhapsody WTA data based on bimodal distribution of lognormal raw read counts 
#of each gene by using mclust package
#Written by Shigeyuki Shichino, ver0.2 20190603
#read raw count table (gene x cell BC table, number of cell BC is 2x of inflection threshold count)

#方針
#1. T cell以外のコンタミと思われるクラスターを除去
#2. Pmel由来の細胞を除去

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
input_name <- "BDtarget_Seurat_Aoki_res0.6.rda"

#再解析するクラスターの指定
cluster_vector <- c(0, 1, 2, 3, 4, 5, 6, 8, 9) #For CD8

#クラスタリングの解像度の指定
resol <- 0.6

#出力名の指定

dir.name=("Seurat_plots/")
sample.name=("Sub1_BDtarget_Seurat_Aoki")
dir.create(dir.name)

#=======================================================================
#データの読み込み
load(input_name)
##BD target panelの遺伝子リストを読み込み
BD_genes <- read.table("BD_genes_plus10.txt", header = TRUE)
BD_genes <- as.vector(BD_genes$Genesymbol)

#subclusterの読み込み, スーラオブジェクトの作成
#T cell以外のクラスターを除外 & Pmel由来の細胞を削除
#BD target panelの遺伝子のみを選択
res1 = mBC@raw.data[, WhichCells(object = mBC, ident = cluster_vector, subset.name = "orig.ident", accept.value = c("Aoki.mouseSampleTag5", "Aoki.mouseSampleTag6"))]
res1 =res1[res1@Dimnames[[1]] %in% BD_genes,]
mBC = CreateSeuratObject(raw.data = res1, min.cells = 5, min.genes = 10)

#スーラの処理（ミトコンドリア, リボソームetc)
colnames(mBC@meta.data)[2]="nReads"
nReads_log = log10(Matrix::colSums(mBC@raw.data))
mBC = AddMetaData(object = mBC, metadata = nReads_log, col.name = "nReads.log")

#plot statistics　（scRNAseqの成績の打ち出し)

file.name=paste(dir.name, sample.name, "nGene.png", sep='')
png(file.name, width = 512, height = 400)
RidgePlot(object = mBC, features.plot = "nGene", group.by="orig.ident", nCol = 1)
dev.off()

file.name=paste(dir.name, sample.name, "nReads_log.png", sep='')
png(file.name, width = 512, height = 400)
RidgePlot(object = mBC, features.plot = "nReads.log", group.by="orig.ident", nCol = 1)
dev.off()
#------------------------------------------------------------------

# normalizing data
mBC = NormalizeData(object = mBC, scale.factor=1000000)
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"),
                 do.par=TRUE, num.cores=8) #細胞ごとのリード数の補正
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)　#変動遺伝子の抽出
length(x = mBC@var.genes)


closeAllConnections()
gc()
# perform PCA, フィルター後はPCs = 50で問題なさそう
mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 50, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)

#JackStrawは処理が重い！！
mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 50, do.par=TRUE, num.cores=8)
file.name=paste(dir.name, sample.name, "Jackstraw.png", sep='')　#PCに濃縮されている遺伝子セットを抽出
png(file.name, width = 750, height = 1000)
mBC = JackStrawPlot(object = mBC, PCs = 1:50)
dev.off()
closeAllConnections()
gc()

#################クラスタリングの実行#######################
#クラスタリングに使うPCの数を指定 (dims.use; 遺伝子セットが濃縮されているPC)
#pvalueの閾値を満たすPCを自動で抽出
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
dims= c(1:(min(tmp1)-1))

#クラスターをどこまで分けるかはマニュアル (resolution)
mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                    resolution =resol, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

#################tSNEの実行#######################
cell_embeddings <- mBC@dr$pca@cell.embeddings[,dims]
n_jobs <- parallel::detectCores()/2
FltSNE_module <- reticulate::import(module = "MulticoreTSNE", delay_load = TRUE)
#pythonを使ってtSNE (rからpythonを呼び出す)
#perplexityとかは自動で最適化してくれる
FltSNE <- FltSNE_module$MulticoreTSNE(n_components = as.integer(2),
                                      perplexity = as.numeric(100.0),
                                      early_exaggeration = as.numeric(12.0),
                                      learning_rate = as.numeric(nrow(cell_embeddings)/12),
                                      n_iter = as.integer(1000),
                                      n_iter_without_progress = as.integer(300),
                                      min_grad_norm = as.numeric(1e-07),
                                      metric = "euclidean",
                                      init = "random",
                                      verbose = as.integer(25),
                                      random_state = as.integer(42),
                                      method = "barnes_hut",
                                      angle = as.numeric(0.5),
                                      #auto_iter = TRUE,
                                      #auto_iter_end = as.integer(5000),
                                      n_jobs = as.integer(n_jobs))

#tSNEで次元圧縮されたデータフレームをスーラオブジェクトに代入
FltSNE_df <- FltSNE$fit_transform(cell_embeddings)
mBC <- PushData(
  object = mBC,
  python_df = FltSNE_df,
  reduction_use = "pca",
  reduction_save = "FltSNE",
  dims_use = dims
)

#####################################################################################
#tSNE結果の書き出し (スーラ)
#pt.size; 細胞あたりのポイントの大きさ (細胞数に応じて調整)
#group.by; メタデータのどのパラメタで色分けするか
#head(mBC@meta.data)でメタデータにアクセス可能 (res.X; クラスタリングのresolutionの値)
#デフォルトでやると最後のクラスタリング結果によって色分けされる
p1 = DimPlot(object = mBC, reduction.use = "FltSNE", do.label = TRUE, label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 0.3, no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"), 
        axis.title.y = element_text(size=10, family = "Arial"), 
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial")) +
  theme(panel.border = element_rect(fill = NA, size = 1.0)) 

#cowplotでlegendだけ抜き出し, save_plotで並べ直す (legendの文字数が長い時用)
legend1 <- cowplot::get_legend(p1)
p1 = p1 + theme(legend.position = 'none')
file.name=paste(dir.name, sample.name, "FltSNE_reso", resol, ".png", sep='')
save_plot(file = file.name, plot_grid(p1, legend1, ncol=2, nrow=1), device="png", 
          units="in", dpi = 600, base_width = 10, base_height = 5, limitsize=FALSE)

#####################################################################################
#スーラオブジェクトの出力
file.name=paste(sample.name, "_res", resol, ".rda", sep='')
save(mBC, file=file.name)

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

###############################################################################
#T cell clusterの解析
dir.name2=paste(dir.name, "/MarkerGenePlots_subset_population", sep='')
dir.create(dir.name2)
tmp_spe = c("Pdcd1", "Tcf7", "Lag3", "Entpd1", "Havcr2", "Tigit")
sub_name <- "Exhausted"
tmp_spe = c("Ccr7", "Il7r", "Bach2", "Sell", "Cxcr5", "Cxcr3")
sub_name <- "Naive_memory"
tmp_spe = c("Cd69", "Gzma", "Gzmb", "Gzmk", "Tnfrsf9", "Klrg1", "Prf1", "Tnf", "Cx3cr1")
sub_name <- "Effector"
tmp_spe = c("Mki67", "Eomes", "Tbx21", "Bcl6", "Maf", "Prdm1")
sub_name <- "Transcription_factor"
tmp_spe = c("Cxcl10", "Cxcl9", "Isg15", "Mx1", "Irf7", "Ifngr1", "Ifnar1")
sub_name <- "IFNg_responder"
tmp_spe = c("Ccl3", "Ccl4", "Xcl1", "Ifng", "Myc", "Irf8")
sub_name <- "IFNg_producer"

dir.name3=paste(dir.name2, sub_name, sep='/')
dir.create(dir.name3)
for (i in c(1:length(tmp_spe))) {
  hoge = sprintf("/%s.png", tmp_spe[i])
  file.name = paste(dir.name3, hoge, sep='')
  png(file.name, width = 512, height = 512)
  FeaturePlot(object = mBC, features.plot = tmp_spe[i], cols.use = c("grey", "red"),
              reduction.use = "FltSNE", no.legend = FALSE,  pt.size = 1.5)
  dev.off()
}

############################################################################
#Tex term, Tex prog, cytotoxic and proliferation singnatures 
load("Sub1_BDtarget_Seurat_Aoki_res0.6.rda")

meta.data <- mBC@meta.data
meta.data <- dplyr::select(meta.data, -c("LCMV.term1", "LCMV.trans1", "LCMV.Trans1"))
mBC@meta.data <- meta.data

files  <- list.files(pattern="Sig.")
dir.name <- "Seurat_plots/Gene.Score"
dir.create(dir.name)
for(i in files){
  gene_list <- read.table(i, header = TRUE)
  gene_list <- list(as.vector(gene_list$Genes))
  sig_name <- str_remove(i, "Sig.")
  sig_name <- str_remove(sig_name, ".txt")
  mBC <- AddModuleScore(object = mBC, 
                        genes.list = gene_list,
                        n.bin = 5,
                        ctrl.size = 5,
                        enrich.name=sig_name)
  
  feature <- paste(sig_name, "1", sep = "")
  
  hoge = paste(dir.name, "feature", sep = "/")
  file.name <- paste(hoge, sig_name, "png", sep = ".")
  png(file.name, width = 512, height = 512)
  FeaturePlot(object = mBC, features.plot = feature, cols.use = c("grey", "red"),
              reduction.use = "FltSNE", no.legend = FALSE,  pt.size = 1)
  dev.off()

  hoge2 = paste(dir.name, "violin", sep = "/")
  file.name2 <- paste(hoge2, sig_name, "png", sep = ".")
  png(file.name2, width = 512, height = 512)
  p <- VlnPlot(object = mBC, features.plot = feature, point.size.use = 0.1)
  plot(p)
  dev.off()
}

#スーラオブジェクトの出力
file.name="Sub1_BDtarget_Seurat_Aoki_res0.6.rda"
save(mBC, file=file.name)

#################################################################################################################
#meta.dataの出力

#raw.dataの取り出し
raw.data <- mBC@raw.data
raw.data.mat <- t(as.matrix(raw.data))

meta.data <- as.data.frame(mBC@meta.data)
write.csv(meta.data, "meta.data.csv", row.names = TRUE)

####################################特定クラスター間のDEGの抽出#######################
load("Sub1_BDtarget_Seurat_Aoki_res0.6.rda")

DEG.subset.1 = FindMarkers(mBC, ident.1 = "0", ident.2 = "3",
                           test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                           features.use = NULL, nthreads = 12,
                           adj.p.val.threshold=0.05)
DEG.subset.2 = FindMarkers(mBC, ident.1 = "3", ident.2 = "0",
                           test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                           features.use = NULL, nthreads = 12,
                           adj.p.val.threshold=0.05)

#plot statistics　（scRNAseqの成績の打ち出し)
dir.name <- "Seurat_plots/Performanca.Cluster/"
dir.create(dir.name)
sample.name=("Sub1_BDtarget_Seurat_Aoki_res0.6")
file.name=paste(dir.name, sample.name, "nGene.png", sep='')
png(file.name, width = 512, height = 400)
RidgePlot(object = mBC, features.plot = "nGene", group.by="res.0.6", nCol = 1)
dev.off()

file.name=paste(dir.name, sample.name, "nReads_log.png", sep='')
png(file.name, width = 512, height = 400)
RidgePlot(object = mBC, features.plot = "nReads.log", group.by="res.0.6", nCol = 1)
dev.off()
