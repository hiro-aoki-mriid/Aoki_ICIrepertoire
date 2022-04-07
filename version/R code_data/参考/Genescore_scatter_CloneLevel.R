#scTCRseq, outputを用いた解析
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

#Tumor中の頻度を%, log10 scaleに変換
d_sub$freq_tumor <- log10(d_sub$freq_tumor*100)

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
    tiff(name_out, width=1.2*ppi, height=1.2*ppi, res=ppi)
    p <- ggplot(d_sub2, aes(x=x_axis, y=y_axis, colour=class)) +
      scale_colour_manual(values=c(Poly="blue", Oligo="red"))+
      geom_jitter(alpha=0.8, size=0.8) +
      theme_bw(base_size = 6) +
      theme(
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(family="Arial"),
        axis.text.y = element_text(family="Arial")) +
      guides(colour=FALSE)
    plot(p)
    dev.off()
    
  }
}
