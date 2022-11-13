library(ggheatmap)
library(ggsci)
library(ggplot2)
library(ggrepel)
library(readr)

# 导入基因表达数据TPM数据
sample_list <- c('N_S_0h_rep1','N_S_0h_rep2','N_S_10m_rep1','N_S_10m_rep2',
                 'N_S_30m_rep1','N_S_30m_rep2','N_S_1h_rep1','N_S_1h_rep2',
                 'N_S_2h_rep1','N_S_2h_rep2','N_S_6h_rep1','N_S_6h_rep2',
                 'N_S_8h_rep1','N_S_8h_rep2','N_S_24h_rep1','N_S_24h_rep2',
                 'N_S_48h_rep1','N_S_48h_rep2')

for (Sample in sample_list) {
  gene_abund <- read_delim(paste("./Quantitative/Stringtie/",
                                 Sample,"_gene_abund.tab",sep = ''),
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)
  gene_abund <- gene_abund[,c(1,9)]
  colnames(gene_abund) <- c("Gene_ID",Sample)
  assign(paste(Sample,"_gene_abund",sep = ''),gene_abund)
}
# 生成基因表达矩阵
library(tidyverse)
gene_abund <- N_S_0h_rep1_gene_abund
sample_list_2 <- c('N_S_0h_rep2','N_S_10m_rep1','N_S_10m_rep2',
                   'N_S_30m_rep1','N_S_30m_rep2','N_S_1h_rep1','N_S_1h_rep2',
                   'N_S_2h_rep1','N_S_2h_rep2','N_S_6h_rep1','N_S_6h_rep2',
                   'N_S_8h_rep1','N_S_8h_rep2','N_S_24h_rep1','N_S_24h_rep2',
                   'N_S_48h_rep1','N_S_48h_rep2')
for (Sample_2 in sample_list_2) {
  gene_abund <- left_join(gene_abund,get(paste(Sample_2,"_gene_abund",sep = '')),
                          by="Gene_ID")
}

# 进行样品相关性分析
gene_abund_matrix <- gene_abund[,2:19]
colnames(gene_abund_matrix) <- c('0h_1','0h_2','10m_1','10m_2',
                                 '30m_1','30m_2','1h_1','1h_2',
                                 '2h_1','2h_2','6h_1','6h_2',
                                 '8h_1','8h_2','24h_1','24h_2',
                                 '48h_1','48h_2')
rownames(gene_abund_matrix) <- gene_abund$Gene_ID
gene_abund_matrix <- as.matrix(gene_abund_matrix)
library(corrplot)
sampel_cor <- cor(gene_abund_matrix)
corrplot(sampel_cor,method = "pie",tl.col = "black",tl.cex = 0.7)

# 样品相关性热图绘制

library (ggplot2)
library (reshape2)#数据转换
require(scales)#数据缩放
library(ggheatmap)
library(pheatmap)
library(aplot)
# sample_col <- data.frame(class = rep(c("4B33","CC-5325","Mut-4"),c(3,3,3)))
# row.names(sample_col) <- colnames(mat)
# ann_colors = list(class = c(`4B33` = "#289B3F", `CC-5325` = "#B3B321",`Mut-4`="#E06726"))





## PCA 
library(ggplot2)
library(ggrepel)
library(ggsci)


colnames(gene_abund_matrix) <- c('0h_1','0h_2','10m_1','10m_2',
                                 '30m_1','30m_2','1h_1','1h_2',
                                 '2h_1','2h_2','6h_1','6h_2',
                                 '8h_1','8h_2','24h_1','24h_2',
                                 '48h_1','48h_2')
gene_matrix <- t(gene_abund_matrix)+0.0001
pca <- prcomp(gene_matrix)
pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
pca_data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca_data$`时间点` <- factor(c("0h","0h","10m","10m","30m","30m",
                           "1h","1h","2h","2h",
                           "6h","6h","8h","8h",
                           "24h","24h","48h","48h"),
                         levels = c("0h","10m","30m","1h","2h","6h","8h",
                                    "24h","48h"))



write.csv(pca_data,file = "./result/03_pca_data.csv")
write.csv(sampel_cor,file = "./result/03_sampel_cor.csv")





## Cre08.g378150


## Cre12.g530650





################################################################################

#01 样本热图

sample_cor <- read_csv(file = "./result/03_sampel_cor.csv")
sample_cor <- as.data.frame(sample_cor)
rownames(sample_cor) <- sample_cor$...1

sample_cor <- sample_cor[-1]

sample_list <- c('0h_1','0h_2','10m_1','10m_2',
                 '30m_1','30m_2','1h_1','1h_2',
                 '2h_1','2h_2','6h_1','6h_2',
                 '8h_1','8h_2','24h_1','24h_2',
                 '48h_1','48h_2')
sample_cor_plot <- ggheatmap(sample_cor,cluster_rows = F,cluster_cols = F,
                             scale = "none",legendName = "PCC",
                             text_position_rows="left",
                             levels_rows=sample_list,levels_cols=sample_list,
                             #levels_rows=sample_list,levels_cols=sample_list,
                             border = "black",text_position_cols = "bottom",
                             color=colorRampPalette(colors = c("#008EA0FF", "white", "#C71000FF"))(50)
)
sample_cor_plot




sample_cor_plot <- sample_cor_plot+theme(axis.text.x = element_text(angle = 60,size = 15,face = "bold",hjust = 1),
                                         axis.text.y = element_text(angle = 0,face = "bold",size = 15),
                                         legend.title = element_text(face = "bold",size = 18),
                                         legend.text = element_text(face = "bold",size = 15))  #+scale_fill_gsea()

sample_cor_plot


#02 PCA 二维散点图

pca_data <- read.csv("./result/03_pca_data.csv", header=T)

#colnames(pca_data) <- as.vector(pca_data[1,])
rownames(pca_data) <- as.vector(pca_data[,1])
#pca_data <- pca_data[-1,]
pca_data <- pca_data[,-1]
#pca_data <- pca_data[,-2]
pca_data$`时间点` <- factor(c("0h","0h","10m","10m","30m","30m",
                           "1h","1h","2h","2h",
                           "6h","6h","8h","8h",
                           "24h","24h","48h","48h"),
                         levels = c("0h","10m","30m","1h","2h","6h","8h",
                                    "24h","48h"))
pca_data$X <- as.numeric(pca_data$X)
pca_data$Y <- as.numeric(pca_data$Y)
ggplot(data=pca_data,aes(x=X,y=Y,fill=`时间点`))+
  geom_point(size=8,shape=22 )+theme_bw(base_size = 44)+
  
  scale_fill_manual(values  = pal_futurama("planetexpress", alpha = 1)(9) )+
  
  xlab(paste("PC1","72.1%","variance",sep=" "))+
  ylab(paste("PC2","16.4%","variance",sep=" "))+
  theme(legend.position = c(0.8,0.9))+
  geom_text_repel(aes(x=X,y=Y,label = rownames(pca_data)),size=5,box.padding = 0.5,
                  point.padding = 0.5,min.segment.length = 0.5,segment.color = NA,color = "black")+
  labs(fill = " ")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = 'right',
        axis.text = element_text(size = 15,face = "bold"),axis.title = element_text(size = 18,face = "bold"),
        legend.text = element_text(size = 15,face = "bold"),legend.title = element_text(size = 20,face = "bold"),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm')
  )

