library(dplyr)
library(reshape2)
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(aplot)
library(stringr)
#write.table(AS_cluster, './result/02_AS_cluster.txt', sep = '\t', col.names = NA, quote = FALSE)
gene_anno_1 <- read_delim(file = "./annotation/gene_anno.tsv",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
gene_anno_2 <- read_delim(file = "./annotation/annotation_info.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
gene_anno_1 <- gene_anno_1[,c(1,5,9)]
gene_anno_2 <- unique(gene_anno_2[,c(2,5)])
colnames(gene_anno_1)[1] <- colnames(gene_anno_2)[1]

Pfam_ID_list <- function(x){
  the_gos <- str_split(x[2]," ",simplify = FALSE)[[1]]
  df_temp <- data.frame(locusName = rep(x[1],length(the_gos)),
                        Pfam_ID = the_gos
                        )
                        
  return(df_temp)
}


merge_pfam <- function(x){
  b = unique(x)
  pfam = str_c(b,collapse=',')
  return(pfam)
}


gene2pfam <- apply(as.matrix(gene_anno_2),1,Pfam_ID_list)
gene2pfam_df <- do.call(rbind.data.frame,gene2pfam)
gene_anno_2 <- gene2pfam_df
gene_anno_2 <- na.omit(gene_anno_2)
gene_anno_2 <- aggregate(Pfam_ID~.,gene_anno_2,FUN = merge_pfam)
gene_anno_1$locusName <- gsub(".v6.1","",gene_anno_1$locusName)
gene_anno <- left_join(gene_anno_1,gene_anno_2,by="locusName")


AS_cluster <- read.csv(file = "./result/02_AS_cluster.txt",sep = "\t",header = T,row.names = 1)
AS_cluster_df <- data.frame(AS_cluster)
AS_cluster_1 <- AS_cluster_df[
  AS_cluster_df["AS_cluster"]==2|AS_cluster_df["AS_cluster"]==1,#|AS_cluster_df["AS_cluster"]==5|AS_cluster_df["AS_cluster"]==6,
]
AS <- data.frame(AS_event = rownames(AS_cluster_1))
AS <- separate(data = AS,col = "AS_event",c("AS_Gene_ID","AS_index"),sep = ".v6.1_")

AS <- separate(data = AS,col = "AS_Gene_ID",c("AS","Gene_ID"),sep = "_",extra = "merge")

gene_AS <- AS[,c(2,1)]
gene_AS <- unique(gene_AS)
gene_AS <- aggregate(AS~.,gene_AS,FUN = merge_pfam)

################################################################################
transcript_factor <- read.csv(file = "./annotation/Cre_TF_list.txt",sep = "\t")
transcript_factor <- unique(transcript_factor[,c(2,3)])
transcript_factor$Gene_ID <- paste(transcript_factor$Gene_ID,rep("_4532",length(transcript_factor$Gene_ID)),sep = "")

colnames(gene_anno)[1] <- colnames(gene_AS)[1]
colnames(transcript_factor)[1] <- colnames(gene_AS)[1]
result_anno <- left_join(transcript_factor,gene_AS,by="Gene_ID") %>% left_join(gene_anno,by="Gene_ID")
result_anno <- result_anno[,c(1,4,3,6,5,2)]
result_anno <- result_anno[!is.na(result_anno$AS),]
write.table(result_anno,file = "./result/04_AS_gene_anno_transcript_factor.tsv",sep = "\t",row.names = F,col.names = T,quote = F)

################################################################################
KEGG_result <- read.csv(file = "./result/01_AS_gene_KEGG_enrich.tsv",header = T,sep = "\t")
KEGG <- KEGG_result[1:50,c(2,10)]
rownames(KEGG) <- KEGG$Description
GO <- read_excel(path ="./result/01_AS_gene_GO_enric.xls",sheet = "Sheet1")



################################################################################
KEGG_pathway_metabolism_list <- c("Glyoxylate and dicarboxylate metabolism","Citrate cycle (TCA cycle)",
                       "Peroxisome","Glycolysis / Gluconeogenesis","Pyruvate metabolism","Fatty acid degradation",
                       "Valine, leucine and isoleucine degradation","Fatty acid biosynthesis",
                       "Fatty acid elongation","Arginine biosynthesis")



KEGG_pathway_metabolism_pd <- KEGG#[KEGG_pathway_metabolism_list,]

ID_list <- function(x){
  the_gos <- str_split(x[2],"/",simplify = FALSE)[[1]]
  df_temp <- data.frame(KEGG = rep(x[1],length(the_gos)),
                        gene_ID = the_gos
  )
  
  return(df_temp)
}

pathway2gene <- apply(as.matrix(KEGG_pathway_metabolism_pd),1,ID_list)
gene2pathway_pd <- do.call(rbind.data.frame,pathway2gene)[,c(2,1)]
gene2pathway_pd <- aggregate(KEGG~.,gene2pathway_pd,FUN = merge_pfam)

colnames(gene_anno)[1] <- colnames(gene_AS)[1]
colnames(gene2pathway_pd)[1] <- colnames(gene_AS)[1]
result_anno <- left_join(gene2pathway_pd,gene_AS,by="Gene_ID") %>% left_join(gene_anno,by="Gene_ID")
result_anno <- result_anno[,c(1,4,3,6,5,2)]


write.table(result_anno,file = "./result/04_AS_gene_anno_metabolism.tsv",sep = "\t",row.names = F,col.names = T,quote = F)


################################################################################

KEGG_pathway_metabolism_pd <- dplyr::filter(.data = KEGG,grepl("signaling",Description))

ID_list <- function(x){
  the_gos <- str_split(x[2],"/",simplify = FALSE)[[1]]
  df_temp <- data.frame(KEGG = rep(x[1],length(the_gos)),
                        gene_ID = the_gos
  )
  
  return(df_temp)
}

pathway2gene <- apply(as.matrix(KEGG_pathway_metabolism_pd),1,ID_list)
gene2pathway_pd <- do.call(rbind.data.frame,pathway2gene)[,c(2,1)]
gene2pathway_pd <- aggregate(KEGG~.,gene2pathway_pd,FUN = merge_pfam)

colnames(gene_anno)[1] <- colnames(gene_AS)[1]
colnames(gene2pathway_pd)[1] <- colnames(gene_AS)[1]
result_anno <- left_join(gene2pathway_pd,gene_AS,by="Gene_ID") %>% left_join(gene_anno,by="Gene_ID")
result_anno <- result_anno[,c(1,4,3,6,5,2)]


write.table(result_anno,file = "./result/04_AS_gene_anno_signal.tsv",sep = "\t",row.names = F,col.names = T,quote = F)


################################################################################

KEGG_pathway_metabolism_pd <- dplyr::filter(.data = KEGG,grepl("ransporters",Description))

ID_list <- function(x){
  the_gos <- str_split(x[2],"/",simplify = FALSE)[[1]]
  df_temp <- data.frame(KEGG = rep(x[1],length(the_gos)),
                        gene_ID = the_gos
  )
  
  return(df_temp)
}

pathway2gene <- apply(as.matrix(KEGG_pathway_metabolism_pd),1,ID_list)
gene2pathway_pd <- do.call(rbind.data.frame,pathway2gene)[,c(2,1)]
gene2pathway_pd <- aggregate(KEGG~.,gene2pathway_pd,FUN = merge_pfam)

colnames(gene_anno)[1] <- colnames(gene_AS)[1]
colnames(gene2pathway_pd)[1] <- colnames(gene_AS)[1]
result_anno <- left_join(gene2pathway_pd,gene_AS,by="Gene_ID") %>% left_join(gene_anno,by="Gene_ID")
result_anno <- result_anno[,c(1,4,3,6,5,2)]


write.table(result_anno,file = "./result/04_AS_gene_anno_transporters.tsv",sep = "\t",row.names = F,col.names = T,quote = F)


################################################################################
KEGG_pathway_metabolism_pd <- dplyr::filter(.data = KEGG,grepl("Spliceosome",Description))

ID_list <- function(x){
  the_gos <- str_split(x[2],"/",simplify = FALSE)[[1]]
  df_temp <- data.frame(KEGG = rep(x[1],length(the_gos)),
                        gene_ID = the_gos
  )
  
  return(df_temp)
}

pathway2gene <- apply(as.matrix(KEGG_pathway_metabolism_pd),1,ID_list)
gene2pathway_pd <- do.call(rbind.data.frame,pathway2gene)[,c(2,1)]
gene2pathway_pd <- aggregate(KEGG~.,gene2pathway_pd,FUN = merge_pfam)

colnames(gene_anno)[1] <- colnames(gene_AS)[1]
colnames(gene2pathway_pd)[1] <- colnames(gene_AS)[1]
result_anno <- left_join(gene2pathway_pd,gene_AS,by="Gene_ID") %>% left_join(gene_anno,by="Gene_ID")
result_anno <- result_anno[,c(1,4,3,6,5,2)]


write.table(result_anno,file = "./result/04_AS_gene_anno_spiceosome.tsv",sep = "\t",row.names = F,col.names = T,quote = F)