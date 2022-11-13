###############################################################################
#             基于JGB的crein v6.1 基因组注释构建 org.db，构建                 #
###############################################################################

library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)
library(readr)

options(stringsAsFactors = F)

gene_anno <- read.table(file = "./annotation/annotation_info.txt",header = TRUE,
                        sep = "\t",quote = ""
                        )
gene_anno[gene_anno==""] <- NA

gene_descript <- read_delim("./annotation/master_annotation_table.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)
gene_descript[gene_descript==""] <- NA

gene_info <- gene_descript %>% dplyr::select(GID = locusName_4532,SYMBOL = geneSymbol,GENENAME = Description) #%>% na.omit()


# GO
gos <- gene_anno %>% dplyr::select(GID=locusName,GOs =GO)

gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())

gos_list <- function(x){
  the_gos <- str_split(x[2]," ",simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1],length(the_gos)),
                        GO = the_gos,
                        EVIDENCE = rep("IEA",length(the_gos))
                        
                        )
  return(df_temp)
}

gene2gol <- apply(as.matrix(gos),1,gos_list)
gene2gol_df <- do.call(rbind.data.frame,gene2gol)

gene2go <- gene2gol_df

gene2go <- na.omit(gene2go)


# KEGG
kos_list <- function(x){
  the_kos <- str_split(x[2]," ",simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1],length(the_kos)),
                        Ko = the_kos
                        )
                        
  return(df_temp)
}

gene2ko <- gene_anno %>% dplyr::select(GID = locusName,Ko = KO)
gene2ko <- na.omit(gene2ko)
gene2kol <- apply(as.matrix(gene2ko),1,kos_list)
gene2kol_df <- do.call(rbind.data.frame,gene2kol)
gene2ko <- gene2kol_df
#gene2ko$Ko <- gsub("Ko:","",gene2ko$Ko)

update_kegg <- function(json) {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  
  save(pathway2name, ko2pathway, file = "./annotation/kegg_info.RData")
}

update_kegg(json = "./annotation/ko00001_Modified.json")
load(file = "./annotation/kegg_info.RData")

gene2pathway <- left_join(gene2ko,ko2pathway, by = "Ko") %>% 
  dplyr::select(GID, Pathway)  %>% na.omit()

tax_id = "3055"
genus = "Chlamydomonas"
species = "reinhardtiisix"


gene2go <- unique(gene2go)
gene2go <- gene2go[!duplicated(gene2go),]
gene2ko <- gene2ko[!duplicated(gene2ko),]
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
gene_info <- gene_info[!duplicated(gene_info),]

gene_info[is.na(gene_info)] <-  "Not applicable"

makeOrgPackage(
  gene_info = gene_info,
  go = gene2go,
  ko = gene2ko,
  pathway = gene2pathway,
  version = "1.36.0",
  maintainer = "yxc <156732531@qq.com>",
  author = "yxc",
  outputDir = "./annotation",
  tax_id = tax_id,
  genus = genus,
  species = species,
  goTable = "go"
  
  
)


## 加载orgdb ，并进行富集分析
install.packages("./annotation/org.Creinhardtiisix.eg.db/",repos = NULL,type="source")
library(org.Creinhardtiisix.eg.db)


columns(org.Creinhardtiisix.eg.db)
keys(org.Creinhardtiisix.eg.db)

select(org.Creinhardtiisix.eg.db,keys = "Cre01.g041800_4532",columns = c("GO"))

## 常规格式化
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]","",pathway2name$Name)
pathway2name <- na.omit(pathway2name)
pathway2gene <- gene2pathway[,c("Pathway","GID")]

write.table(pathway2name,file = "./annotation/pathway2name.tsv",sep = "\t",quote = F,row.names = F)
write.table(pathway2gene,file = "./annotation/pathway2gene.tsv",sep = "\t",quote = F,row.names = F)








################################################################################
#                      日常使用代码                                            #
################################################################################
##富集分析 不需要org.DB 
library(clusterProfiler)
library(org.Creinhardtiisix.eg.db)
pathway2gene <- read.table("./annotation/pathway2gene.tsv",header = T,sep = "\t")
pathway2name <- read.table("./annotation/pathway2name.tsv",header = T,sep = "\t")


gene <- read.table("./tmp/genelist.txt",header = T)
gene_list <- gsub(".v6.1","",gene[,3]) 

#gene <- read.table("./tmp/tmp_gene.txt",header = T)
#gene_list <- gsub(".[1-9].v6.1","",gene[,1]) 

gene_list <- unique(gene_list)

# kegg
ekp <- enricher(gene_list,
                TERM2NAME = pathway2name,
                TERM2GENE = pathway2gene,
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1
                )
dotplot(ekp,showCategory = 50,color = "pvalue")
ekp_result <- ekp@result
ekp_result <- tidyr::separate(data = ekp_result,col = "GeneRatio",c("GeneRatio_gene","GeneRatio_bg"),sep = "/")
ekp_result <- tidyr::separate(data = ekp_result,col = "BgRatio",c("BgRatio_gene","BgRatio_bg"),sep = "/")
ekp_result[,3:6] <- lapply(ekp_result[,3:6],as.numeric)
ekp_result$richFactor<-ekp_result$GeneRatio_gene/ekp_result$BgRatio_gene
ekp_result[c(1:50),]["Description"]
write.table(ekp_result,file = "./result/01_AS_gene_KEGG_enrich.tsv",sep = "\t",
            quote = F, col.names = T,row.names = F)

# go
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Creinhardtiisix.eg.db,
                keyType = "GID",
                ont = "ALL", # CC BP MF ALL
                qvalueCutoff = 1,
                pvalueCutoff = 1
                
                )

dotplot(ego)
ego_result <- ego@result
ego_result <- tidyr::separate(data = ego_result,col = "GeneRatio",c("GeneRatio_gene","GeneRatio_bg"),sep = "/")
ego_result <- tidyr::separate(data = ego_result,col = "BgRatio",c("BgRatio_gene","BgRatio_bg"),sep = "/")
ego_result[,4:7] <- lapply(ego_result[,4:7],as.numeric)
ego_result$richFactor<-ego_result$GeneRatio_gene/ego_result$BgRatio_gene
write.table(ego_result,file = "./result/01_AS_gene_GO_enrich.tsv",sep = "\t",
            quote = F, col.names = T,row.names = F)

################################################################################
####################################KEGG 绘图###################################
################################################################################
###差异转录本功能富集可视化

###KEGG 富集结果
library(dplyr)
library(reshape2)
library(ggplot2)
library(readxl)

DTUG_KEGG <- read.csv(file = "./result/01_AS_gene_KEGG_enrich.tsv",header = T,sep = "\t")
DTUG_KEGG <- DTUG_KEGG[c(2,7,11,12)][1:50,]
DTUG_KEGG_order <- DTUG_KEGG[order(DTUG_KEGG$Count,decreasing = F),]
DTUG_KEGG_order$`-log10(pvalue)` <- -log10(DTUG_KEGG_order$pvalue)
DTUG_KEGG_order$Description <- factor(DTUG_KEGG_order$Description,
                                      levels = DTUG_KEGG_order$Description
)
colnames(DTUG_KEGG_order)[3] <- "Gene Count"

library(ggheatmap)
library(ggsci)
library(ggplot2)
library(ggrepel)
library(readr)
color=colorRampPalette(colors = c("#008EA0FF", "white", "#C71000FF"))(50)

p <- ggplot(DTUG_KEGG_order,aes(x=`RichFactor`,y=Description)) +
  #geom_bar(position = 'dodge',stat = 'identity',width = 0.5)+
  geom_point(aes(size=`Gene Count`,color=`-log10(pvalue)`)) +
  coord_cartesian(clip="off")+
  scale_color_gradient(low = "#008EA0E5",high = "red")+
  theme_bw(base_size = 22)+
  scale_size_continuous(
                        )+
  labs(x="RichFactor",y="KEGG TERM") + 
  scale_size_area(max_size = 6 , breaks=c(2,10,18,26,35))+
  theme(axis.text.x = element_text(angle = 0,size = 15,face = "bold",hjust = 1),
        axis.text.y = element_text(angle = 0,face = "bold",size = 10),
        axis.title.x = element_text(angle = 0,face = "bold",size = 20),
        axis.title.y = element_text(angle = 90,face = "bold",size = 20),
        legend.title = element_text(face = "bold",size = 15),
        legend.text = element_text(face = "bold",size = 13))
p



#color <- c("#FF62BC","#E76BF3", "#9590FF", "#00B0F6",
#           "#00BFC4", "#00BF7D", "#39B600", "#A3A500", "#D89000", "#F8766D")





###GO 富集结果

library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(aplot)

DTUG_GO <- read_excel(path ="./result/01_AS_gene_GO_enric.xls",sheet = "Sheet1")
DTUG_GO <- DTUG_GO[c(1,2,3,8,12,13)]
DTUG_GO <- tidyr::unite(DTUG_GO,"GO_term","ID","Description",sep = " ")
#DTUG_GO <- melt(DTUG_GO,variable.name = "value_type",value.name = "-log10(P_value)")
#DTUG_GO$`-log10(P_value)` <- -log10(DTUG_GO$`-log10(P_value)`)
#DTUG_GO <- melt(DTUG_GO,variable.name = "value_type",value.name = "基因数目")
#DTUG_GO$`-log10(P_value)` <- -log10(DTUG_GO$`-log10(P_value)`)
DTUG_GO_order <- DTUG_GO
#DTUG_GO_order <- DTUG_GO[order(DTUG_GO$基因数目,decreasing = F),]
#DTUG_GO_order <- DTUG_GO[order(DTUG_GO,decreasing = F),]
DTUG_GO_order$ONTOLOGY <- factor(DTUG_GO_order$ONTOLOGY,levels = c("Biological Process","Cellular Component","Molecular Function"))
#DTUG_GO_order <- DTUG_GO[order(DTUG_GO$基因本体类别,decreasing = F),]
DTUG_GO_order$GO_term<- factor(DTUG_GO_order$GO_term,levels = DTUG_GO$GO_term[30:1])

DTUG_GO_order$`-log10(pvalue)` <- -log10(DTUG_GO_order$pvalue)


GO_plot <- 
  #ggplot(DTUG_GO_order, aes(x = `RichFactor`,y = `GO_term`,fill = `-log10(pvalue)`))+
  
  ggplot(DTUG_GO_order, aes(x = `RichFactor`,y = `GO_term`))+
  #geom_bar(position = 'dodge',stat = 'identity',width = 0.9) + 
  geom_point(aes(size=`Gene Count`,color=`-log10(pvalue)`))+
  coord_cartesian(clip="off")+
  #scale_fill_gradient(low = "#008EA0E5",high = "red")+
  scale_color_gradient(low = "#008EA0E5",high = "red")+
  #scale_fill_manual(values = c("#C71000E5","#008EA0E5","#FF6F00E5") )+
  #facet_wrap(~Comparison,nrow = 8,strip.position = "top")+
  scale_size_area(max_size = 10 , breaks=c(3,14,25,47,69))+
  
  labs(y ="GO Term",x = "RichFactor")+
  theme_bw(base_size = 22) +
  theme(
  axis.text.x = element_text(angle = 0,size = 15,face = "bold",hjust = 1),
  axis.text.y = element_text(angle = 0,face = "bold",size = 15),
  axis.title.x = element_text(angle = 0,face = "bold",size = 20),
  axis.title.y = element_text(angle = 90,face = "bold",size = 20),
  legend.title = element_text(face = "bold",size = 15),
  legend.text = element_text(face = "bold",size = 12))

DTUG_GO_order$x <- rep(1,30)
GO_plot_1 <- ggplot(DTUG_GO_order,aes(x=x,y=`GO_term`))+
  geom_tile(aes(fill=`ONTOLOGY`))+
  labs(fill = "GO ONTOLOGY")+
  theme_bw(base_size = 22)+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
  )+
  scale_fill_manual(values = pal_futurama("planetexpress", alpha = 1)(10)[c(1,2,3)]) +
  theme(
    legend.title = element_text(face = "bold",size = 15),
    legend.text = element_text(face = "bold",size = 12))
#GO_plot_1
GO_plot_2 <- GO_plot%>%
  insert_right(GO_plot_1,width = 0.07) 


GO_plot_2
