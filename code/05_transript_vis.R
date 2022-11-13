##
library(magrittr)
library(dplyr)
library(ggtranscript)
library(ggplot2)
library(rtracklayer)
library(ggheatmap)
library(patchwork)
library(paletteer)
library(ggsci)

gtf_path <- "./annotation/CreinhardtiiCC_4532_707_v6.1.gene_exons.gtf"
gtf <- rtracklayer::import(gtf_path)
gtf <- gtf %>% dplyr::as_tibble()

ASlevel <- read.table(file="./result/02_AS_cluster.txt",sep = "\t",header = T)
colnames(ASlevel)[1] <- "AS"
ASlevel <- ASlevel[,c(1:10)]

ASlevel_pd <- dplyr::filter(.data = ASlevel,grepl("SRS1",AS)|grepl("SRE1",AS)| grepl("PGI1",AS)|grepl("PKS3",AS))
ASlevel_pd <- ASlevel_pd[c(2,9,8,7,1,3,4,5,6,10),]
ASlevel_pd$AS <- gsub(".v6.1","",ASlevel_pd$AS)

#color.1 <- colorRampPalette(rev(c("#C71000FF","white","#008EA0FF")))(100)

##mat <- matrix( c(1,1,1,1,2,3,4,5),ncol=4,nrow=2,byrow=TRUE)
#par(oma = c(2,2,2,2))
#l   <- layout(mat,widths = c(1,1,1,1), heights = c(1, 1))


#sample_cor <- read_csv(file = "./result/03_sampel_cor.csv")
sample_AS <- ASlevel_pd
#sample_cor <- tidyr::separate(sample_cor,col = "AS",into = c("AS_gene","name_index"),sep = ".v6.1_")
#sample_cor <- tidyr::separate(sample_cor,col = "name_index",into = c("name","index"),sep = "_",extra ="merge")
#sample_cor <- tidyr::unite(sample_cor,col = "AS",c("AS_gene","name"),sep = "_")
#sample_cor <- sample_cor[,c(1,3:11)]
rownames(sample_AS) <- sample_AS$AS

sample_AS <- sample_AS[-1]
colnames(sample_AS) <- gsub("N_S_","" ,gsub("_ASlevel","", colnames(sample_AS))    )
sample_list <- colnames(sample_AS) #gsub("N_S_","" ,gsub("_ASlevel","", colnames(sample_cor))    )
AS_list <- factor(rownames(sample_AS)[10:1],levels = rownames(sample_AS)[10:1])
sample_AS_S <- t(scale(t(sample_AS),scale = T,center = T))
sample_AS_plot <- ggheatmap(sample_AS_S,cluster_rows = F,cluster_cols = F,
                             scale = "none",legendName = "AS_level",
                             text_position_rows="left",
                             levels_rows=AS_list,levels_cols=sample_list,
                             #levels_rows=sample_list,levels_cols=sample_list,
                             border=NA,text_position_cols = "bottom",
                             color=colorRampPalette(colors = c("#008EA0FF", "white", "#C71000FF"))(50)
) + 
  theme(axis.text.x = element_text(angle = 0,size = 15,face = "bold",hjust = 1),
                        axis.text.y = element_text(angle = 0,face = "bold",size = 8),
                        legend.title = element_text(face = "bold",size = 10),
                        legend.text = element_text(face = "bold",size = 10))  #+scale_fill_gsea()
sample_AS_plot

################################################################################
trancript_exp <- read.csv(file = "./Quantitative/Stringtie_Count/transcript_TPM_marix.tsv",
                          sep = "\t")
colnames(trancript_exp)[1] <- "transcript_id_ref"
trancript_exp <- tidyr::separate(trancript_exp,col = "transcript_id_ref",
                                 into = c("transcript_id","ref"),sep = ".v6.1")
trancript_exp$ref <- gsub("_","",trancript_exp$ref)

for (sam in colnames(sample_AS)) {
  trancript_exp[sam] <- apply(trancript_exp[,c( paste("N_S",sam,"rep1",sep = "_"), paste("N_S",sam,"rep2",sep = "_"))],1,mean)
}
trancript_exp <- trancript_exp[,c(colnames(trancript_exp)[1:2],colnames(sample_AS))]

#trancript_exp <- dplyr::filter(.data = trancript_exp,
#              grepl("SRS1",ref)|grepl("SRE1",ref)| grepl("PGI1",ref)|grepl("PKS3",ref))
#trancript_exp <- trancript_exp[c(1:3,8,9,4:7,13:16),]
#write.table(x=trancript_exp$transcript_id,file = "./annotation/transcript_id.txt",quote = F,
#            row.names = F,col.names = F)

trancript_exp <- tidyr::unite(trancript_exp,col = "transcript",c("transcript_id","ref"),sep = "_")

trancript_exp_pd <- dplyr::filter(.data = trancript_exp,
                    grepl("SRS1",transcript)|grepl("SRE1",transcript)| grepl("PGI1",transcript)|grepl("PKS3",transcript))
trancript_exp_pd <- trancript_exp_pd[c(1:3,8,9,4:7,13:16),]
sample_cor <- trancript_exp_pd

rownames(sample_cor) <- sample_cor$transcript

sample_cor <- sample_cor[-1]
#colnames(sample_cor) <- gsub("N_S_","" ,gsub("_ASlevel","", colnames(sample_cor))    )
sample_list <- colnames(sample_cor) #gsub("N_S_","" ,gsub("_ASlevel","", colnames(sample_cor))    )
AS_list <- factor(rownames(sample_cor)[13:1],levels = rownames(sample_cor)[13:1])
sample_cor_s <- t(scale(t(sample_cor),scale = T,center = T))
sample_exp_plot <- ggheatmap(sample_cor_s,cluster_rows = F,cluster_cols = F,
                             scale = "none",legendName = "Expression_level",
                             text_position_rows="left",
                             levels_rows=AS_list,levels_cols=sample_list,
                             #levels_rows=sample_list,levels_cols=sample_list,
                             border=NA,text_position_cols = "bottom",
                             color=colorRampPalette(colors = c("#008EA0FF", "white", "#C71000FF"))(50)
) + theme(axis.text.x = element_text(angle = 0,size = 15,face = "bold",hjust = 1),
          axis.text.y = element_text(angle = 0,face = "bold",size = 10),
          legend.title = element_text(face = "bold",size = 10),
          legend.text = element_text(face = "bold",size = 10))  #+scale_fill_gsea()
sample_exp_plot
################################################################################
sample_AS_plot / sample_exp_plot 
################################################################################
# 选取特定基因绘制转录本
gene_of_interest <-"SRS1" # "PGI1" #,"PKS3"
gene_annotation_from_gtf <- gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 
gene_annotation_from_gtf <- gene_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_id
  )

gene_annotation_from_gtf$transcript_id <- gsub(".v6.1","",gene_annotation_from_gtf$transcript_id)
gene_annotation_from_gtf$transcript_id <- factor(gene_annotation_from_gtf$transcript_id,levels = sort(unique(gene_annotation_from_gtf$transcript_id))[3:1])

gene_annotation_from_gtf$CDS_CHANGE_TYPE <- rep("NOT CHANGE",length(gene_annotation_from_gtf$transcript_id))

gene_annotation_from_gtf$CDS_CHANGE_TYPE <- factor(gene_annotation_from_gtf$CDS_CHANGE_TYPE,levels = c("NOT CHANGE","CHANGE"))
gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
gene_cds <- gene_annotation_from_gtf %>% dplyr::filter(type == "CDS")

SRS1 <- gene_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = CDS_CHANGE_TYPE),
    height = 0.25,
    #fill = "white", 
  ) +
  
  geom_range(
    data = gene_cds
    ,aes(fill = CDS_CHANGE_TYPE)
    #,fill = "black"
    
  ) + 
  scale_fill_manual(values = pal_futurama("planetexpress", alpha = 1)(10)[c(3,2)])+
  geom_intron(
    data = to_intron(gene_exons, "transcript_id"),
    aes(strand = strand)
  )+ labs(y ="SRS1",x = "")+theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 10,face = "bold",hjust = 1),
        axis.text.y = element_text(angle = 0,face = "bold",size = 12),
        axis.title.y = element_text(angle = 90,face = "bold",size = 15),
        legend.position = 'none')

  
SRS1

################################################################################
# 选取特定基因绘制转录本
gene_of_interest <-"SRE1"#"SRS1" # "PGI1" #,"PKS3"
gene_annotation_from_gtf <- gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 
gene_annotation_from_gtf <- gene_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_id
  )

gene_annotation_from_gtf %>% head()
gene_annotation_from_gtf$transcript_id <- gsub(".v6.1","",gene_annotation_from_gtf$transcript_id)
gene_annotation_from_gtf$transcript_id <- factor(gene_annotation_from_gtf$transcript_id,levels = sort(unique(gene_annotation_from_gtf$transcript_id))[2:1])
gene_annotation_from_gtf$CDS_CHANGE_TYPE <- c(rep("CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre12.g498350_4532.2", na.rm = T)),
                                  rep("NOT CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre12.g498350_4532.1", na.rm = T))
                                  )
gene_annotation_from_gtf$CDS_CHANGE_TYPE <- factor(gene_annotation_from_gtf$CDS_CHANGE_TYPE,levels = c("NOT CHANGE","CHANGE"))
gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
gene_cds <- gene_annotation_from_gtf %>% dplyr::filter(type == "CDS")

SRE1 <- gene_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = CDS_CHANGE_TYPE),
    height = 0.25,
    #fill = "white", 
  ) +
  
  geom_range(
    data = gene_cds
    ,aes(fill = CDS_CHANGE_TYPE)
    #,fill = "black"
    
  ) + 
  scale_fill_manual(values = pal_futurama("planetexpress", alpha = 1)(10)[c(3,2)])+
  geom_intron(
    data = to_intron(gene_exons, "transcript_id"),
    aes(strand = strand)
  )+ labs(y ="SRE1",x = "")+theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 10,face = "bold",hjust = 1),
        axis.text.y = element_text(angle = 0,face = "bold",size = 12),
        axis.title.y = element_text(angle = 90,face = "bold",size = 15),
    legend.position = 'right')
SRE1 

################################################################################
# 选取特定基因绘制转录本
gene_of_interest <-"PGI1" #,"PKS3"
gene_annotation_from_gtf <- gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 
gene_annotation_from_gtf <- gene_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_id
  )

gene_annotation_from_gtf %>% head()
gene_annotation_from_gtf$transcript_id <- gsub(".v6.1","",gene_annotation_from_gtf$transcript_id)
gene_annotation_from_gtf$transcript_id <- factor(gene_annotation_from_gtf$transcript_id,levels = sort(unique(gene_annotation_from_gtf$transcript_id))[4:1])
gene_annotation_from_gtf$CDS_CHANGE_TYPE <- c(rep("CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre03.g175400_4532.3", na.rm = T)),
                                  rep("CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre03.g175400_4532.4", na.rm = T)),
                                  rep("NOT CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre03.g175400_4532.2", na.rm = T)),
                                  rep("NOT CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre03.g175400_4532.1", na.rm = T))
)
gene_annotation_from_gtf$CDS_CHANGE_TYPE <- factor(gene_annotation_from_gtf$CDS_CHANGE_TYPE,levels = c("NOT CHANGE","CHANGE"))
gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
gene_cds <- gene_annotation_from_gtf %>% dplyr::filter(type == "CDS")

PGI1 <- gene_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = CDS_CHANGE_TYPE),
    height = 0.25,
    #fill = "white", 
  ) +
  
  geom_range(
    data = gene_cds
    ,aes(fill = CDS_CHANGE_TYPE)
    #,fill = "black"
    
  ) + 
  scale_fill_manual(values = pal_futurama("planetexpress", alpha = 1)(10)[c(3,2)])+
  geom_intron(
    data = to_intron(gene_exons, "transcript_id"),
    aes(strand = strand)
  )+ labs(y ="PGI1",x = "")+theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 10,face = "bold",hjust = 1),
        axis.text.y = element_text(angle = 0,face = "bold",size = 12),
        axis.title.y = element_text(angle = 90,face = "bold",size = 15),
        legend.position = 'none')
PGI1 

################################################################################
# 选取特定基因绘制转录本
gene_of_interest <-"PKS3" #"SRE1"#"SRS1" # "PGI1" #,"PKS3"
gene_annotation_from_gtf <- gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 
gene_annotation_from_gtf <- gene_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_id
  )

gene_annotation_from_gtf %>% head()
gene_annotation_from_gtf$transcript_id <- gsub(".v6.1","",gene_annotation_from_gtf$transcript_id)
gene_annotation_from_gtf$transcript_id <- factor(gene_annotation_from_gtf$transcript_id,levels = sort(unique(gene_annotation_from_gtf$transcript_id))[4:1])
gene_annotation_from_gtf$CDS_CHANGE_TYPE <- c(rep("CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre17.g722150_4532.2", na.rm = T)),
                                  rep("CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre17.g722150_4532.3", na.rm = T)),
                                  rep("NOT CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre17.g722150_4532.1", na.rm = T)),
                                  rep("CHANGE",sum(gene_annotation_from_gtf$transcript_id == "Cre17.g722150_4532.4", na.rm = T))
)
gene_annotation_from_gtf$CDS_CHANGE_TYPE <- factor(gene_annotation_from_gtf$CDS_CHANGE_TYPE,levels = c("NOT CHANGE","CHANGE"))
gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
gene_cds <- gene_annotation_from_gtf %>% dplyr::filter(type == "CDS")

PKS3 <- gene_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = CDS_CHANGE_TYPE),
    height = 0.25,
    #fill = "white", 
  ) +
  
  geom_range(
    data = gene_cds
    ,aes(fill = CDS_CHANGE_TYPE)
    #,fill = "black"
    
  ) + 
  scale_fill_manual(values = pal_futurama("planetexpress", alpha = 1)(10)[c(3,2)])+
  geom_intron(
    data = to_intron(gene_exons, "transcript_id"),
    aes(strand = strand)
  )+ labs(y ="PKS3",x = "")+theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 10,face = "bold",hjust = 1),
        axis.text.y = element_text(angle = 0,face = "bold",size = 12),
        axis.title.y = element_text(angle = 90,face = "bold",size = 15),
        legend.position = 'none')
PKS3 
################################################################################
sample_AS_plot / (SRS1 | SRE1 )/( PGI1 |  PKS3)
layout <- "
#111
3456
"
sample_cor_plot + SRS1 + SRE1 + PGI1 +  PKS3 + plot_layout(design = layout)
################################################################################
# extract exons and cds for the MANE-select transcript

gene_cds <- gene_annotation_from_gtf %>% dplyr::filter(type == "CDS")

gene_transcipt_exons <- gene_exons %>% dplyr::filter(transcript_id == "Cre17.g722150_4532.2.v6.1")
gene_transcipt_cds <- gene_cds %>% dplyr::filter(transcript_id == "Cre17.g722150_4532.2.v6.1")




# add transcript name column to junctions for plotting
# sod1_junctions <- sod1_junctions %>% dplyr::mutate(transcript_name = "SOD1-201")

gene_transcipt_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    fill = "white", 
    height = 0.25
  ) +
  geom_range(
    data = gene_transcipt_cds
  ) + 
  geom_intron(
    data = to_intron(gene_transcipt_exons, "transcript_id")
  ) + 
  scale_size_continuous(range = c(0.1, 1))
  #+ 
  geom_junction(
    data = sod1_junctions,
    aes(size = mean_count),
    junction.y.max = 0.5
  ) +
  geom_junction_label_repel(
    data = sod1_junctions,
    aes(label = round(mean_count, 2)),
    junction.y.max = 0.5
  ) + 
  scale_size_continuous(range = c(0.1, 1))

# extract exons and cds for the MANE-select transcript
sod1_201_exons <- sod1_exons %>% dplyr::filter(transcript_name == "SOD1-201")
sod1_201_cds <- sod1_cds %>% dplyr::filter(transcript_name == "SOD1-201")



