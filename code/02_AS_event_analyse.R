##富集分析 不需要org.DB 

library(readr)
library(tidyverse)

timepoint_list <- list.dirs(path ="./rmats",full.names = F, recursive = F)[1:9][c(1,2,6,3,5,8,9,4,7)]

AStype_list <- c("A3SS","A5SS","SE","RI","MXE")


rmatsdata_process <- function(timepoint, AStype){
  
  rmatsdata <- read_delim (paste("rmats/",timepoint,"/",AStype,".MATS.JC.txt",sep = ""),  #"rmats/N_S_0h/A3SS.MATS.JC.txt",
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  MATS_JC <- rmatsdata
  
  MATS_JC["AS_type"] <- rep(AStype,nrow(MATS_JC))
  
  MATS_JC <- separate(data = MATS_JC,col = IJC_SAMPLE_1,
                      into = c( paste( timepoint,"IJC_SAMPLE_1",sep = "_"),paste( timepoint,"IJC_SAMPLE_2",sep = "_")),
                      sep = ",")
  
  MATS_JC <- separate(data = MATS_JC,col = SJC_SAMPLE_1,
                      into = c( paste( timepoint,"SJC_SAMPLE_1",sep = "_"),paste( timepoint,"SJC_SAMPLE_2",sep = "_")),
                      sep = ",")
  
  MATS_JC <- separate(data = MATS_JC,col = IncLevel1,
                      into = c( paste( timepoint,"IncLevel1",sep = "_"),paste( timepoint,"IncLevel2",sep = "_")),
                      sep = ",")
  if(AStype != "MXE"){
    MATS_JC <- MATS_JC[,c(2:11,13,14,15,16,23,24,27)]
    MATS_JC <- unite(data = MATS_JC,AS,c("AS_type",1:10),sep = "_")
  }else{
    MATS_JC <- MATS_JC[,c(2:13,15,16,17,18,25,26,29)]
    MATS_JC <- unite(data = MATS_JC,AS,c("AS_type",1:12),sep = "_")
  }
  
  
  MATS_JC[,2:7]<-lapply(MATS_JC[,2:7],as.numeric)
  
  MATS_JC[is.na(MATS_JC)] <- 0
  
  if (AStype == "RI"){
   MATS_JC[paste( timepoint,"ASlevel",sep = "_")] = apply(MATS_JC[c( paste( timepoint,"IncLevel1",sep = "_"),paste( timepoint,"IncLevel2",sep = "_"))],1,mean)
  }else{
    MATS_JC[paste( timepoint,"ASlevel",sep = "_")] = 1 - apply(MATS_JC[c( paste( timepoint,"IncLevel1",sep = "_"),paste( timepoint,"IncLevel2",sep = "_"))],1,mean)
  }
  
  
  return(MATS_JC)
}




#for (timepoint in timepoint_list) {
  
  
#  for (AStype in AStype_list) {
    
#      tmp = rmatsdata_process(timepoint,AStype)
      
#     assign(paste(AStype,timepoint,sep = "_"),tmp)
        
#  }
  
#}

for (AStype in AStype_list) {
  
  AS_temp = list()
  
  
  for (timepoint in timepoint_list) {
    tmp = rmatsdata_process(timepoint,AStype)
    tmp["min"] = apply(tmp[,4:5],1,min)
    tmp <- tmp[tmp["min"]>0,]
    tmp <-tmp[,c("AS",paste( timepoint,"ASlevel",sep = "_"))]
    
    AS_temp[[timepoint]] <- tmp
    
  }
  assign(paste(AStype,"list",sep = "_"),AS_temp)
  
  
}


AS_ALL_list = list()
for (AStype in AStype_list) {
  
  tmp_list = get(paste(AStype,"list",sep = "_"))
  
  tmp_all <-  tmp_list %>% reduce(full_join, by = "AS")
  #tmp_all[is.na(tmp_all)] <- 0.0001
  tmp_all <- na.omit(tmp_all)
  
  #tmp_all["min"] <- apply(tmp_all[2:10], 1,min)
  #tmp_all <- tmp_all[tmp_all["min"]>0,]
  
  AS_ALL_list[[AStype]] <- tmp_all
  
  #assign(paste(AStype,"all",sep = "_"),tmp_all)
}


AS_ALL <- do.call("rbind", AS_ALL_list)
write.table(AS_ALL,file = "./tmp/AS_all_nona.tsv",sep = "\t",row.names = F,col.names = T,quote = F)

# 时序分析 Mfuzz
library(Mfuzz)


AS_LEVEL <- read.table(file = "./tmp/AS_all_nona.tsv",sep = "\t",row.names = 1, header = T)
AS_LEVEL <- AS_LEVEL[,1:9]
AS_LEVEL <- as.matrix(AS_LEVEL)

mfuzz_class <- new('ExpressionSet',exprs = AS_LEVEL)

mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

mfuzz_class <- standardise(mfuzz_class)

set.seed(314)
cluster_num <- 6

mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

#作图，详情 ?mfuzz.plot2
#time.labels 参数设置时间轴，需要和原基因表达数据集中的列对应
#颜色、线宽、坐标轴、字体等细节也可以添加其他参数调整，此处略，详见函数帮助
color.2 <- colorRampPalette(rev(c(   "#C71000FF","white","#008EA0FF")))(100)
label_mfuzz <- gsub("_ASlevel","",gsub("N_S_","",colnames(AS_LEVEL)))


mat <- matrix( c(1,2,3,7,4,5,6,0),ncol=4,nrow=2,byrow=TRUE)
#par(oma = c(2,2,2,2))
l   <- layout(mat,widths = c(2,2,2,0.5), heights = c(1, 1))


for (i in c(1:cluster_num)) {
  mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster,
              colo = color.2, 
              mfrow = NA, 
              time.labels = label_mfuzz,
              centre=TRUE,
              x11=F,
              col.main="black",
              #centre.col="yellow",
              ylab="ASlevel changes",
              #colo="red",
              single = i,
              centre.lwd=3,cex.main = 1.5,cex.lab = 1.5
              #plot.main = "asadsa"
  )
  
}
mfuzzColorBar(col=color.2,main="Mvlaue",cex.main=1,horizontal = F)




cluster_size <- mfuzz_cluster$size

names(cluster_size) <- 1:cluster_num


AS_cluster <- mfuzz_cluster$cluster
AS_cluster <- cbind(AS_LEVEL[names(AS_cluster), ], AS_cluster)
head(AS_cluster)
write.table(AS_cluster, './result/02_AS_cluster_nona.txt', sep = '\t', col.names = NA, quote = FALSE)

AS_cluster_df <- data.frame(AS_cluster)
AS_cluster_1 <- AS_cluster_df
[
  AS_cluster_df["AS_cluster"]==1|AS_cluster_df["AS_cluster"]==2,#|AS_cluster_df["AS_cluster"]==5|AS_cluster_df["AS_cluster"]==6,
  ]



#RI_all["max"] <- apply(RI_all[,7:50],1,max)

#RI_all_1 <- RI_all[RI_all["min"]>0,] #%>% na.omit()
RI_AS <- data.frame(AS = rownames(AS_cluster_1))
#RI_AS["AS"] <- rownames(AS_cluster_1)
RI_AS <- separate(data = RI_AS,col = "AS",c("AS_Gene_ID","AS_index"),sep = ".v6.1_")
RI_AS$`gene` <- gsub("RI_","",RI_AS$AS_Gene_ID)
RI_AS$`gene` <- gsub("A3SS_","",RI_AS$`gene`)
RI_AS$`gene` <- gsub("A5SS_","",RI_AS$`gene`)
RI_AS$`gene` <- gsub("MXE_","",RI_AS$`gene`)
RI_AS$`gene` <- gsub("SE_","",RI_AS$`gene`)

write.table(RI_AS,file = "./tmp/genelist.txt",sep = "\t",
            quote = F, col.names = T,row.names = F)


###############################################################################

#03 可变事件统计柱形图
library(ggplot2)
library(readxl)
library(reshape2)
library(readxl)
library(paletteer)
library(ggsci)
#AS_statistics <- read_excel("./AS_事件检测结果统计.xls",sheet = "表3 以IJC>=2&&SJC>=2&&Skiplevel>=")
AS_statistics <- read_excel("./result/01_AS_gene_GO_enric.xls",sheet = "Sheet2")
#AS_statistics <- t(AS_statistics)
AS_statistics <- as.data.frame(AS_statistics)
#AS_statistics$`样本` <- c("N_S_0h","N_S_10m","N_S_30m","N_S_1h","N_S_2h",
#                        "N_S_6h","N_S_8h","N_S_24h","N_S_48h"
#)
AS_statistics$`样本` <- c("0h","10m","30m","1h","2h",
                        "6h","8h","24h","48h"
)
colnames(AS_statistics)[1] <- "AS_event"
AS_statistics_mat = melt(AS_statistics,variable.name = "Sample",value.name = "Numbers of AS events")
colnames(AS_statistics_mat) <- c("Sample","AS_event","Numbers of AS events")


AS_statistics_mat_order = AS_statistics_mat[sort(c(1:45),decreasing = T),]
rownames(AS_statistics_mat_order) <- NULL


AS_statistics_mat_order$Sample <- factor(AS_statistics_mat_order$Sample,
                                         levels = c("0h","10m","30m",
                                                    "1h","2h","6h",
                                                    "8h","24h","48h"))

AS_statistics_mat_order$AS_event <- factor(AS_statistics_mat_order$AS_event,
                                           levels = c("RI","SE","A3SS","A5SS","MXE"))

AS_statistics_plot <- ggplot(AS_statistics_mat_order, aes(x = AS_event,y = `Numbers of AS events`,fill = Sample,Sample=1))+
  geom_bar(position = 'dodge',stat = 'identity',width = 0.9,colour="black")+
  ##scale_fill_manual(values =c( "#A4A4A4","#8E8E8E", "#696969","#2F2F2F","#232323" ))  +   
  #########设定颜色
  #scale_fill_discrete(labels=c("SE","A3SS","A5SS","MXE","RI5")) +
  #scale_fill_manual(values =pal_futurama("planetexpress", alpha = 0.9)(5))  +
  scale_fill_manual(labels = c("0h","10m","30m",
                               "1h","2h","6h",
                               "8h","24h","48h"), 
                    values = pal_futurama("planetexpress", alpha = 0.5)(9)) +
  # scale_fill_brewer( values=pal_futurama("planetexpress", alpha = 0.9)(5), labels = c("SE","A3SS","A5SS","MXE","RI5")) + 
  guides(fill = guide_legend(reverse = F))+
  xlab("Alternative splicing type")+
  ylab("NO. of the AS events")+
  labs(fill = " ")+
  coord_cartesian(clip="off")+
  theme_bw(base_size = 25)+
  #theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        legend.position = 'right',
        axis.text = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(angle = 0,size = 15,face = "bold",hjust = 1),
        axis.title = element_text(size = 18,face = "bold"),
        legend.text = element_text(size = 13,face = "bold"),
        legend.title = element_text(size = 15,face = "bold"),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm')
  )+
  geom_text(aes(label=`Numbers of AS events`),
            angle = 90,
            position = position_dodge(0.9))



AS_statistics_plot

