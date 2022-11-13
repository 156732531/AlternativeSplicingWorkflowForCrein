library(IsoformSwitchAnalyzeR)


### Make design matrix

#样品信息
sample_list <- c('N_S_0h_rep1','N_S_0h_rep2','N_S_10m_rep1','N_S_10m_rep2',
                 'N_S_30m_rep1','N_S_30m_rep2','N_S_1h_rep1','N_S_1h_rep2',
                 'N_S_2h_rep1','N_S_2h_rep2','N_S_6h_rep1','N_S_6h_rep2',
                 'N_S_8h_rep1','N_S_8h_rep2','N_S_24h_rep1','N_S_24h_rep2',
                 'N_S_48h_rep1','N_S_48h_rep2')

### Import Salmon example data in R package
Stringtie_quant <- importIsoformExpression(
  parentDir = "./Quantitative/Stringtie/",
  readLength = 100
)


myDesign <- data.frame(
  sampleID = colnames(Stringtie_quant$abundance)[-1],
  condition = gsub('_rep\\d', '', colnames(Stringtie_quant$abundance)[-1])
)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = Stringtie_quant$counts,
  isoformRepExpression = Stringtie_quant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = './annotation/CreinhardtiiCC_4532_707_v6.1.gene_exons.gtf',#./raw_datas/06.Transcripts/stringtie/assembly.gtf',
  isoformNtFasta       = './annotation/CreinhardtiiCC_4532_707_v6.0.transcript.fa.tmp',
  showProgress = FALSE
)

###筛选
exampleSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 1,
  removeSingleIsoformGenes = TRUE
)

###Testing for Isoform Switches via DEXSeq
exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = exampleSwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)

extractSwitchSummary(exampleSwitchListAnalyzed)

colnames(exampleSwitchListAnalyzed$isoformFeatures)

### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
mySwitchList <- addORFfromGTF( exampleSwitchListAnalyzed,pathToGTF = './annotation/CreinhardtiiCC_4532_707_v6.1.gene_exons.gtf')
mySwitchList <- analyzeNovelIsoformORF( mySwitchList,analysisAllIsoformsWithoutORF = TRUE)

### Extracting Nucleotide and Amino Acid Sequences
mySwitchList <- extractSequence( mySwitchList )

### Add CPC2 analysis
exampleSwitchListAnalyzed <- analyzeCPC2(
  switchAnalyzeRlist   = mySwitchList,
  pathToCPC2resultFile = './raw_datas/06.Transcripts/CPC2/isoformSwitchAnalyze.txt',
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

### Add PFAM analysis
exampleSwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = exampleSwitchListAnalyzed,
  pathToPFAMresultFile = './raw_datas/06.Transcripts/Pfam/pfam_result',
  showProgress=FALSE
)

### Add SignalP analysis
exampleSwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = exampleSwitchListAnalyzed,
  pathToSignalPresultFile  = './raw_datas/06.Transcripts/signalP/output_protein_type.txt'
)




### This example relies on the example data from the 'Importing External Sequence Analysis' section above 

exampleSwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = exampleSwitchListAnalyzed,
  quiet=TRUE
)

# the consequences highlighted in the text above
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

exampleSwitchListAnalyzed <- analyzeIntronRetention(exampleSwitchListAnalyzed,
                                                    onlySwitchingGenes=TRUE,
                                                    alpha=0.05,
                                                    dIFcutoff = 0.1,
                                                    showProgress=TRUE,
                                                    quiet=FALSE)

exampleSwitchListAnalyzed_filtered <- analyzeSwitchConsequences(
  exampleSwitchListAnalyzed,
  #consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.1, # very high cutoff for fast runtimes - you should use the default (0.1)
  showProgress=FALSE
)
################################################################################
result <- exampleSwitchListAnalyzed$isoformFeatures

result <- dplyr::filter(result,dIF >= 0.1)

gene_result <- unique(result[,c(3,4)])
write.table(gene_result,file = "./tmp/tmp_gene.txt",quote = F,row.names = F,col.names = T,sep = "\t")
################################################################################

extractSplicingSummary(
  exampleSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE
)


b <- extractSplicingGenomeWide(
  exampleSwitchListAnalyzed,
  )

extractSwitchSummary(exampleSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)

exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = exampleSwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)
extractSwitchSummary(exampleSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

# 以0h为控制组
exampleSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)

extractSwitchSummary(exampleSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

extractTopSwitches(
  exampleSwitchListAnalyzedSubset, 
  filterForConsequences = TRUE, 
  n = 50, 
  sortByQvals = TRUE
)

switchingIso <- extractTopSwitches( 
  exampleSwitchListAnalyzedSubset, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

subset(switchingIso, gene_name == 'GLN2')

switchPlot(exampleSwitchListAnalyzedSubset, gene = 'GLN2')
switchPlotTranscript(exampleSwitchListAnalyzedSubset, gene = 'POA6')

###0h VS 10m
exampleSwitchListAnalyzed_10m <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_10m <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_10m, 
  exampleSwitchListAnalyzed_10m$isoformFeatures$condition_2 == 'N_S_10m'
)

exampleSwitchListAnalyzed_10m <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_10m, 
  subset =  abs(exampleSwitchListAnalyzed_10m$isoformFeatures$dIF)>=0.1)


exampleSwitchListAnalyzed_10m <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_10m, 
  subset =  exampleSwitchListAnalyzed_10m$isoformFeatures$switchConsequencesGene==TRUE)


exampleSwitchListAnalyzed_10m <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = exampleSwitchListAnalyzed_10m,
  quiet=TRUE
)

extractSwitchSummary(exampleSwitchListAnalyzed_10m, dIFcutoff = 0.1, filterForConsequences = TRUE)

a <- extractTopSwitches(
  exampleSwitchListAnalyzed_10m, 
  filterForConsequences = TRUE, 
  n = 50, 
  sortByQvals = TRUE
)



switchingIso_10m <- extractTopSwitches( 
  exampleSwitchListAnalyzed_10m, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

subset(switchingIso_24h, gene_name == 'GLN2')

switchPlot(exampleSwitchListAnalyzed_30m, gene = 'GLN2')
switchPlotTranscript(exampleSwitchListAnalyzed_10m, gene = 'FAP391')
extractSwitchSummary(exampleSwitchListAnalyzed_10m, dIFcutoff = 0.1, filterForConsequences = TRUE)

###0h vs 30m
exampleSwitchListAnalyzed_30m <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_30m <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_30m, 
  exampleSwitchListAnalyzed_30m$isoformFeatures$condition_2 == 'N_S_30m'
)
switchingIso_30m <- extractTopSwitches( 
  exampleSwitchListAnalyzed_30m, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

###0h vs 1h
exampleSwitchListAnalyzed_1h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_1h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_1h, 
  exampleSwitchListAnalyzed_1h$isoformFeatures$condition_2 == 'N_S_1h'
)
switchingIso_1h <- extractTopSwitches( 
  exampleSwitchListAnalyzed_1h, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

###0h vs 2h
exampleSwitchListAnalyzed_2h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_2h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_2h, 
  exampleSwitchListAnalyzed_2h$isoformFeatures$condition_2 == 'N_S_2h'
)
switchingIso_2h <- extractTopSwitches( 
  exampleSwitchListAnalyzed_2h, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

subset(switchingIso_10m, gene_id == 'Cre10.g433600.v5.5')

switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre06.g258733.v5.5')
switchPlotTranscript(exampleSwitchListAnalyzed_2h, gene = 'FAP391')
extractSwitchSummary(exampleSwitchListAnalyzed_2h, dIFcutoff = 0.1, filterForConsequences = TRUE)

###0h vs 6h
exampleSwitchListAnalyzed_6h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_6h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_6h, 
  exampleSwitchListAnalyzed_6h$isoformFeatures$condition_2 == 'N_S_6h'
)
switchingIso_6h <- extractTopSwitches( 
  exampleSwitchListAnalyzed_6h, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

###0h vs 8h
exampleSwitchListAnalyzed_8h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_8h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_8h, 
  exampleSwitchListAnalyzed_8h$isoformFeatures$condition_2 == 'N_S_8h'
)
switchingIso_8h <- extractTopSwitches( 
  exampleSwitchListAnalyzed_8h, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

###0h vs 24h
exampleSwitchListAnalyzed_24h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_24h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_24h, 
  exampleSwitchListAnalyzed_24h$isoformFeatures$condition_2 == 'N_S_24h'
)
switchingIso_24h <- extractTopSwitches( 
  exampleSwitchListAnalyzed_24h, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

###0h vs 48h
exampleSwitchListAnalyzed_48h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed, 
  exampleSwitchListAnalyzed$isoformFeatures$condition_1 == 'N_S_0h'
)
exampleSwitchListAnalyzed_48h <- subsetSwitchAnalyzeRlist(
  exampleSwitchListAnalyzed_48h, 
  exampleSwitchListAnalyzed_48h$isoformFeatures$condition_2 == 'N_S_48h'
)

switchingIso_48h <- extractTopSwitches( 
  exampleSwitchListAnalyzed_48h, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

AS_type <- as.data.frame(exampleSwitchListAnalyzed$AlternativeSplicingAnalysis)
AS_detail <- as.data.frame(exampleSwitchListAnalyzed$isoformFeatures)
AS_detail <- subset(AS_detail,condition_1 == "N_S_0h")
AS_detail <- subset(AS_detail,abs(dIF) >= 0.1)

write.csv(AS_type,file = "./组别AS/AS_type.csv")
write.csv(switchingIso,file = "./组别AS/switchingIso.csv")

write.csv(switchingIso_10m,file = "./组别AS/switchingIso_10m.csv")

write.csv(switchingIso_30m,file = "./组别AS/switchingIso_30m.csv")
write.csv(switchingIso_1h,file = "./组别AS/switchingIso_1h.csv")
write.csv(switchingIso_2h,file = "./组别AS/switchingIso_2h.csv")
write.csv(switchingIso_6h,file = "./组别AS/switchingIso_6h.csv")
write.csv(switchingIso_8h,file = "./组别AS/switchingIso_8h.csv")
write.csv(switchingIso_24h,file = "./组别AS/switchingIso_24h.csv")
write.csv(switchingIso_48h,file = "./组别AS/switchingIso_48h.csv")



###特定基因差异表达转录情况

## Cre02.g095126


switchPlot(exampleSwitchListAnalyzed_10m, gene = 'Cre02.g095126.v5.5')
switchPlot(exampleSwitchListAnalyzed_30m, gene = 'Cre02.g095126.v5.5')
switchPlot(exampleSwitchListAnalyzed_2h, gene = 'Cre02.g095126.v5.5')
switchPlot(exampleSwitchListAnalyzed_6hm, gene = 'Cre02.g095126.v5.5')

switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre02.g095126.v5.5') #

switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre02.g095126.v5.5') #




## Cre03.g149100

switchPlot(exampleSwitchListAnalyzed_10m, gene = 'Cre03.g149100.v5.5')
switchPlot(exampleSwitchListAnalyzed_30m, gene = 'Cre03.g149100.v5.5')
switchPlot(exampleSwitchListAnalyzed_1h, gene = 'Cre03.g149100.v5.5') #
switchPlot(exampleSwitchListAnalyzed_2h, gene = 'Cre03.g149100.v5.5')
switchPlot(exampleSwitchListAnalyzed_6h, gene = 'Cre03.g149100.v5.5')
switchPlot(exampleSwitchListAnalyzed_8h, gene = 'Cre03.g149100.v5.5') #
switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre03.g149100.v5.5')
switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre03.g149100.v5.5')

## Cre05.g232002 
switchPlot(exampleSwitchListAnalyzed_10m, gene = 'Cre05.g232002.v5.5') # 
switchPlot(exampleSwitchListAnalyzed_30m, gene = 'Cre05.g232002.v5.5') # 
switchPlot(exampleSwitchListAnalyzed_1h, gene = 'Cre05.g232002.v5.5')
switchPlot(exampleSwitchListAnalyzed_2h, gene = 'Cre05.g232002.v5.5')
switchPlot(exampleSwitchListAnalyzed_6h, gene = 'Cre05.g232002.v5.5')
switchPlot(exampleSwitchListAnalyzed_8h, gene = 'Cre05.g232002.v5.5') #
switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre05.g232002.v5.5')
switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre05.g232002.v5.5')

## Cre05.g241850##
switchPlot(exampleSwitchListAnalyzed_10m, gene = 'Cre05.g241850.v5.5')
switchPlot(exampleSwitchListAnalyzed_30m, gene = 'Cre05.g241850.v5.5')
switchPlot(exampleSwitchListAnalyzed_1h, gene = 'Cre05.g241850.v5.5') #

switchPlot(exampleSwitchListAnalyzed_2h, gene = 'Cre05.g241850.v5.5') #
switchPlot(exampleSwitchListAnalyzed_6h, gene = 'Cre05.g241850.v5.5') #
switchPlot(exampleSwitchListAnalyzed_8h, gene = 'Cre05.g241850.v5.5') #
switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre05.g241850.v5.5')
switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre05.g241850.v5.5')



## Cre08.g359350 ##
switchPlot(exampleSwitchListAnalyzed_10m, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_30m, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_1h, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_2h, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_6h, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_8h, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre08.g359350.v5.5')
switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre08.g359350.v5.5')


## Cre08.g378150
switchPlot(exampleSwitchListAnalyzed_1h, gene = 'Cre08.g378150.v5.5')

switchPlot(exampleSwitchListAnalyzed_6h, gene = 'Cre08.g378150.v5.5')

switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre08.g378150.v5.5')

switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre08.g378150.v5.5')

## Cre12.g530650##
switchPlot(exampleSwitchListAnalyzed_10m, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_30m, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_1h, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_2h, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_6h, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_8h, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_24h, gene = 'Cre12.g530650.v5.5')
switchPlot(exampleSwitchListAnalyzed_48h, gene = 'Cre12.g530650.v5.5')

switchPlotTranscript(exampleSwitchListAnalyzed_24h, gene = 'Cre08.g359350.v5.5')+
  labs(y ="条目名称",x = "基因数目")


