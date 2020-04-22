## ----setup_ggplot2, include=FALSE-------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(rtracklayer)
require(Rsubread)
require(profileplyr)
require(magrittr)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(rtracklayer)
require(rio)
require(ggplot2)
require(Rsamtools)
AsSlides <- TRUE
TESTORNOT <- FALSE


## ----setwd_introtoR,eval=F--------------------------------------------------------------------------------------
## setwd("/PathToMyDownload/GenomicHeatmapsAndProfiles/r_course")
## ## For Fuchs Lab
## # e.g. setwd("/Workshop_Data/TestFolder/")


## ----GetFasta,cache=TRUE,warning=FALSE,message=FALSE,eval=FALSE-------------------------------------------------
## require(BSgenome.Hsapiens.UCSC.hg19)
## mainChromosomes <- paste0("chr",c(1:22,"X","Y","M"))
## mainChrSeq <- lapply(mainChromosomes,
##                      function(x)BSgenome.Hsapiens.UCSC.hg19[[x]])
## names(mainChrSeq) <- mainChromosomes
## mainChrSeqSet <- DNAStringSet(mainChrSeq)
## writeXStringSet(mainChrSeqSet,
##                 "BSgenome.Hsapiens.UCSC.hg19_majorChr.fa")


## ----index,cache=TRUE,warning=FALSE,message=FALSE,eval=FALSE----------------------------------------------------
## require(Rsubread)
## buildindex("BSgenome.Hsapiens.UCSC.hg19",
##            "BSgenome.Hsapiens.UCSC.hg19_majorChr.fa",
##            indexSplit = TRUE)


## ----align,cache=TRUE,warning=FALSE,message=FALSE,eval=FALSE----------------------------------------------------
## align("BSgenome.Hsapiens.UCSC.hg19",
##       "ENCFF001NQP.fastq.gz",
##       output_file="ENCFF001NQP.bam",
##       type="dna")


## ----countw,cache=TRUE,warning=FALSE,message=FALSE,eval=FALSE---------------------------------------------------
## require(Rsamtools)
## sortBam("ENCFF001NQP.bam","Sorted_ENCFF001NQP")
## indexBam("Sorted_ENCFF001NQP.bam")


## ----countw2a,cache=TRUE,warning=FALSE,message=FALSE,eval=FALSE,echo=TRUE---------------------------------------
## MappedReads <- idxstatsBam("Sorted_ENCFF001NQP.bam")
## TotalMappedReads <- sum(MappedReads$mapped)


## ----countw2,cache=TRUE,warning=FALSE,message=FALSE,eval=TRUE,echo=FALSE----------------------------------------

MappedReads <- idxstatsBam("~/Desktop/Projects/brc/renvTest/AnotherTest2/BAMs/Sorted_shPTBP153_1.bam")
TotalMappedReads <- 10302212
TotalMappedReads


## ----bw,cache=TRUE,warning=FALSE,message=FALSE,eval=FALSE-------------------------------------------------------
## require(rtracklayer)
## coverageOverGenome <- coverage("Sorted_ENCFF001NQP.bam",
##                                weight = (10^6)/TotalMappedReads)
## export.bw(coverageOverGenome,con="Sorted_ENCFF001NQP.bw")


## ---- gettingData,eval=FALSE------------------------------------------------------------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## BiocManager::install("GEOquery")
## 


## ----download_Data1,cache=TRUE,warning=FALSE,message=FALSE,echo=TRUE,eval=FALSE---------------------------------
## require(GEOquery)
## 
## GSM4332935_BigWig <- getGEOSuppFiles("GSM4332935",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332934_BigWig <- getGEOSuppFiles("GSM4332934",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332944_BigWig <- getGEOSuppFiles("GSM4332944",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332940_BigWig <- getGEOSuppFiles("GSM4332940",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332950_BigWig <- getGEOSuppFiles("GSM4332950",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")


## ----download_Data2,cache=TRUE,warning=FALSE,message=FALSE,echo=FALSE,eval=FALSE--------------------------------
## require(GEOquery)
## 
## GSM4332935_BigWig <- getGEOSuppFiles("GSM4332935",fetch_files = FALSE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332934_BigWig <- getGEOSuppFiles("GSM4332934",fetch_files = FALSE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332944_BigWig <- getGEOSuppFiles("GSM4332944",fetch_files = FALSE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332940_BigWig <- getGEOSuppFiles("GSM4332940",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")
## GSM4332950_BigWig <- getGEOSuppFiles("GSM4332950",fetch_files = TRUE,
##                                      makeDirectory = FALSE,
##                                      filter_regex = "*.bw")


## ----echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE------------------------------------------------------------
require(rtracklayer)
ZBTB1_GR <- import.bed("Beds/ZBTB1Peaks.bed")
ZBTB1_GR


## ----echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE------------------------------------------------------------
ZBTB1_GR


## ----echo=FALSE,eval=TRUE---------------------------------------------------------------------------------------
promoterPos <- TxDb.Hsapiens.UCSC.hg19.knownGene %>% genes %>% promoters
rtracklayer::export(promoterPos,con="Beds/Promoter_Positions.bed")


## ----echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE------------------------------------------------------------
require(rtracklayer)
promoter_positions <- import.bed("Beds/Promoter_Positions.bed")
promoter_positions


## ----echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE------------------------------------------------------------
ZBTB1_In_Promoters <- ZBTB1_GR[ZBTB1_GR %over% promoter_positions]
ZBTB1_Outside_Promoters <- ZBTB1_GR[!ZBTB1_GR %over% promoter_positions]
export.bed(ZBTB1_In_Promoters,
           con="Beds/ZBTB1_In_Promoters.bed")



## computeMatrix reference-point --referencePoint center

## -bs binSize -a BP_after_BED -b BP_before_BED

## -S  BIGWIGSofINTEREST

## --regionsFileName BEDFILEOFINTEREST

## --outFileName OUTFILE

## 

## computeMatrix reference-point --referencePoint center

## -bs 50 -a 2000 -b 2000 -p 4

## -S  GSM4332935_Sorted_GFP_MinusN_FLAGNormalised.bw GSM4332940_Sorted_GFP_PlusN_FLAGNormalised.bw GSM4332945_Sorted_ZBTB_MinusN_FLAGNormalised.bw GSM4332950_Sorted_ZBTB_PlusN_FLAGNormalised.bw

## --regionsFileName ZBTB1Peaks.bed

## --outFileName Deeptools_ZBTB1Peaks.MAT


## plotHeatmap -m MATRIXFILE

## -o HEATMAP

## --colorList lowColour,highColour


## plotHeatmap -m Deeptools_ZBTB1Peaks.MAT

## -o Deeptools_ZBTB1Peaks.MAT.png

## --colorList white,darkred


## computeMatrix reference-point --referencePoint center -bs 50 -a 2000 -b 2000 -p 4 -S  /Users/thomascarroll/Desktop/software/Github/Genomic_HeatmapsAndProfiles//ppBw//GSM4332935_Sorted_GFP_MinusN_FLAGNormalised.bw /Users/thomascarroll/Desktop/software/Github/Genomic_HeatmapsAndProfiles//ppBw//GSM4332940_Sorted_GFP_PlusN_FLAGNormalised.bw /Users/thomascarroll/Desktop/software/Github/Genomic_HeatmapsAndProfiles//ppBw//GSM4332945_Sorted_ZBTB_MinusN_FLAGNormalised.bw /Users/thomascarroll/Desktop/software/Github/Genomic_HeatmapsAndProfiles//ppBw//GSM4332950_Sorted_ZBTB_PlusN_FLAGNormalised.bw --regionsFileName //Users/thomascarroll/Desktop/software/Github/Genomic_HeatmapsAndProfiles/ZBTB1Peaks.bed --outFileName /Users/thomascarroll/Desktop/software/Github/Genomic_HeatmapsAndProfiles/Deeptools_ZBTB1Peaks.MAT


## ----zbtb1_1,echo=TRUE,eval=TRUE,cache=TRUE,warning=FALSE,message=FALSE-----------------------------------------
require(profileplyr)
bigWigs <- dir(full.names=TRUE,pattern = "FLAG")
zbtb1_profile <- BamBigwig_to_chipProfile(bigWigs,
                                   testRanges = "Beds/ZBTB1Peaks.bed",
                                   style = "point",
                                   format = "bigwig",
                                   distanceAround = 2000)

zbtb1_profile <- as_profileplyr(zbtb1_profile)
zbtb1_profile


## ----"zbtb1_1plot",eval=TRUE,cache=TRUE,dependson="zbtb1_1",fig.width = 10,fig.height = 5-----------------------
generateEnrichedHeatmap(zbtb1_profile)


## ----"zbtb1_1Deeptools",eval=TRUE,cache=TRUE,fig.width = 10,fig.height = 5--------------------------------------
zbtb1_deeptools <- import_deepToolsMat("Deeptools_ZBTB1Peaks.MAT")
generateEnrichedHeatmap(zbtb1_deeptools)


## ---------------------------------------------------------------------------------------------------------------
blacklist <- rtracklayer::import("Beds/ENCFF001TDO.bed.gz")
blacklist


## ----grouped,echo=TRUE,eval=TRUE,cache=TRUE,fig.width = 10,fig.height = 5,message=FALSE,warning=FALSE-----------
zbtb1_profile_Grouped <- groupBy(zbtb1_profile,
                                 group = blacklist,
                                 include_nonoverlapping = TRUE)
generateEnrichedHeatmap(zbtb1_profile_Grouped)


## ----bled,echo=TRUE,eval=TRUE,cache=TRUE,fig.width = 10,fig.height = 5------------------------------------------
zbtb1_profile_bl <- zbtb1_profile[!zbtb1_profile %over% blacklist]
generateEnrichedHeatmap(zbtb1_profile_bl)


## ----echo=TRUE,eval=FALSE,dependson="bled"----------------------------------------------------------------------
## export_deepToolsMat(zbtb1_profile_bl,
##                     con="zbtb1_Profileplyr.MAT",
##                     overwrite = TRUE)


## 
## plotHeatmap -m zbtb1_Profileplyr.MAT.gz

##             -o zbtb1_Profileplyr.MAT.png

##             --colorList white,darkred


## ----zbtb_Sub,echo=TRUE,eval=TRUE,cache=TRUE,dependson="plotPromoters",fig.width = 10,fig.height = 5------------
zbtb1_profile_subset <- zbtb1_profile_bl[,,3:4]
generateEnrichedHeatmap(zbtb1_profile_subset)



## ----ATF4,echo=TRUE,eval=TRUE,cache=TRUE------------------------------------------------------------------------
bigWigs <- dir(full.names=TRUE,pattern = "ATF4")
atf4_profile <- BamBigwig_to_chipProfile(bigWigs,
                                   testRanges = "Beds/ZBTB1Peaks.bed",
                                   style = "point",
                                   format = "bigwig",distanceAround = 2000)
atf4_profile <- as_profileplyr(atf4_profile)
atf4_profile_bl <- atf4_profile[!atf4_profile %over% blacklist]


## ----JoinAndPlot,echo=TRUE,eval=TRUE,cache=TRUE,dependson="ATF4",fig.width = 10,fig.height = 5------------------
ztbtb1Andatf4_profile <- c(zbtb1_profile_subset,atf4_profile_bl)
generateEnrichedHeatmap(ztbtb1Andatf4_profile)


## ----sampleUpdate,echo=TRUE,eval=TRUE,cache=TRUE,dependson="ATF4",fig.width = 10,fig.height = 5-----------------
rownames(sampleData(ztbtb1Andatf4_profile)) <- c("ZBTB1_1","ZBTB1_2","ATF4_1")
generateEnrichedHeatmap(ztbtb1Andatf4_profile)


## ----sampleUpdate22,echo=TRUE,eval=TRUE,cache=TRUE,dependson="sampleUpdate",fig.width = 10,fig.height = 5-------
generateEnrichedHeatmap(ztbtb1Andatf4_profile,
                        sample_names=c("zbtb_A","zbtb_B","atf4_A"))


## ----sampleUpdate2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="sampleUpdate",fig.width = 10,fig.height = 5--------
sampleData(ztbtb1Andatf4_profile)$Antibody <- c("ZBTB1","ZBTB1","ATF4")
sampleData(ztbtb1Andatf4_profile)$Grade <- c("Favourite","Good","Good")


## ----sampleUpdatePlot,echo=TRUE,eval=TRUE,cache=TRUE,dependson="sampleUpdate2",fig.width = 10,fig.height = 5----
generateEnrichedHeatmap(ztbtb1Andatf4_profile,
                        color_by_sample_group = "Antibody")


## ----clusterAndPlot,echo=TRUE,eval=TRUE,cache=TRUE,dependson="sampleUpdate"-------------------------------------
ztbtb1Andatf4_profile <- clusterRanges(ztbtb1Andatf4_profile,
                        fun = rowSums,
                        cutree_rows = 2)


## ----clusterAndPlot2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="clusterAndPlot",fig.width = 10,fig.height = 5----
generateEnrichedHeatmap(ztbtb1Andatf4_profile,
                        color_by_sample_group = "Antibody")


## ----groupAndPlot,echo=TRUE,eval=TRUE,cache=TRUE,dependson=c("sampleUpdate","makeBed2")-------------------------
ATF4ANDZBTB1_Peaks <- import.bed("Beds/ATF4ANDZBTB1_Within2000kb.bed")
ztbtb1Andatf4_profile <- groupBy(ztbtb1Andatf4_profile,
                                 group = ATF4ANDZBTB1_Peaks,
                                 GRanges_names = "ProximalPeaks",
                                 include_nonoverlapping = TRUE)


## ----groupAndPlot2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 5--------
generateEnrichedHeatmap(ztbtb1Andatf4_profile,
                        color_by_sample_group = "Antibody",
                        extra_annotation_columns = "cluster")


## ----summarise,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot"------------------------------------------
ztbtb1Andatf4_summarized <- summarize(ztbtb1Andatf4_profile,
                        rowSums,
                        output="long",
                       sampleData_columns_for_longPlot = "Antibody")
ztbtb1Andatf4_summarized[1:4,]


## ----summarisePlot,echo=TRUE,eval=TRUE,cache=TRUE,dependson="summarise",fig.width = 10,fig.height = 4.5---------
require(ggplot2)
ggplot(ztbtb1Andatf4_summarized,
       aes(y=Signal,x=Sample,fill=GR_overlap_names))+
  facet_wrap(Antibody~.,scales = "free")+
  geom_violin()+scale_y_log10()


## ----annotate,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 4.5-----------
ztbtb1Andatf4_profile <- annotateRanges_great(ztbtb1Andatf4_profile,
                                              species = "hg19")


## ----annotated,echo=TRUE,eval=TRUE,cache=TRUE,dependson="annotate"----------------------------------------------
annotatedRanges <- rowRanges(ztbtb1Andatf4_profile)
annotatedRanges[1,]


## ----annotated2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="annotated"--------------------------------------------
write.table(annotatedRanges,"annotatedRanges.csv",sep=",",row.names=FALSE)
export.bed(annotatedRanges,"annotatedRanges.bed")


## ----annotated3,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 5-----------
generateEnrichedHeatmap(ztbtb1Andatf4_profile,
                        color_by_sample_group = "Antibody",
                        extra_annotation_columns = "cluster",
                        genes_to_label = "ASNS",
                        gene_label_font_size = 18)


## ----annotated4,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 5,message=FALSE,warning=FALSE----
require(GSEABase)
aa_meta <- getGmt("GeneSets/aminoacid_metabolism.gmt")
aa_meta_Genes <- geneIds(aa_meta)


## ----annotated5,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 5,message=FALSE,warning=FALSE----
ztbtb1Andatf4_profile2 <- groupBy(ztbtb1Andatf4_profile,
                                 group = aa_meta_Genes,
                                 include_nonoverlapping = TRUE)



## ----annotated8,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 5,message=FALSE,warning=FALSE----

generateEnrichedHeatmap(ztbtb1Andatf4_profile2,
                        color_by_sample_group = "Antibody",
                        genes_to_label = "ASNS",
                        gene_label_font_size = 18)


## ----annotated9,echo=TRUE,eval=TRUE,cache=TRUE,dependson="groupAndPlot",fig.width = 10,fig.height = 5,message=FALSE,warning=FALSE----
ztbtb1Andatf4_profile2 <- groupBy(ztbtb1Andatf4_profile2,group="cluster")
generateEnrichedHeatmap(ztbtb1Andatf4_profile2,
                        color_by_sample_group = "Antibody",genes_to_label = "ASNS",
                        extra_annotation_columns = "GL_overlap_names",
                        gene_label_font_size = 18)

