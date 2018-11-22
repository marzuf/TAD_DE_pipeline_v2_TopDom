SSHFS <- TRUE
setDir <- ifelse(SSHFS, "/media/electron", "")

# file with assignment from entrez to all regions
gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 

gene2tadDT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
head(gene2tadDT)

pipDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom")

#entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2position/filter_entrez_map.Rdata") # previous also held symbols
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, stringsAsFactors = F, header=T)
head(entrezDT)

outFold <- file.path(pipDir, "LOOK_COORD_TADs")
system(paste0("mkdir -p ", outFold))


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- ifelse(plotType == "png", 300, 7)

curr_ds <- "TCGAcrc_msi_mss"

#################################################################################################################################
################################################################################################################################# find some interesting TADs ?
#################################################################################################################################


nucleoGenesEntrez <- setNames(entrezDT$entrezID[grepl("^NAP", entrezDT$symbol)], 
                       entrezDT$symbol[grepl("^NAP", entrezDT$symbol)])

nucleoTADs <- setNames(gene2tadDT$region[gene2tadDT$entrezID %in% nucleoGenesEntrez], 
                       gene2tadDT$entrezID[gene2tadDT$entrezID %in% nucleoGenesEntrez])

table(nucleoTADs)
nucleoTADs
# chr10_TAD120 chr11_TAD179   chr11_TAD5 chr12_TAD109  chr13_TAD54  chr18_TAD17  chr19_TAD60  chr19_TAD66  chr1_TAD190  chr20_TAD34  chr3_TAD167 
# 1            1            1            1            1            1            2            2            1            1            1 
# chr4_TAD114   chr6_TAD56 chr7_BOUND41 chr8_BOUND49 chrX_BOUND36 chrX_BOUND38   chrX_TAD72 
# 1            1            1            1            1            1            2 


nucleoGenesEntrez <- setNames(entrezDT$entrezID[grepl("^HMGN", entrezDT$symbol)], 
                              entrezDT$symbol[grepl("^HMGN", entrezDT$symbol)])

nucleoTADs <- setNames(gene2tadDT$region[gene2tadDT$entrezID %in% nucleoGenesEntrez], 
                       gene2tadDT$entrezID[gene2tadDT$entrezID %in% nucleoGenesEntrez])

table(nucleoTADs)
nucleoTADs
# chr10_TAD112  chr10_TAD136   chr10_TAD35   chr10_TAD98  chr11_TAD109   chr11_TAD23   chr11_TAD27   chr11_TAD64    chr11_TAD8   chr12_TAD26 
# 1             1             1             1             1             1             1             1             1             1 
# chr13_TAD101   chr13_TAD66 chr14_BOUND14 chr14_BOUND41    chr14_TAD4   chr14_TAD52   chr14_TAD61   chr14_TAD80    chr14_TAD9 chr15_BOUND23 
# 1             1             1             1             1             1             1             1             1             1 
# chr15_TAD10  chr15_TAD123   chr15_TAD40   chr15_TAD80   chr15_TAD96   chr16_TAD46   chr16_TAD60  chr17_TAD121   chr17_TAD70   chr17_TAD71 
# 1             1             1             1             1             1             1             1             1             1 
# chr17_TAD93   chr17_TAD94   chr18_TAD42   chr18_TAD75   chr18_TAD80 chr19_BOUND39    chr1_TAD15   chr1_TAD207   chr1_TAD212   chr1_TAD220 
# 1             1             1             1             1             1             1             1             1             1 
# chr1_TAD269   chr1_TAD335    chr1_TAD54   chr21_TAD20   chr21_TAD37   chr22_TAD44   chr22_TAD48   chr2_TAD106   chr2_TAD138   chr2_TAD144 
# 1             1             1             1             1             1             1             1             1             1 
# chr2_TAD164   chr2_TAD278    chr2_TAD36  chr3_BOUND11   chr3_TAD137   chr3_TAD195   chr3_TAD198   chr3_TAD205   chr3_TAD221   chr3_TAD236 
# 1             1             1             1             1             1             1             1             1             2 
# chr3_TAD64  chr4_BOUND31   chr5_TAD122   chr5_TAD130   chr5_TAD132   chr5_TAD176    chr5_TAD54    chr5_TAD85   chr6_TAD118   chr6_TAD163 
# 1             1             2             1             1             1             1             1             2             1 
# chr6_TAD4    chr6_TAD46   chr7_TAD114   chr7_TAD166    chr7_TAD70    chr7_TAD76    chr9_TAD10   chr9_TAD115    chr9_TAD19    chr9_TAD77 
# 1             1             1             1             1             1             1             1             1             1 
# chrX_TAD135    chrX_TAD19    chrX_TAD66    chrX_TAD74    chrX_TAD76     chrX_TAD9 
# 1             1             1             1             1             1 
# 

olfactoGenesEntrez <- setNames(entrezDT$entrezID[grepl("^OR", entrezDT$symbol)], 
                              entrezDT$symbol[grepl("^OR", entrezDT$symbol)])

olfactoTADs <- setNames(gene2tadDT$region[gene2tadDT$entrezID %in% olfactoGenesEntrez], 
                       gene2tadDT$entrezID[gene2tadDT$entrezID %in% olfactoGenesEntrez])

olfactoTADs <- olfactoTADs[grep("_TAD", olfactoTADs)]
table(olfactoTADs)
olfactoTADs
# chr10_TAD188  chr10_TAD26  chr10_TAD68  chr11_TAD10 chr11_TAD105 chr11_TAD108  chr11_TAD11 chr11_TAD112  chr11_TAD12 chr11_TAD127  chr11_TAD13 
# 3            3            2           19            2            1           15            3            1            2            7 
# chr11_TAD14 chr11_TAD145  chr11_TAD15 chr11_TAD174 chr11_TAD175  chr11_TAD29  chr11_TAD32  chr11_TAD46   chr11_TAD6   chr11_TAD7  chr11_TAD71 
# 3            1            2           13           30            1            1            1            1            1           28 
# chr11_TAD72  chr11_TAD73  chr11_TAD75  chr11_TAD76  chr11_TAD77   chr11_TAD9  chr11_TAD96  chr12_TAD67  chr12_TAD68  chr12_TAD83  chr12_TAD84 
# 11          106           26            5            9           26            2            1            8            1           24 
# chr12_TAD85  chr13_TAD43  chr13_TAD69  chr13_TAD70  chr13_TAD71   chr14_TAD1   chr14_TAD2  chr14_TAD28   chr14_TAD3  chr14_TAD35   chr14_TAD4 
# 2            3            1            1            2           27            1            1            2            3            2 
# chr14_TAD9 chr15_TAD139   chr15_TAD3  chr16_TAD56  chr16_TAD61   chr16_TAD7  chr17_TAD10  chr17_TAD65  chr17_TAD88  chr19_TAD20  chr19_TAD21 
# 1            8            3            1            1            3           18            1            2           12            2 
# chr19_TAD73  chr1_TAD105  chr1_TAD177  chr1_TAD199  chr1_TAD220  chr1_TAD221  chr1_TAD361  chr1_TAD362  chr1_TAD363  chr21_TAD22  chr2_TAD106 
# 1            1            1            2           21           10            3           43            7            1            2 
# chr2_TAD107  chr2_TAD196  chr2_TAD209  chr2_TAD210  chr2_TAD268  chr2_TAD336  chr3_TAD124  chr3_TAD140  chr3_TAD155  chr3_TAD186    chr3_TAD7 
# 1            1            2            1            1            4            4           18            1            2            1 
# chr4_TAD16   chr4_TAD56   chr4_TAD99  chr5_TAD115  chr5_TAD229  chr6_TAD128  chr6_TAD176   chr6_TAD47   chr6_TAD48   chr6_TAD49   chr6_TAD50 
# 1            1            1            1            3            1            1            4            4            1            1 
# chr6_TAD51   chr7_TAD10  chr7_TAD129  chr7_TAD133  chr7_TAD142  chr7_TAD144   chr7_TAD15  chr7_TAD197  chr7_TAD199  chr7_TAD200  chr7_TAD201 
# 24            1            2            1            1            1            1            4            4            2           14 
# chr8_TAD12   chr8_TAD13   chr8_TAD26  chr9_TAD111  chr9_TAD121  chr9_TAD129  chr9_TAD140  chr9_TAD141   chr9_TAD44   chr9_TAD45   chr9_TAD88 
# 3            3            1           12            1            2           14            1            3            4            4 
# chrX_TAD107  chrX_TAD135 
# 9            1 


#################################################################################################################################
################################################################################################################################# retrieve dataset data
#################################################################################################################################

interestCol <- "firebrick1"

geneList <- eval(parse(text = load(file.path(pipDir, "OUTPUT_FOLDER", curr_ds, "0_prepGeneData/pipeline_geneList.Rdata"))))
regionList <- eval(parse(text = load(file.path(pipDir, "OUTPUT_FOLDER", curr_ds, "0_prepGeneData/pipeline_regionList.Rdata"))))

tads_meanLogFC <- eval(parse(text = load(file.path(pipDir, "OUTPUT_FOLDER", curr_ds, "3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))))

### INTRA TAD CORRELATION
tads_meanCorr <- eval(parse(text = load(file.path(pipDir, "OUTPUT_FOLDER", curr_ds, "4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))))
tads_meanCorr <- sort(tads_meanCorr, decreasing=T)
plot(x = tads_meanCorr, type ="h", col=ifelse(names(tads_meanCorr) %in% olfactoTADs, interestCol , "grey"))

### RATIO DOWN
tads_ratioDown <- eval(parse(text = load(file.path(pipDir, "OUTPUT_FOLDER", curr_ds, "8c_runAllDown/all_obs_ratioDown.Rdata"))))
tads_ratioDown <- sort(tads_ratioDown, decreasing=T)
plot(x = tads_ratioDown, type ="h", col=ifelse(names(tads_ratioDown) %in% olfactoTADs, interestCol , "grey"))








