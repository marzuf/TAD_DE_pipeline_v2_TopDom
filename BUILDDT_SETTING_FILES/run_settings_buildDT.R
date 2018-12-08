
pipOutFold <- "BUILDDT_OUTPUT_FOLDER"



inputFolder <- "OUTPUT_FOLDER"

# list the folder to use for the cross-comparison
# ls -1 -Q <folder>
# cannot run EZH2 data with the other because not the same gene list files !

comp_folders <- list.files(inputFolder, full.names=TRUE)
#comp_folders <- comp_folders[1:3]
#comp_folders <- comp_folders[1:3]

#comp_folders <- comp_folders[21:22]
#comp_folders <- comp_folders[23:length(comp_folders)]

#comp_folders <- c(
#  "GSE40419_normal_cancer",
#  "GSE73765_noninf_list",
#  "TCGAbrca_lum_bas",
#  "TCGAstad_msi_gs",
#  "TCGAucec_msi_cnl"
#)


#comp_folders <- c(
#  "GSE40419_normal_cancer",
##  "GSE73765_noninf_list",
#  "TCGAbrca_lum_bas")
#  "TCGAstad_msi_gs",
#  "TCGAucec_msi_cnl"
#)

#comp_folders <- file.path("OUTPUT_FOLDER", comp_folders)

upToRank <- 100

#allDown <- c("ratioDown", "FCdown", "meanConcordRatioFC" , "prodConcord", "prodMeanConcord", "prodSignedRatio")
#toPlotRanked <- c("ratioDown", "FCdown", "meanConcordRatioFC" , "prodConcord", "prodMeanConcord", "prodSignedRatio")


toPlotRanked <- c("ratioDown", "prodSignedRatio")

toPlotLolli <- c("ratioDown", "prodSignedRatio")
nTopLolli <- 10


