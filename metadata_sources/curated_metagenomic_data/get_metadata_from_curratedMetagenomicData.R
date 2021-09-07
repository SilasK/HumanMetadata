# BiocManager::install("curatedMetagenomicData")

library(curatedMetagenomicData)
#dataset= curatedMetagenomicData("LeChatelierE_2013.marker_abundance.stool",dryrun = F)

L= combined_metadata@listData


n <- length(L[[1]])
DF <- structure(L, row.names = c(NA, -n), class = "data.frame")

write.table(DF,sep='\t',file = 'curatedMetagenomicData.tsv')


table(DF$dataset_name, !is.na(DF$BMI))


