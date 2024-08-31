adjustDatasets <- function(fstTableName, secTableName){
  newFstTableName <- gsub("pre_norm", "after_EB", fstTableName)
  newSecTableName <- gsub("pre_norm", "after_EB", secTableName)
  fstTable <- read.delim(file=fstTableName, header=TRUE, sep=",")
  #row.names(fstTable)=fstTable[,'tracking_id'] #the genes col name is tracking_id. change tracking_id to be the row names
  numSamplesFstTable <- ncol(fstTable) - 1

  secTable <- read.delim(file=secTableName, header=TRUE, sep=",")
  #row.names(secTable)=secTable[,'tracking_id'] #the genes col name is tracking_id. change tracking_id to be the row names
  numSamplesSecTable <- ncol(secTable) - 1 
  mergedTables <- merge(fstTable,secTable,by="tracking_id")
  row.names(mergedTables)=mergedTables[,'tracking_id'] #the genes col name is tracking_id. change tracking_id to be the row names
  mergedTables <- subset(mergedTables, select = -c(tracking_id)) #remove the col tracking_id
  filteredMergedTables <- mergedTables [which(rowSums(mergedTables) > 0), ]
  mergedTablesMatrix <- as.matrix(sapply(filteredMergedTables, as.numeric))
  group <- c(rep(1, numSamplesFstTable), rep(2, numSamplesSecTable))
  adjustedMatrix <- ComBat(mergedTablesMatrix, batch=group)
  mergedTables [which(rowSums(mergedTables) != 0), ] <- adjustedMatrix
  newFstTable <- mergedTables[,1:numSamplesFstTable]
  newSecTable <- mergedTables[,(numSamplesFstTable+1):ncol(mergedTables)]
  newFstTable <- cbind(tracking_id = rownames(newFstTable), newFstTable)
  write.table(newFstTable, file = newFstTableName, sep ="\t", dec =".", row.names = FALSE, col.names = TRUE)
  newSecTable <- cbind(tracking_id = rownames(newSecTable), newSecTable)
  write.table(newSecTable, file = newSecTableName, sep ="\t", dec =".", row.names = FALSE, col.names = TRUE)
}