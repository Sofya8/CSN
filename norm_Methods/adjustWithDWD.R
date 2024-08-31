adjustWithDWD <- function(fstTableName, secTableName){
  newFstTableName <- gsub("pre_norm", "after_DWD", fstTableName)
  newSecTableName <- gsub("pre_norm", "after_DWD", secTableName)
  
  fstTable <- read.delim(file=fstTableName, header=TRUE, sep=",")
  #row.names(fstTable)=fstTable[,'tracking_id'] #the genes col name is tracking_id. change tracking_id to be the row names
  numSamplesFstTable <- ncol(fstTable)
  
  secTable <- read.delim(file=secTableName, header=TRUE, sep=",")
  #row.names(secTable)=secTable[,'tracking_id'] #the genes col name is tracking_id. change tracking_id to be the row names
  numSamplesSecTable <- ncol(secTable) 
  
  dwdInput1 <- fstTable[,2:numSamplesFstTable]
  dwdInput2 <- secTable[,2:numSamplesSecTable]
  dwd.output <- dwd(dwdInput1, dwdInput2, skip.match=TRUE)
  newFstTable <- fstTable
  newSecTable <- secTable
  newFstTable[,2:numSamplesFstTable] <- dwd.output$x
  newSecTable[,2:numSamplesSecTable] <- dwd.output$y
  write.table(newFstTable, file = newFstTableName, sep ="\t", dec =".", row.names = FALSE, col.names = TRUE)
  write.table(newSecTable, file = newSecTableName, sep ="\t", dec =".", row.names = FALSE, col.names = TRUE)
}