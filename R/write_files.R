write.OTU <- function(data, file) {
  
  file <- .ensure.filepath(file, ".csv")
  
  write.csv(data, file=file, quote=FALSE)
}