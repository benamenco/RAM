read.OTU <- function(file, sep=",") {
  
  input <- read.table(file, header=TRUE, row.names=1, sep=sep)
  valid.OTU(input)
  
  input$taxonomy <- as.character(input$taxonomy)
  
  return(input)
}

read.meta <- function(file, sep=",") {
  
  input <- read.table(file, header=TRUE, row.names=1, sep=sep)
  ### PUT META DATA VALIDATION HERE
  
  
  return(input)
}