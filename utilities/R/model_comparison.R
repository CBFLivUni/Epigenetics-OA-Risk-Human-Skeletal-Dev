

GetDiffVsBaseline <- function(data, baselineVals) {
  
  # Create DELTA IC
  data <- append(data, list(DELTA_IC=data.frame(data$IC)))
   
  # Calculate gene-wise differences in AIC, relative to model 1
  for (i in 1:ncol(data$IC)) {
    
    # Calculate difference
    data$DELTA_IC[,i] <- data$DELTA_IC[,i] - baselineVals
    
  }
  
  # Return
  return(data)
  
}