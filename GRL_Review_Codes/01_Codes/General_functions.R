# This codes include general functions for additional analysis

# This function calculates daily mean variable
Daily_mean <- function(data){
  DM <- matrix(data=data,nrow=window_size)
  DM <- colMeans(DM,na.rm=T)
  return(DM)
}

# Get required variable after QC
# Input is the variable name
Var_QC <- function(variable){
  varlist <- colnames(AMF)
  # Variable name in the dataset
  var <- varlist[grepl(variable,varlist)&!grepl("QC",varlist)]
  var <- AMF[var][,1]
  # QC for this variable
  var_QC <- varlist[grepl(variable,varlist)&grepl("QC",varlist)]
  if(length(var_QC!=0)){
    # Apply QC, only keep QC = 0
    var_QC <- AMF[var_QC][,1]
    var[var_QC!=0] <- NA
  }
  return(var)
}

