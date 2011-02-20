##
##  PURPOSE:   Create a list with individual longitudinal profiles of a given variable
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   07/05/2008 (as a stand alone function)
##             26/11/2010:  added to the mixAK package
##
##  FUNCTIONS: getProfiles
##
## ==========================================================================

## *************************************************************
## getProfiles
## *************************************************************
##
getProfiles <- function(t, y, id, data)
{

## It is assumed that values for one subject follow each other
  
  T <- data[,t]
  ID <- data[,id]

  tabID <- table(ID)
  nID <- length(tabID)
  cumID <- c(0, cumsum(tabID))

  RET <- list()
  for (i in 1:nID){
    RET[[i]] <- data.frame(T[(cumID[i]+1):cumID[i+1]])
    colnames(RET[[i]]) <- t
    RET[[i]] <- cbind(RET[[i]], data[(cumID[i]+1):cumID[i+1], y])
  }
  names(RET) <- names(tabID)
  return(RET)
}
