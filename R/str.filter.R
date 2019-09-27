#' @title str.filter
#' @description this function filters out individuals or loci with missing data according to a cut off percentage
#' @param structure.file a matrix or data frame object representing the structure-like genotype file
#' @param cutoff float. A percentage cutoff between 0 and 1 to exclude loci or individuals with missing data. Defaut NULL.
#' @param filter.loci logical. If TRUE loci with missing samples will be filtered out.
#' @param filter.samp logical. If TRUE samples with missing loci will be filtered out.
#' @author Marcelo Gehara
#' @export 
str.filter<-function(structure.file, cutoff=NULL, filter.loci=T, filter.samp=F){

  if(filter.samp==T){  
    missdata<-NULL
    for(i in 1:nrow(structure.file))
        missdata<-c(missdata,length(grep(-9,structure.file[i,]))/ncol(structure.file))
      
    hist(missdata, breaks = 200, xlab = "percentage of missing loci in a sample")
    
    if(is.null(cutoff)){
      cutoff<-readline("Cutoff to filter out the samples: ")
    }
    
    structure.file<-structure.file[-which(missdata>=cutoff),]

    } else if (filter.loci==T){  
    
      missdata<-NULL
    
      for(i in 1:ncol(structure.file))
          missdata<-c(missdata,length(grep(-9,structure.file[,i]))/nrow(structure.file))
    
        hist(missdata, breaks = 200, xlab = "percentage of missing samples in a locus")
        
        if(is.null(cutoff)){
          cutoff<-readline("Cutoff to filter out the loci: ")
        } 
    
    structure.file<-structure.file[-which(missdata>=cutoff),]
  }
  return(structure.file)
}

