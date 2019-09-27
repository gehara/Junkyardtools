#' alignerTrimmer
#' @description This function aligns and trims all fasta sequences present in a directory
#' @author Marcelo Gehara and Edward Myers
#' @param path.to.fasta character string. path to the directory of fasta files. 
#' @param samples.cutoff integer. Exclude loci that have less samples than the sample.cutoff.
#' @param quantile.cutoff float. Exclude all sequences with length bellow the quantile 
#'                               of the distribution of lengths across all samples. This is to exclude samples with short sequences before trimming.
#' @param segsites.cutoff float. Exclude locus if the number segsites is larger than the segsites.cutoff of total length of the alignment.
#' @param remove.gaps logical. If TRUE remove gaps and trim alignment.
#' @param input.id.name character string. Text indicating a wildcard to find the loci to process. 
#' @param output.prefix character string. Text prefix to identify the output.
#' @param output.dir character string. Path to the output directory.
#' @return processed sequence alignmnents to the to the output directory. 
#' @export
alignerTrimmer<-function(path.to.fasta = getwd(), samples.cutoff=10, 
                         quantile.cutoff = 0.25, 
                         segsites.cutoff = 0.5, 
                         remove.gaps = T,
                         input.id.name = "uce-",
                         output.prefix = "ALTRIM_",
                         output.dir = "OUT")
  {
  
  #path<-getwd()
  
  if(dir.exists(output.dir)==F){
    dir.create(output.dir)
  }
  setwd(path.to.fasta)
  x <- list.files(pattern = input.id.name)
  x <- x[grep(input.id.name,x,fixed=T)]
  
  for(j in 1:length(x)){
    
    seq<-paste(path.to.fasta,"/",x[j],sep="")
    
    seq <- readDNAStringSet(seq,format="fasta")
    
    if(length(names(seq))<samples.cutoff) next
    
    len <- width(seq)
    misdata <- which(len<quantile(len, probs = quantile.cutoff))
    seq <- seq[-misdata,]
    
    if(length(names(seq))<samples.cutoff) next
    
    if(length(seq)==0) next
    
    seq<-AlignSeqs(seq, verbose = F)
    
    if(width(seq[1])==0) next
    
    seq<-DNAMultipleAlignment(seq)
    
    if(remove.gaps==T) seq<-maskGaps(seq, min.fraction=0.001, min.block.width=1)
    
    if(is.null(ncol(as.matrix(seq)))) next
    
    if(ncol(as.matrix(seq))==0) next
    
    seq<-as.DNAbin(seq)
    
    if(length(seg.sites(seq))>=(ncol(seq)*0.5)) next
    #checkAlignment(seq)
    #readline(prompt="Press [enter] to continue")
    
    print(paste(j,"pi",nuc.div(seq),"| seg sites",length(seg.sites(seq)),"| length",ncol(seq)))
    setwd(paste(path.to.fasta,"/",output.dir,sep=""))
    write.dna(seq,paste(output.prefix,x[j],sep=""),format = "fasta",colw = 1000)
    setwd(path.to.fasta)
  }
  #setwd(path)
}

