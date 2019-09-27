#' fasta2VCF
#' @description  Converts a colection of phased fasta alignments into a VCF file  
#' @author Gehara, M 
#' @param path path to the directory where all fasta alignments are.
#' @param samples2include one column data frame with the name of the samples to include in the first column  
#' @return out.vcf file to the path directory
#' @export
fasta2VCF <- function(path, samples2include) {
  
  library(ape)
  setwd(path)
  samples<-as.character(sort(samples2include[,1]))
  
  x<-list.files(pattern=".fa")
  x<-x[grep(".fa",x,fixed=T)]
  
  VCF<-matrix(ncol=9+length(samples),nrow=0)
  colnames(VCF)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",samples)
  write.table(VCF,"out.vcf",col.names=T,row.names=F,quote=F,sep="\t")
  
  for(u in 1:length(x)){
    
    chrom<-x[u]
    
    seq <- read.dna(x[u],format = "fasta")
    seq <- seq[match(rownames(seq),sort(rownames(seq))),]
      
    fas <- as.character(seq)
    bin <- NULL
    pos <- NULL
    for (i in 1:ncol(fas)) {
      a <- length(grep("a", fas[, i]))
      c <- length(grep("c", fas[, i]))
      g <- length(grep("g", fas[, i]))
      t <- length(grep("t", fas[, i]))
      r <- length(grep("r", fas[, i]))
      y <- length(grep("y", fas[, i]))
      m <- length(grep("m", fas[, i]))
      k <- length(grep("k", fas[, i]))
      s <- length(grep("s", fas[, i]))
      w <- length(grep("w", fas[, i]))
      h <- length(grep("h", fas[, i]))
      b <- length(grep("b", fas[, i]))
      v <- length(grep("v", fas[, i]))
      d <- length(grep("d", fas[, i]))
      n <- length(grep("n", fas[, i]))
      if (n > 0) {
        next
      }
      if(!nrow(fas) %in% c(a, c, g, t, r, y, m, k, s, w, h, 
                           b, v, d)) {
        bin <- cbind(bin, fas[, i])
        pos <- c(pos, i)
      }
    }
    
    if(is.null(bin)){next}
    
    for(z in 1:ncol(bin)){
      unico<-unique(bin[,z])
      
      if(length(unico)>2){
        warning(cat(paste("locus",x[u],"has more than two mutations in base position",pos[z]),
                    paste("Position excluded!"),sep="\n"))
        next
        }
      
      for(j in 1:length(unico)){
        bin[,z] <- gsub(unico[j],j-1,bin[,z])
      }
      
      bin2 <- bin[row.names(bin),z]
      bin2 <- data.frame(bin2)
      samp <- rownames(bin2)
      
      snps<-NULL
      for(w in 1:length(samples)){
        ss<-grep(samples[w],samp)
        if(length(ss)==0){
          snps<-c(snps,"./.")
        } else {
          snps<-c(snps,paste(bin2[ss,],collapse="|"))
        }
      }
      write.table(t(data.frame(c(chrom,pos[z],".",unico,20,"PASS",".","GT",snps))),"out.vcf",
                  append=T,quote=F,col.names = F,row.names = F,sep="\t")
      
    }
    print(u)
  }
  
  #return(VCF)
}
