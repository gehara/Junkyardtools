#' fasta2struc
#' @description THIS SCRIPT PULLS OUT A RANDOM SNP FROM THE ALIGNMENTS TO MAKE A STRUCTURE FILE.
#' @author Gehara, M.
#' @param path path to the directory where all fasta alignments are.
#' @param samples2include one column data frame with the name of the samples to include in the first column.
#' @param write.output logical. If TRUE writes the structure file to working directory.
#' @return Structure-like genotype matrix.
#' @export
fasta2struc<-function(path, samples2include, write.output=TRUE) {
  
  library(ape)
  
  samples<-as.character(samples2include[,1])
  samples<-sort(c(samples,samples))
    
  x <- list.files(pattern=".fa")
  x <- x[grep(".fa", x, fixed=T)]

  dat <- rep(-9,(length(samples)*length(x)))
  structure.file <- matrix(dat,nrow=length(samples),ncol=length(x))
  rownames(structure.file) <- samples

  for(u in 1:length(x)){

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
  
    for(i in 1:ncol(bin)){
      unico <- unique(bin[,i])
      for(j in 1:length(unico)){
        bin[,i]<-gsub(unico[j],j-1,bin[,i])
        }
      }
  
      bin <- data.frame(bin)
      samp <- rownames(bin)

      zzz<-lapply(as.character(samples2include[,1]), function(y) grep(y, samp, fixed=T))
      len<-lapply(zzz,length)
      if(length(which(unlist(len)>=3))>0){
        warning("sample names match more than two alleles!")
      }      
      
      str<-NULL
      for(i in 1:length(zzz)){
        if(length(zzz[[i]])==0){
          str<-c(str,c(-9,-9))
        }
        str<-c(str,as.character(bin[c(zzz[[i]]),sample(1:ncol(bin),1)]))
      }
    
      structure.file[,u]<-str
  
      cat(paste("loci",">", u, x[u], "|", length(x)-u,"to go"), sep="\n")
      }

  cols <- NULL
    for(i in 1:ncol(structure.file)){
      if(length(grep(-9,structure.file[,i]))==nrow(structure.file)){
        cols<-c(cols,i)
      }
      }
  structure.file <- structure.file[,-cols]
  
  if(write.output==T)
  write.table(structure.file, "out.str",quote = F, row.names = T, col.names = F,sep=" ")

  return(structure.file)
}
