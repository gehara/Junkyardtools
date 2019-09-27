#' @title excludeFromFasta
#' @description this is to remove samples from final fasta alignments
#' @param path character string. Path to the fasta alignments
#' @param samples2exclude one column data frame with the name of the samples to exclude from fasta alignments.
#' @param out.prefix.name character string. A character string to use as prefix of new fasta alignments. default is "out" 
#' @author Marcelo Gehara
#' @export
excludeFromFasta<-function(path, samples2exclude, out.prefix.name="out") {
  
    setwd(path)
    out.dir<-dir.exists("output")
    if(out.dir==F) dir.create("output")
    #find your alignments using some unique feature of the names
    fasta<-list.files(pattern=".fa")
    fasta<-fasta[grep(".fa",fasta,fixed=T)]

    #run this bullshit.
    mess<-NULL
    for (i in 1:length(fasta)){
      x<-read.dna(fasta[i], format="fasta")
      y<-which(rownames(x) %in% samples2exclude[,1])
      if(length(y)==0){
        setwd("./output")
        write.dna(x,paste(out.prefix.name, "_", fasta[i], sep=""), forma="fasta", colw = 1000)
        setwd("../")
        next
      }else{
      x<-x[-y,]
      if(length(x)==0){
        mess<-c(mess,paste("zero samples in locus",fasta[i],", locus excluded!"))
                next
      }
      setwd("./output")
      write.dna(x,paste(out.prefix.name, "_", fasta[i], sep=""), forma="fasta", colw = 1000)
      setwd("../")
      cat(paste("writing locus > ",fasta[i],"|",length(fasta)-i,"loci to go"),sep="\n")
      }
      
    }
    warning(cat(mess,sep="\n"))
}
