Nia_Jax<-function(x)
{list.files(x)->list
  lapply(c(1:length(list)), function(x) read.dna(list[x], format="fasta"))->dat
  site_wise<-function(x)
  {x <- as.matrix(x)
  nms <- dimnames(x)[[1]]
  n <- dim(x)
  s <- n[2]
  n <- n[1]
  keep <- .C(GlobalDeletionDNA, x, n, s, rep(1L, s))[[4]]
  x <- x[, as.logical(keep)]
  s <- dim(x)[2]
  return(x)}
  
  lapply(c(1:length(dat)), function(x) site_wise(dat[[x]]))->genes
  
  
  convert_fasta<-function(genes)
    for (i in 1:length(genes)) 
    {gene <- genes[[i]]
    name <- paste("gene", i, sep="_")
    write.dna(gene,format="fasta", file=name, colsep = "", indent = "\t")}
  convert_fasta(genes)}