#data is in nexus format, folder= ".", make sure you have navigated to that folder or write the link ahead of time, type="absolute", "fraction" of PIS.
library(ape)
library(ips)

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


JYD<-function(folder, type)
{
  list.files(folder)->files
  lapply(c(1:length(files)), function(x) read.dna(files[[x]]))->dat
  lapply(c(1:length(dat)), function(x) as.alignment(dat[[x]]))->dat2
  lapply(c(1:length(dat2)), function(x) as.DNAbin(dat2[[x]]))->dat3
  lapply(c(1:length(dat3)), function(x) site_wise(dat3[[x]]))->dat3
  lapply(c(1:length(dat3)), function(x) pis(dat3[[x]], type))->dat3
  unlist(dat3)->dat3
  density(dat3)->dens
  plot(dens,ylab="density",xlab="pars inf sites", main="Density of Parsimony Informative Sites Up in Here", col="blue")
  polygon(dens, col="indianred4")
  files->loci
  dat3->pis
  data.frame(files, pis)->out
  return(out)}