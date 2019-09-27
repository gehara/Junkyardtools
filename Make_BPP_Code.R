
##This will add a and b to your sequences for your IMAP file

add_post<-function(data)
  
{paste(data,"a",sep="")->one
  paste(data,"b",sep="")->two
  c(one,two)->three
  as.matrix(three)->four
  return(four)
}


#read in your table, with the first column (species names) 
lapply(c(1:length(tab[,1])), function(x) add_post(tab$Species[[x]]))->stuff

write.table(unlist(stuff), "names.txt")


###This part will convert the names for BPP
#name converter
bpp_namer<-function(genes)
{row.names(genes)<-paste("^",row.names(genes), sep="")
return(genes)}

##get your genes
list.files()->list
lapply(c(1:length(list)), function(x) read.dna(list[x]))->genes

#loop it

lapply(c(1:length(genes)), function(x) bpp_names(genes[[x]]))->newgenes

#write gene function
write.genes<-function(genes)
  for (i in 1:length(genes)) 
  {gene <- genes[[i]]
  name <- paste("gene", i, sep="_")
  write.phy(gene, file=name, interleave=F,strict=F)}

#change your directory
setwd("~/Downloads/Pan_phy_T60S5M40/rat_genes")

write.genes(newgenes)