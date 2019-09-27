install.packages("introgress")
library(introgress)
data("AdmixDataSim1")
help("AdmixDataSim1")
data("LociDataSim1")

#### get individuals assignment
setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/")
samp <- read.table("rat_imap.txt")

samp <- samp[grep("b",samp[,1]),]
samp[,1] <- gsub("b","",samp[,1])
samp12 <- data.frame(samp[which(samp[,2] %in% c("one","two")),1])
samp23 <- data.frame(samp[which(samp[,2] %in% c("three","two")),1])
samp34 <- data.frame(samp[which(samp[,2] %in% c("three","four")),1])

samp1 <- data.frame(samp[which(samp[,2] %in% c("one")),1])
samp4 <- data.frame(samp[which(samp[,2] %in% c("four")),1])
samp3 <- data.frame(samp[which(samp[,2] %in% c("three")),1])
samp2 <- data.frame(samp[which(samp[,2] %in% c("two")),1])


### calculate FSTs

##### pops 3 and 4
getwd()
library(PopGenome)
data <- readData("./fasta", big.data = T, FAST = T, parallized=F, include.unknown = T)
data <- diversity.stats(data)
fastas <- rownames(get.sum.data(data))

x <- read.csv("P_obsoletus_threshold_60_locale.csv")
pops <- list(as.character(samp3[,1]),as.character(samp4[,1]))
pops[[3]]<-as.character(x$INDV[which(!(as.character(x$INDV) %in% unlist(pops)))])
pops[[1]]<-c(paste(pops[[1]],"a", sep=""),paste(pops[[1]],"b", sep=""))
pops[[2]]<-c(paste(pops[[2]],"a", sep=""),paste(pops[[2]],"b", sep=""))
pops[[3]]<-c(paste(pops[[3]],"a", sep=""),paste(pops[[3]],"b", sep=""))

data <- set.populations(data,pops)
data <- F_ST.stats(data,FAST=F)

fsts<-get.F_ST(data, pairwise = F)

fsts[which(fsts[,1]>0.8),]

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/fasta")
y<-fastas[which(fsts[,1]>0.8)]
dir.create("../highfsts_34")
for(i in 1:length(y)){
  system(paste("cp",y[i],"../highfsts_34"))
}         

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/highfsts_34")
fasta2struc(samples2include = samp34, write.output = T, name = "34")

###################
############# pops 1 and 2

pops <- list(as.character(samp1[,1]),as.character(samp2[,1]))
pops[[3]]<-as.character(x$INDV[which(!(as.character(x$INDV) %in% unlist(pops)))])
pops[[1]]<-c(paste(pops[[1]],"a", sep=""),paste(pops[[1]],"b", sep=""))
pops[[2]]<-c(paste(pops[[2]],"a", sep=""),paste(pops[[2]],"b", sep=""))
pops[[3]]<-c(paste(pops[[3]],"a", sep=""),paste(pops[[3]],"b", sep=""))

data <- set.populations(data,pops)
data <- F_ST.stats(data,FAST=F)

fsts<-get.F_ST(data, pairwise = F)

fsts[which(fsts[,6]>.95),]

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/fasta")
y<-fastas[which(fsts[,6]>0.95)]
getwd()
dir.create("../highfsts_12")
for(i in 1:length(y)){
  system(paste("cp",y[i],"../highfsts_12"))
}       

setwd("../highfsts_12")
fasta2struc(samples2include = samp12, write.output = T, name="12")


###############33
################## pops 2 and 3


pops <- list(as.character(samp2[,1]),as.character(samp3[,1]))
pops[[3]]<-as.character(x$INDV[which(!(as.character(x$INDV) %in% unlist(pops)))])
pops[[1]]<-c(paste(pops[[1]],"a", sep=""),paste(pops[[1]],"b", sep=""))
pops[[2]]<-c(paste(pops[[2]],"a", sep=""),paste(pops[[2]],"b", sep=""))
pops[[3]]<-c(paste(pops[[3]],"a", sep=""),paste(pops[[3]],"b", sep=""))

data <- set.populations(data,pops)
data <- F_ST.stats(data,FAST=F)

fsts<-get.F_ST(data, pairwise = F)

fsts[which(fsts[,6]>0.8),]

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/fasta")
y<-fastas[which(fsts[,6] > 0.8)]
dir.create("../highfsts_23")
for(i in 1:length(y)){
  system(paste("cp",y[i],"../highfsts_23"))
}

setwd("../highfsts_23")
fasta2struc(samples2include = samp23, write.output = T, name="23")

#####################################################
#####################################################
####################################################

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/")
x<-read.csv("P_obsoletus_threshold_60_locale.csv")

admix.data <- read.table("12out.str")
ncol(admix.data)
nrow(admix.data)

P1 <- admix.data[which(admix.data[,1] %in% as.character(x$INDV[which(x$pop1>0.90)])),]
P2 <- admix.data[which(admix.data[,1] %in% as.character(x$INDV[which(x$pop2>0.90)])),]
ncol(P1)
ncol(P2)
nrow(P1)
nrow(P2)

admix.data <- admix.data[which(admix.data[,1]
                                %in% as.character(x$INDV[c(which(x$pop1<0.95 & x$pop1>0.5),which(x$pop2<0.95 & x$pop2>0.5))])),]

ncol(admix.data)
nrow(admix.data)

admix.data <- t(admix.data)
P1<-t(P1)
P2<-t(P2)

locusinfo <- data.frame(cbind(rownames(admix.data)[2:nrow(admix.data)],rep("C",nrow(admix.data)-1)))
colnames(locusinfo)<-c("locus","type")

nrow(locusinfo)
nrow(admix.data)
ncol(locusinfo)


write.table(admix.data,"admix.data12", col.names = F, row.names = F)
admix.data <- read.table("admix.data12", header = T, row.names = as.character(locusinfo[,1]))

write.table(P1,"P1", col.names = F, row.names = F)
P1 <- read.table("P1", header = T,row.names = as.character(locusinfo[,1]))

write.table(P2,"P2", col.names = F, row.names = F)
P2 <- read.table("P2", header = T, row.names = as.character(locusinfo[,1]))

format.this<-function(admix){
  new<-NULL
  i<-1
  while(i <= ncol(admix)-1){
    x<-admix[,i:(i+1)]
    
    x[,1]<-gsub(0,"A",x[,1])
    x[,2]<-gsub(0,"A",x[,2])
    
    x[,1]<-gsub(1,"D",x[,1])
    x[,2]<-gsub(1,"D",x[,2])
    
    x[,1]<-gsub(-9,NA,x[,1])
    x[,2]<-gsub(-9,NA,x[,2])
    
    
    y<-apply(x,1,paste,collapse="/")
    new<-cbind(new,y)
    i<-i+2
    print(i)
  }
  return(new)
}
admix.data<-format.this(admix.data)
P1<-format.this(P1)
P2<-format.this(P2)

admix.data <- admix.data[-47,]
P1 <- P1[-47,]
P2 <- P2[-47,]
locusinfo <- locusinfo[-47,]

row.names(admix.data)=="V106"

admix.data <- admix.data[-104,]
P1 <- P1[-104,]
P2 <- P2[-104,]
locusinfo <- locusinfo[-104,]


count.matrix <- prepare.data(admix.gen = admix.data, loci.data = locusinfo,
                             parental = P1, parental2 = P2, pop.id = F,
                             ind.id = F, fixed = F)

hi.index.sim <- est.h(introgress.data = count.matrix, loci.data = locusinfo,
                       fixed = F)

mk.image(introgress.data = count.matrix, loci.data = locusinfo,
         marker.order = NULL, hi.index = hi.index.sim, ylab.image = "individuals",
         xlab.h = "population 2 ancestry", pdf = TRUE, out.file = "plot12_final.pdf")
getwd()

gen.out <- genomic.clines(introgress.data = count.matrix, hi.index = hi.index.sim, loci.data = locusinfo,
                          sig.test = T, method="permutation")

dput(gen.out, "genout12.txt")

gen.out$Summary.data

clines.plot(cline.data=gen.out, rplots=3, cplots=3, pdf=T, out.file="clines12_highfst.pdf")


int.het <- calc.intersp.het(introgress.data=count.matrix)

introgress:::triangle.plot(hi.index=hi.index.sim, int.het=int.het, pdf=FALSE)



##################################################
##################################################
##################################################

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/")
x<-read.csv("P_obsoletus_threshold_60_locale.csv")

admix.data <- read.table("23out.str")
ncol(admix.data)
nrow(admix.data)

P2 <- admix.data[which(admix.data[,1] %in% as.character(x$INDV[which(x$pop2>0.90)])),]
P3 <- admix.data[which(admix.data[,1] %in% as.character(x$INDV[which(x$pop3>0.90)])),]
ncol(P2)
ncol(P3)
nrow(P2)
nrow(P3)

admix.data <- admix.data[which(admix.data[,1]
                               %in% as.character(x$INDV[c(which(x$pop3<0.90 & x$pop3>0.5),which(x$pop2<0.90 & x$pop2>0.5))])),]

ncol(admix.data)
nrow(admix.data)

admix.data <- t(admix.data)
P2<-t(P2)
P3<-t(P3)

locusinfo <- data.frame(cbind(rownames(admix.data)[2:nrow(admix.data)],rep("C",nrow(admix.data)-1)))
colnames(locusinfo)<-c("locus","type")

nrow(locusinfo)
nrow(admix.data)
ncol(locusinfo)


write.table(admix.data,"admix.data23", col.names = F, row.names = F)
admix.data <- read.table("admix.data23", header = T, row.names = as.character(locusinfo[,1]))

write.table(P2,"P2", col.names = F, row.names = F)
P2 <- read.table("P2", header = T,row.names = as.character(locusinfo[,1]))

write.table(P3,"P3", col.names = F, row.names = F)
P3 <- read.table("P3", header = T, row.names = as.character(locusinfo[,1]))

format.this<-function(admix){
  new<-NULL
  i<-1
  while(i <= ncol(admix)-1){
    x<-admix[,i:(i+1)]
    
    x[,1]<-gsub(0,"A",x[,1])
    x[,2]<-gsub(0,"A",x[,2])
    
    x[,1]<-gsub(1,"D",x[,1])
    x[,2]<-gsub(1,"D",x[,2])
    
    x[,1]<-gsub(-9,NA,x[,1])
    x[,2]<-gsub(-9,NA,x[,2])
    
    
    y<-apply(x,1,paste,collapse="/")
    new<-cbind(new,y)
    i<-i+2
    print(i)
  }
  return(new)
}
admix.data<-format.this(admix.data)
P2<-format.this(P2)
P3<-format.this(P3)

count.matrix <- prepare.data(admix.gen = admix.data, loci.data = locusinfo,
                             parental = P2, parental2 = P3, pop.id = F,
                             ind.id = F, fixed = F)

hi.index.sim <- est.h(introgress.data = count.matrix, loci.data = locusinfo,
                      fixed = F)

mk.image(introgress.data = count.matrix, loci.data = locusinfo,
         marker.order = NULL, hi.index = hi.index.sim, ylab.image = "individuals",
         xlab.h = "population 2 ancestry", pdf = TRUE, out.file = "plot23_final.pdf")
getwd()

gen.out <- genomic.clines(introgress.data = count.matrix, hi.index = hi.index.sim, loci.data = locusinfo,
                          sig.test = T, method="permutation")

dput(gen.out, "genout23.txt")

gen.out$Summary.data

clines.plot(cline.data=gen.out, rplots=3, cplots=3, pdf=T, out.file="clines23_highfst.pdf")


int.het<-calc.intersp.het(introgress.data=count.matrix)

introgress:::triangle.plot(hi.index=hi.index.sim, int.het=int.het, pdf=FALSE)



##################################################
##################################################
##################################################

setwd("/home/marcelo/Dropbox/AMNH/Pantherophis/")
x<-read.csv("P_obsoletus_threshold_60_locale.csv")

admix.data <- read.table("34out.str")
ncol(admix.data)
nrow(admix.data)

P3 <- admix.data[which(admix.data[,1] %in% as.character(x$INDV[which(x$pop3>0.90)])),]
P4 <- admix.data[which(admix.data[,1] %in% as.character(x$INDV[which(x$pop4>0.90)])),]
ncol(P3)
ncol(P4)
nrow(P3)
nrow(P4)

admix.data <- admix.data[which(admix.data[,1]
                               %in% as.character(x$INDV[c(which(x$pop3<0.90 & x$pop3>0.5),which(x$pop4<0.90 & x$pop4>0.5))])),]

ncol(admix.data)
nrow(admix.data)

admix.data <- t(admix.data)
P3<-t(P3)
P4<-t(P4)

locusinfo <- data.frame(cbind(rownames(admix.data)[2:nrow(admix.data)],rep("C",nrow(admix.data)-1)))
colnames(locusinfo)<-c("locus","type")

nrow(locusinfo)
nrow(admix.data)
ncol(locusinfo)


write.table(admix.data,"admix.data34", col.names = F, row.names = F)
admix.data <- read.table("admix.data34", header = T, row.names = as.character(locusinfo[,1]))

write.table(P3,"P3", col.names = F, row.names = F)
P3 <- read.table("P3", header = T,row.names = as.character(locusinfo[,1]))

write.table(P4,"P4", col.names = F, row.names = F)
P4 <- read.table("P4", header = T, row.names = as.character(locusinfo[,1]))

format.this<-function(admix){
  new<-NULL
  i<-1
  while(i <= ncol(admix)-1){
    x<-admix[,i:(i+1)]
    
    x[,1]<-gsub(0,"A",x[,1])
    x[,2]<-gsub(0,"A",x[,2])
    
    x[,1]<-gsub(1,"D",x[,1])
    x[,2]<-gsub(1,"D",x[,2])
    
    x[,1]<-gsub(-9,NA,x[,1])
    x[,2]<-gsub(-9,NA,x[,2])
    
    
    y<-apply(x,1,paste,collapse="/")
    new<-cbind(new,y)
    i<-i+2
    print(i)
  }
  return(new)
}
admix.data<-format.this(admix.data)
P3<-format.this(P3)
P4<-format.this(P4)

count.matrix <- prepare.data(admix.gen = admix.data, loci.data = locusinfo,
                             parental = P3, parental2 = P4, pop.id = F,
                             ind.id = F, fixed = F)

hi.index.sim <- est.h(introgress.data = count.matrix, loci.data = locusinfo,
                      fixed = F)

mk.image(introgress.data = count.matrix, loci.data = locusinfo,
         marker.order = NULL, hi.index = hi.index.sim, ylab.image = "individuals",
         xlab.h = "population 2 ancestry", pdf = TRUE, out.file = "plot34_final.pdf")
getwd()

gen.out <- genomic.clines(introgress.data = count.matrix, hi.index = hi.index.sim, loci.data = locusinfo,
                          sig.test = T, method="permutation")

dput(gen.out, "genout34.txt")

gen.out$Summary.data

clines.plot(cline.data=gen.out, rplots=3, cplots=3, pdf=T, out.file="clines34_highfst.pdf")


int.het<-calc.intersp.het(introgress.data=count.matrix)

introgress:::triangle.plot(hi.index=hi.index.sim, int.het=int.het, pdf=FALSE)










