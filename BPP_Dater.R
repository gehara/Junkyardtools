##will date given ratograms from treePL (here using a bootstrap dist of trees), the mcmc output from BPP 4.0, the taxa to drop from the tree to get a localized rate, and generation time (gen)
library(fitdistrplus)
library(ape)
dater<-function (trees, drop, mcmc, genL, genU)
{   lapply(c(1:length(trees)), function(x) drop.tip(trees[[x]], drop))->t2
 
  mcmc2<- gsub("([.-])|[[:punct:]]", "\\1", mcmc)
 tau<- as.numeric(gsub(",","",mcmc2))
 
 ###estimates a prior for substition rate from trees 
 lapply(c(1:length(t2)), function(x) t2[[x]]$edge.length)->All_edge
 unlist(All_edge)->edge
 edge/1e6->rate
 mean(rate)->rateMean
 sd(rate)->rateSD
 
 subra = (  rateMean + sqrt(  rateMean^2 + 4*rateSD^2 ) ) / ( 2 * rateSD^2 )
 subsh = 1 +  rateMean * subra
 substRatesPr<-rgamma(length(tau),  subra, shape=subsh)
 
 
 #fitdist(rate, "gamma")->gamrate

 #substRatesPr<-rgamma(length(tau),  gamrate$estimate[1], gamrate$estimate[2])
 
 
###estimates a prior for generation time from data
 midpoint<-((genL-genU)/2)+genU
Gensd<- (genU-midpoint)/2
ragen = (  midpoint + sqrt(  midpoint^2 + 4*Gensd^2 ) ) / ( 2 * Gensd^2 )
shgen = 1 +  midpoint * ragen
GenPr<-rgamma(length(tau), ragen, shape=shgen)
 
substRatesPr/GenPr->substRateGen
  
tau/substRateGen->date

  result<-list()
  unlist(date)->date


  return(date)
}


  #example Langaha, first read all of your r8s trees from treePL in.
dropsy<-function(trees, taxa)
{
trees[[1]]$tip.label->tips

grep(taxa, tips)->d
trees[[1]]$tip.label[-c(d)]->drops
return(drops)}

dropsy(trees, "Langaha")->dropLang

read.table("Langaha.mcmc.txt"   )->Langaha

dater(trees, dropLang, Langaha$V2, 2,3)->date_Langaha


mg
use

val<-c(13545.55,17299, 39783, 11961, 33027,7000, 54819, 18527, 12623, 20819, 355, 14557, 7569, 30321,33932, 2851, 9049
       )