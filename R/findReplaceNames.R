#' @title findReplaceNames
#' @description  This script will find and replace expressions in sequence names of fasta files.
#' @author Edward Myers and Marcelo Gehara
#' @param path.to.fasta character string. path to the directory of fasta files. 
#' @param input.id.name character string. Text indicating a wildcard to find the loci to process. 
#' @param expr.to.replace charactees to exclude from samples names
#' @param replacement  
#' @return processed sequence alignmnents to the to the output directory. 
#' @details Ex. if the names starts as DBS_355_uce-8046a, it will end up as DBS_355a (note the regex in 'gsub')
#' @export
findReplaceNames<-function(path.to.fasta,
                       input.id.name,
                       expr.to.replace = "(_uce-\\d+)",
                       replacement = "")
  {

  x<-list.files(pattern=input.id.name)
  for(i in 1:length(x)){
    fas <- read.dna(x[i],format="fasta")
    nam <- row.names(fas)
    gsub(expr.to.replace, replacement, nam, perl = TRUE) -> nam
    rownames(fas) <- nam
    write.dna(fas,x[i],format="fasta")
  }
  
  }
