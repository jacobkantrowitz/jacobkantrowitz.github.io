# Author: Joe Perez-Rogers
# Date: 2014-10-08
# Purpose: Clean up a directory by sorting files into file-type folders

cleanDir <- function(dir){
  f <- list.files(dir)
  extensions <- as.vector(sapply(f,function(x){strsplit(x,"[.]")[[1]][length(strsplit(x,"[.]")[[1]])]}))
  dirs <- names(table(extensions))
  dirs <- dirs[!c(dirs%in%c("Rmd","html"))]
  sapply(dirs,function(x){system(paste("mkdir ",dir,"/",x,sep=""))})
  for(fname in f){
    ext <- strsplit(fname,"[.]")[[1]][length(strsplit(fname,"[.]")[[1]])]
    if(!c(ext%in%c("Rmd","html"))){
      system(paste("mv ",dir,"/",fname," ",dir,"/",ext,"/",fname,sep=""))
    } else {
      cat(paste(fname,"was not moved because of its file extension.\n",sep=" "))
    }
  }
}