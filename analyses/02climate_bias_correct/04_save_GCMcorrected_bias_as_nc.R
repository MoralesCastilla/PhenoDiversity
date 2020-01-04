#######################################################
# Script to run mean replacement bias-correction globally in RC:
#' Read in climate data from GCMs daily projections
#' Get climatology (mean for each day from 1950-1990 for each site)
#' bias correction for coastal cells in BEST data
#'
#'
#'
#' Started 20th January 2017
# Ignacio Morales-Castilla, Benjamin Cook & Elizabeth M. Wolkovich
#######################################################

## remove objects
rm(list=ls())
options(stringsAsFactors = FALSE)

## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','doParallel','abind')
lapply(packs.to.extract,require, character.only=T)

target.folder="/n/wolkovich_lab/climatefuture/corrected_GCMs/"

seqs<-list(seq(1,99,1),seq(100,198,1),seq(199,297,1),seq(298,365,1))
list.years<-list(seq(1950,2005,1),seq(2006,2080,1),seq(2081,2100,1))
#setwd("~/MEGA/Work_Harvard_postdoc/vitis/data/R/")
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## loading
#for(t in 1:1){
#for(j in 2:4){


  
  t=2
  j=16
  GCM.i<-paste("GCM.",j,'.',sep="")
  storing.list.mins<-list()
  storing.list.maxs<-list()
  
  for(i in 1:4){
    print(i)
    temp.i<-paste(min(seqs[[i]]),max(seqs[[i]]),sep="-")
    years.i<-paste(min(list.years[[t]]),max(list.years[[t]]),sep="-")
    name.rda<-paste(target.folder,dir(target.folder)[which(grepl(GCM.i,dir(target.folder),fixed=T) 
                                                           & grepl(temp.i,dir(target.folder),fixed=T)
                                                           & grepl(years.i,dir(target.folder),fixed=T))],sep="")
    
    
    
    load(name.rda)
    namesgcm<-names(GCM.bias.corr)
    last.char<-substrRight(namesgcm,1)
    storing.list.mins[[i]]<-GCM.bias.corr[[which(last.char=="1")]]
    storing.list.maxs[[i]]<-GCM.bias.corr[[which(last.char=="2")]]
    
  }
  
  GCM.bias.corr.min<-stack(unlist(storing.list.mins))
  writeRaster(GCM.bias.corr.min,paste(target.folder,GCM.i,"Tmin",
                                      min(list.years[[t]]),"-",max(list.years[[t]]),".nc",sep=""),format="CDF",overwrite=TRUE)
  rm(GCM.bias.corr.min,storing.list.mins)
  
  GCM.bias.corr.max<-stack(unlist(storing.list.maxs))
  writeRaster(GCM.bias.corr.max,paste(target.folder,GCM.i,"Tmax",
                                      min(list.years[[t]]),"-",max(list.years[[t]]),".nc",sep=""),format="CDF",overwrite=TRUE)
  rm(GCM.bias.corr.max,storing.list.maxs)
  
  
  