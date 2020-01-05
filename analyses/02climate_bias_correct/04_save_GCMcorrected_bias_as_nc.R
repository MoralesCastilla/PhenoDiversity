#'######################################################
#' Script to save bias-corrected GCMs in 03_function_correcting_GCM_bias.R as .nc:
#'
#' Started 20th January 2017
#' Ignacio Morales-Castilla & Elizabeth M. Wolkovich
#'######################################################

## remove objects
rm(list=ls())
options(stringsAsFactors = FALSE)

## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','doParallel','abind')
lapply(packs.to.extract,require, character.only=T)


## auxiliary functions and vectors
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# vectors to divide up the tasks across years and day of the year
seqs<-list(seq(1,99,1),seq(100,198,1),seq(199,297,1),seq(298,365,1))
list.years<-list(seq(1950,2005,1),seq(2006,2080,1),seq(2081,2100,1))

#' Function to save as .nc files the results produced in the script
#' '03_function_correcting_GCM_bias.R'
#' @param GCMnum - number of GCM simulation to be corrected 
#' @param t - each year period in vector list.years: value 1 = 1950-2005; 
#'          value 2 = 2006-2080; value 3 = 2081-2100 
#' @param target.folder - local folder where results from '03_function_correcting_GCM_bias.R'
#'          are stored          
#' @return list with two rasterStacks: GCM.bias.corr.min and GCM.bias.corr.max
#'
saveGCMcorr.as.nc <- function(j,GCMnum,target.folder){  
  
  
  GCM.i<-paste("GCM.",GCMnum,'.',sep="")
  
  
  ## declare lists to store results
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
  GCM.bias.corr.max<-stack(unlist(storing.list.maxs))
 return(list(GCM.bias.corr.min,GCM.bias.corr.max))  
}
 

#### example of application ####

target.folder <- "~/your.local.dir.here"

a<-saveGCMcorr.as.nc(t,GCMnum)
 
GCM.bias.corr.min<-a[[1]]
GCM.bias.corr.max<-a[[2]]

writeRaster(GCM.bias.corr.min,paste(target.folder,GCM.i,"Tmin",
                                    min(list.years[[t]]),"-",
                                    max(list.years[[t]]),".nc",sep=""),
            format="CDF",overwrite=TRUE)
rm(GCM.bias.corr.min,storing.list.mins)

writeRaster(GCM.bias.corr.max,paste(target.folder,GCM.i,"Tmax",
                                    min(list.years[[t]]),"-",
                                    max(list.years[[t]]),".nc",sep=""),
            format="CDF",overwrite=TRUE)
rm(GCM.bias.corr.max,storing.list.maxs)

