#######################################################
# Script and functions to:
# * Read in climate data from GCMs daily projections
# * Extract values for coordinates of wine growing regions
# * Assess bias comparing GCMs vs. BEST - correct bias in GCMs
#
# Started 28th November 2016
# Ignacio Morales-Castilla & Elizabeth M. Wolkovich
#######################################################

## remove objects (activate if needed)
#rm(list=ls())
options(stringsAsFactors = FALSE)

## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','doParallel')
lapply(packs.to.extract,require, character.only=T)


## load functions
is.leapyear=function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}


#' function to extract data from GCMs 
#' (download from http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html)
#'
#' @param GCM.folder - directory with GCM simulations daily max & min temp data
#' @param GCMi - simulation of interest (from 1 to 30)
#' @param coords - data.frame with 3 columns: Site | Long | Lat
#' @param output.folder - where results are to be stored

#' @return list(climate.data.storing.min,climate.data.storing.max) - list 
#' with min and max temperatures for each selected location


Extract.GCM.fx<-function(GCM.folder,GCMi,coords,output.folder){
  
  k<-GCMi
  dfiles<-dir(GCM.folder)[1:50]
  names1<-c(rep("_g16.00",9),rep("_g16.0",21))
  names2<-seq(1,30,1)
  namesi<-paste(names1,names2,sep="")
  #namesi<-namesi[-c(16:20)]

  ## setting files name
  nc.name.min<-dfiles[which(grepl(namesi[k],dfiles) & grepl("TREFHTMN",dfiles))]
  nc.name.max<-dfiles[which(grepl(namesi[k],dfiles) & grepl("TREFHTMX",dfiles))]

  gcm.tempmn <- nc_open(paste(GCM.folder,nc.name.min,sep=""))
  gcm.tempmx <- nc_open(paste(GCM.folder,nc.name.max,sep=""))
  var.GCM <- names(gcm.tempmn$var)
  var.GCM.max <- names(gcm.tempmx$var)

  # open GCM as brick and rotate to -180 - 180 longitude
  daily.temp.GCM.min <- brick(paste(GCM.folder,nc.name.min,sep=""), varname=var.GCM[2], layer="time")
  daily.temp.GCM.max <- brick(paste(GCM.folder,nc.name.max,sep=""), varname=var.GCM.max[2], layer="time")

  # dates in the raster layer of our interest
  dates.in.brick<-as.Date(names(daily.temp.GCM.min),format="X%Y.%m.%d")
  dates.of.interest<-which(format.Date(dates.in.brick,"%Y")>=1950 &
                             format.Date(dates.in.brick,"%Y")<=2000)
  days.per.year<-seq(1,length(dates.of.interest),365)

  climate.data.storing.min<-list()
  climate.data.storing.max<-list()
  
  for(i in 1:2){
    print(i)

    GCM.day.imin<-daily.temp.GCM.min[[dates.of.interest[days.per.year[i]:(days.per.year[i+1]-1)]]]
    GCM.day.imax<-daily.temp.GCM.max[[dates.of.interest[days.per.year[i]:(days.per.year[i+1]-1)]]]

    climate.data.storing.min[[i]]<-t(extract(GCM.day.imin[[1:10]],coords[,2:3])-273)
    climate.data.storing.max[[i]]<-t(extract(GCM.day.imax[[1:10]],coords[,2:3])-273)

  }
  
  return(list(climate.data.storing.min,climate.data.storing.max))

}



