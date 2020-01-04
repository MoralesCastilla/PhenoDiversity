#'######################################################
#' Script to run mean-replacement bias-correction for GCMs:
#' 
#' Read in climate data from GCMs daily projections
#' Get climatology (mean for each day from 1950-1990 for each site)
#' Bias correction of GCMs with respect to BEST data
#' 
#' Started 20th January 2017
#' Ignacio Morales-Castilla, Benjamin I. Cook & Elizabeth M. Wolkovich
#'######################################################

## remove objects
rm(list=ls())
options(stringsAsFactors = FALSE)

## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','doParallel','abind')
lapply(packs.to.extract,require, character.only=T)

#### Functions ####

#' Bias.correction.GCM.mean.rep - this function corrects daily temperature
#' data in GCMs by replacing the climatology of GCM by that of BEST data.
#' @param day - escalar indicating day of the year (1 to 365)
#' @param GCMi - escalar indicating simulation number of GCM (1 to 30)
#' @param GCM.folder - local folder containing GCM raw data (from
#'        http://www.cesm.ucar.edu/projects/community-projects/LENS/)    
#' @param output.folder - local folder to store results
#' @param corrGCM.folder - local folder containing GCM daily climatologies (60 files*)
#' 
#' @return GCM.corr - raster stack with minimum and maximum temperatures corrected 
#' 
#' Note that: 
#' 1) the below funtion needs the following files to be in the working directory:
#' -"BEST.mean.clim_min.nc" - multi-layer raster file with 365 layers corresponding
#' to each day of the year, with values for minimum daily temperatures averaged 
#' across the 1950-1990 period, from BEST climate data 
#' -"BEST.mean.clim_max.nc" - multi-layer raster file with 365 layers corresponding
#' to each day of the year, with values for maximum daily temperatures averaged 
#' across the 1950-1990 period, from BEST climate data
#' 
#' 2) *the below funtion needs the following files to be in the corrGCM.folder:
#' -"GCM.climatology.tmin_XX.nc" - multi-layer raster file with 365 layers corresponding
#' to each day of the year, with values for minimum daily temperatures averaged 
#' across the 1950-1990 period, from GCM climate data 
#' -"GCM.climatology.tmax_XX.nc" - multi-layer raster file with 365 layers corresponding
#' to each day of the year, with values for maximum daily temperatures averaged 
#' across the 1950-1990 period, from GCM climate data
#' XX in file names should be replaced by values 1-30 corresponding to each
#' GCM simulation, which yields a total of 60 files.
#' 

Bias.correction.GCM.mean.rep<-function(day,GCMi,GCM.folder,
                                       output.folder,corrGCM.folder){

k<-as.numeric(GCMi)
dfiles<-dir(GCM.folder)#[1:40]
names1<-c(rep("_g16.00",9),rep("_g16.0",21))
names2<-seq(1,30,1)
names3<-c("20051231","20801231","21001231")
namesi<-paste(names1,names2,sep="")
list.years<-list(seq(1950,2005,1),seq(2006,2080,1),seq(2081,2100,1))
dcorr.files<-dir(corrGCM.folder)
dcorr.names<-paste("_",k,".nc",sep="")

## reading in BEST climatologies from nc
BEST.mean.clim.min <- brick("BEST.mean.clim_min.nc")
BEST.mean.clim.max <- brick("BEST.mean.clim_max.nc")


## setting files name
nc.name.min <- dfiles[which(grepl(namesi[k],dfiles) & 
                            grepl("TREFHTMN",dfiles) & grepl(names3[t],dfiles))]
nc.name.max <- dfiles[which(grepl(namesi[k],dfiles) & 
                            grepl("TREFHTMX",dfiles) & grepl(names3[t],dfiles))]

daily.temp.GCM.min <- brick(paste(GCM.folder,nc.name.min,sep=""), 
                            varname="TREFHTMN", layer="time")
daily.temp.GCM.max <- brick(paste(GCM.folder,nc.name.max,sep=""), 
                            varname="TREFHTMX", layer="time")

## setting files name for climatology
nc.name.min.corr<-dcorr.files[which(grepl(dcorr.names,dcorr.files) & 
                                      grepl("tmin",dcorr.files))]
nc.name.max.corr<-dcorr.files[which(grepl(dcorr.names,dcorr.files) & 
                                      grepl("tmax",dcorr.files))]

# open climatology as brick #
mean.each.min <- brick(paste(corrGCM.folder,nc.name.min.corr,sep=""))
mean.each.max <- brick(paste(corrGCM.folder,nc.name.max.corr,sep=""))

## names of each layer
name.layer<-names(daily.temp.GCM.min)
name.layer2<-names(daily.temp.GCM.max)
date.layer<-as.Date(names(daily.temp.GCM.min),format="X%Y.%m.%d")
date.layer2<-as.Date(names(daily.temp.GCM.max),format="X%Y.%m.%d")
dates<-seq.Date(from=as.Date(paste(min(list.years[[t]]),"/1/1",sep="")),
                to=as.Date(paste(max(list.years[[t]]),"/12/31",sep="")),by=1)


## add layer if 1st jan not present and remove from last year
if(!dates[1]%in%date.layer){
  ## re-ordering dates
  first.jans<-format(date.layer[format.Date(date.layer,"%m")=="01" 
                                             & format.Date(date.layer,"%d")=="01"][1:4],"X%Y.%m.%d")
  first.jan.min<-calc(daily.temp.GCM.min[[first.jans]],mean) #first of jan as mean two following years
  
  if(length(date.layer[!date.layer%in%dates])>0){
  daily.temp.GCM.min<-dropLayer(daily.temp.GCM.min,nlayers(daily.temp.GCM.min))
  }
  
  daily.temp.GCM.min<-addLayer(daily.temp.GCM.min,first.jan.min)
  names(daily.temp.GCM.min[[nlayers(daily.temp.GCM.min)]])<-dates[1]

}

if(!dates[1]%in%date.layer2){
  ## re-ordering dates
  first.jans<-format(date.layer2[format.Date(date.layer2,"%m")=="01" 
                                & format.Date(date.layer2,"%d")=="01"][1:4],"X%Y.%m.%d")
  if(length(date.layer2[!date.layer2%in%dates])>0){
    daily.temp.GCM.max<-dropLayer(daily.temp.GCM.max,nlayers(daily.temp.GCM.max))
  }
  
  first.jan.max<-calc(daily.temp.GCM.max[[first.jans]],mean) #first of jan as mean two following years
  daily.temp.GCM.max<-addLayer(daily.temp.GCM.max,first.jan.max)
  names(daily.temp.GCM.max[[nlayers(daily.temp.GCM.max)]])<-dates[1]
}


## drop layers (if needed only)
not.in.dates.min<-which(!names(daily.temp.GCM.min)%in%format(dates,"X%Y.%m.%d"))
not.in.dates.max<-which(!names(daily.temp.GCM.max)%in%format(dates,"X%Y.%m.%d"))

if(length(not.in.dates.min)>0){
daily.temp.GCM.min<-dropLayer(daily.temp.GCM.min,not.in.dates.min)
}
if(length(not.in.dates.max)>0){
daily.temp.GCM.max<-dropLayer(daily.temp.GCM.max,not.in.dates.max)
}


## correct layer names for leap years 
## remove 29th feb
feb29s<-dates[format.Date(dates,"%m")=="02" & format.Date(dates,"%d")=="29"]
dates<-dates[-which(dates%in%feb29s)]

## missing layers (leap years have changed 12.31 by 2.29)
missing.layers.min<-dates[which(!format(dates,"X%Y.%m.%d")%in%names(daily.temp.GCM.min))]
missing.layers.max<-dates[which(!format(dates,"X%Y.%m.%d")%in%names(daily.temp.GCM.max))]


## ordering and naming
if(length(dates)==nlayers(daily.temp.GCM.min)){
daily.temp.GCM.min<-daily.temp.GCM.min[[order(names(daily.temp.GCM.min))]]
daily.temp.GCM.max<-daily.temp.GCM.max[[order(names(daily.temp.GCM.max))]]
names(daily.temp.GCM.min)<-dates
names(daily.temp.GCM.max)<-dates
}


## replacing means: take daily value, substract GCM-mean for that date, add BEST-mean for that date
## function to replace mean
    daily.index<-names(daily.temp.GCM.min)[seq(day,nlayers(daily.temp.GCM.min),365)]
    GCM.corr.min<-rotate(daily.temp.GCM.min[[daily.index]]-273)-rotate(mean.each.min[[day]]-273)+BEST.mean.clim.min[[day]]
    GCM.corr.max<-rotate(daily.temp.GCM.max[[daily.index]]-273)-rotate(mean.each.max[[day]]-273)+BEST.mean.clim.max[[day]]


## naming the layers
nyears<-list.years[[t]]
DOY<-day
names.layers<-list()
for (j in 1:length(DOY)){
DAY<-DOY[j]  

names.layers.year<-array(NA,dim=c(length(nyears),1))
for(h in 1:length(nyears)){
  year.i<-nyears[h]
names.layers.year[h,1]<-paste(year.i,DAY,sep=".")
    
}
names.layers[[j]]<-names.layers.year
}
names.layers<-unlist(names.layers)

names(GCM.corr.min)<-names.layers[1:nlayers(GCM.corr.min)]
names(GCM.corr.max)<-names.layers[1:nlayers(GCM.corr.max)]


GCM.corr<-stack(GCM.corr.min,GCM.corr.max)
return(GCM.corr)

}


