#'######################################################
#' Script to run mean-replacement bias-correction for GCMs:
#' 
#' Read in climate data from GCMs daily projections
#' Get climatology (mean for each day from 1950-1990 for each site)
#' bias correction with respect to BEST data
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

## setting raster options to optimize
#rasterOptions(format="CDF",overwrite=TRUE,maxmemory = 1e09, chunksize=1e08,progress="text") 
#rasterTmpFile("clean_this_after_")
#rasterOptions()

## source function
#setwd("~/MEGA/Work_Harvard_postdoc/vitis/data/R/")
#source("Meansubs_GCMs.R")

## load avg BEST climatologies
#BEST.folder="~/MEGA/Work_Harvard_postdoc/vitis/data/Climate Data/BEST/"
#setwd(BEST.folder)
#BEST.mean.clim<-readRDS("BEST_avg_climatologies_interpolated.rds")
#BEST.mean.clim.min<-brick("BEST.mean.clim_min.nc")
#BEST.mean.clim.max<-brick("BEST.mean.clim_max.nc")


## load avg GCM climatologies 
#setwd(GCM.folder)
#corrGCM.folder<-"D:/Vitis_Data/climatefuture/GCMs_corrected/"
corrGCM.folder<-"/n/wolkovich_lab/climatefuture/GCM_climatology/"
#GCM.climatology<-GCM.climatology

## function to extract and compare data
GCM.folder="/n/wolkovich_lab/climatefuture/"
#output.folder="/n/wolkovich_lab/Lab/Nacho/"
output.folder="/n/wolkovich_lab/climatefuture/corrected_GCMs/"
#GCM.folder="D:/Vitis_Data/climatefuture/GCMs/"
#output.folder="D:/Vitis_Data/climatefuture/GCMs_corrected/"
#GCMi<-2


# open GCM as brick #
Bias.correction.GCM.mean.rep<-function(day,GCMi,GCM.folder,output.folder,corrGCM.folder){

#for(k in 1:20){
    
#GCMi<-1
k<-as.numeric(GCMi)
dfiles<-dir(GCM.folder)[1:40]
#dfiles.fut<-dir(GCM.folder)[1:40]
names1<-c(rep("_g16.00",9),rep("_g16.0",21))
#names1<-c(rep("_g16.00",4))
names2<-seq(1,30,1)
#names2<-seq(1,3,1)
names3<-c("20051231","20801231","21001231")
namesi<-paste(names1,names2,sep="")
#namesi<-namesi[-c(16:20)]
list.years<-list(seq(1950,2005,1),seq(2006,2080,1),seq(2081,2100,1))
dcorr.files<-dir(corrGCM.folder)
dcorr.names<-paste("_",k,".nc",sep="")

## reading in BEST climatologies from nc
BEST.mean.clim.min<-brick("BEST.mean.clim_min.nc")
BEST.mean.clim.max<-brick("BEST.mean.clim_max.nc")


#for(t in 1:1){ #t=1 
## setting files name
nc.name.min<-dfiles[which(grepl(namesi[k],dfiles) & grepl("TREFHTMN",dfiles) & grepl(names3[t],dfiles))]
nc.name.max<-dfiles[which(grepl(namesi[k],dfiles) & grepl("TREFHTMX",dfiles) & grepl(names3[t],dfiles))]

daily.temp.GCM.min <- brick(paste(GCM.folder,nc.name.min,sep=""), varname="TREFHTMN", layer="time")
daily.temp.GCM.max <- brick(paste(GCM.folder,nc.name.max,sep=""), varname="TREFHTMX", layer="time")

## setting files name for climatology
nc.name.min.corr<-dcorr.files[which(grepl(dcorr.names,dcorr.files) & grepl("tmin",dcorr.files))]
nc.name.max.corr<-dcorr.files[which(grepl(dcorr.names,dcorr.files) & grepl("tmax",dcorr.files))]

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


## unlisting and stacking layers
#GCM.corrected.min<-stack(unlist(GCM.corr.min))
#GCM.corrected.max<-stack(unlist(GCM.corr.max))
#GCM.corrected.min<-GCM.corr.min
#GCM.corrected.max<-GCM.corr.max
    

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

#writeRaster(GCM.corrected.min,paste(output.folder,"GCM.",k,".Tmin",
#            min(list.years[[t]]),"-",max(list.years[[t]]),".nc",sep=""),format="CDF",overwrite=TRUE)
#writeRaster(GCM.corrected.max,paste(output.folder,"GCM.",k,".Tmax",
#            min(list.years[[t]]),"-",max(list.years[[t]]),".nc",sep=""),format="CDF",overwrite=TRUE)


#}
GCM.corr<-stack(GCM.corr.min,GCM.corr.max)
return(GCM.corr)

}


# paralellize code
seqs<-list(seq(1,99,1),seq(100,198,1),seq(199,297,1),seq(298,365,1))
list.years<-list(seq(1950,2005,1),seq(2006,2080,1),seq(2081,2100,1))
cl <- makeCluster(24)
registerDoParallel(cl)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

for (j in c(1,16,17,18,19)){
  GCMi<-j
  #list.store.mins<-list()
  #list.store.maxs<-list()
  
  
  for (i in 1:4){
    print(i)
    t=1
    #k= #i
    
    Sys.time()
    GCM.bias.corr<-foreach(day = seqs[[i]], .packages=c("raster","ncdf4"),
                           .verbose=T,.errorhandling="pass",.combine=stack)  %dopar%  
      Bias.correction.GCM.mean.rep(day,GCMi,GCM.folder,output.folder,corrGCM.folder)
    
    Sys.time()
    #namesgcm<-names(GCM.bias.corr)
    #last.char<-substrRight(namesgcm,1)
    #list.store.mins[[i]]<-GCM.bias.corr[[which(last.char=="1")]]
    #list.store.maxs[[i]]<-GCM.bias.corr[[which(last.char=="2")]]
    save(GCM.bias.corr,file=paste(output.folder,"GCM.",GCMi,".Tempboth",
                                  min(list.years[[t]]),"-",max(list.years[[t]]),"_day",
                                  min(seqs[[i]]),"-",max(seqs[[i]]),".RData",sep=""))
    
  }
  
  #GCM.bias.corr.min<-stack(unlist(list.store.mins))
  #GCM.bias.corr.max<-stack(unlist(list.store.maxs))
  
  #save(GCM.bias.corr[[which(last.char=="1")]],file=paste(output.folder,"GCM.",GCMi,".Tempmin",
  #                                  min(list.years[[t]]),"-",max(list.years[[t]]),"_day",
  #                                  min(seqs[[i]]),"-",max(seqs[[t]]),".RData",sep=""))
  #save(GCM.bias.corr[[which(last.char=="2")]],file=paste(output.folder,"GCM.",GCMi,".Tempmax",
  #                                  min(list.years[[t]]),"-",max(list.years[[t]]),"_day",
  #                                  min(seqs[[i]]),"-",max(seqs[[t]]),".RData",sep=""))
  
  
  Sys.time()
  #rm(GCM.bias.corr.min,GCM.bias.corr.max,GCM.bias.corr,list.store.mins,list.store.maxs)
  rm(GCM.bias.corr)
}

stopCluster(cl)
