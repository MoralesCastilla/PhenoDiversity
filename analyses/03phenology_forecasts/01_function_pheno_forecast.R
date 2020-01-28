#'######################################################
#' Function to forecast phenological models:
#'
#' Started 20th January 2017
#' Ignacio Morales-Castilla & Elizabeth M. Wolkovich
#'######################################################

## remove objects (should not be needed if the function is sourced)
#rm(list=ls())
#options(stringsAsFactors = FALSE)



#' Function to forecast phenological models 
#' 
#' @param Year.pheno - Year for which forecasting phenology 
#' @param target.folder -  local folder containing .nc files with maximum and minimum daily temperatures
#'          resulting from file 04_save_GCMcorrected_bias_as_nc.R
#' @param GCMi - number of simulation from the climate forecast data from the LENS model (1-30)          
#' @param period - each year period in vector list.years: value 1 = 1950-2005; 
#'          value 2 = 2006-2080; value 3 = 2081-2100          
#' @param varieties - data frame with parameters for varieties to model from PMP v5.0 
#'        see file phenology_parameterization/varieties_parameters_best.csv         
#' @details pheno.seq.fx is used within the function. It can be found in 01_functions_pheno_models.R
#' @return Phenol.year.i - list of rasterStacks: containing values in DOY for Budbreak, Flowering and Veraison
#'

## function to project phenology for a given year
Future.pheno.fx<-function(Year.pheno,target.folder,GCMi,period,varieties){
  
  ## names of layers to load
  dfiles<-dir(target.folder)
  year0_max<-paste(min(period),max(period),sep="-")
  if(GCMi<10){GCM.name<-paste("GCM.",0,GCMi,sep="")} else {GCM.name<-paste("GCM",GCMi,sep=".")}
  
  to.load.max<-paste(target.folder,GCM.name,".Tmax",year0_max,".nc",sep="")
  to.load.min<-paste(target.folder,GCM.name,".Tmin",year0_max,".nc",sep="")
  
  ## loading layers
  daily.temp.GCM.min <- brick(to.load.min)
  daily.temp.GCM.max <- brick(to.load.max)
  
  ## re-naming
  names(daily.temp.GCM.min)<-rename.GCM(daily.temp.GCM.min,period)
  names(daily.temp.GCM.max)<-rename.GCM(daily.temp.GCM.max,period)
  
  
  #for(i in 2:length(period)){}
  ## subset first years
  #Year.pheno=1951
  which.year<-which(period==Year.pheno)
  subs.GCM.mins<-daily.temp.GCM.min[[c(seq(which.year-1,nlayers(daily.temp.GCM.min),length(period)),
                                       seq(which.year,nlayers(daily.temp.GCM.min),length(period)))]]
  subs.GCM.max<-daily.temp.GCM.max[[c(seq(which.year-1,nlayers(daily.temp.GCM.max),length(period)),
                                      seq(which.year,nlayers(daily.temp.GCM.max),length(period)))]]
  
  ## get coords and values
  mat.min<-rasterToPoints(subs.GCM.mins,progress="text")
  mat.max<-rasterToPoints(subs.GCM.max,progress="text")
  ## 4mins to extract points
  
  X<-list(mat.max,mat.min)
  Y<-do.call(cbind,X)
  Y<-array(Y,dim=c(dim(X[[1]]),length(X)))
  mat.avgs<-apply(Y,c(1,2),mean,na.rm=T)
  ## 5+mins to compute avg.Temp
  mats<-as.data.frame(t(mat.avgs))
  rm(mat.min,mat.max,mat.avgs)
  
  ## apply function each 
  pheno.stage<-list()
  #for(v in 1:2){#v=1
    
  for(v in 1:nrow(varieties)){#v=1
  print(v)
  pheno.stage[[v]]<-t(sapply(mats,pheno.seq.fx,variety=varieties[v,]))
  ## 8mins each
  }
  
  ## translate into rasters
  ras.BB=ras.FL=ras.VER<-subs.GCM.mins[[1]]
  Phenol.year.i<-list()
  for(ras in 1:length(pheno.stage)){
  ras.BB[which(!is.na(subs.GCM.mins[[1]][]))]<-pheno.stage[[ras]][,1]
  ras.FL[which(!is.na(subs.GCM.mins[[1]][]))]<-pheno.stage[[ras]][,2]
  ras.VER[which(!is.na(subs.GCM.mins[[1]][]))]<-pheno.stage[[ras]][,3]
  
  ## stack rasters
  Phenol.year.i[[ras]]<-stack(ras.BB,ras.FL,ras.VER)
  
  }
  
  
  return(Phenol.year.i)
  
}



