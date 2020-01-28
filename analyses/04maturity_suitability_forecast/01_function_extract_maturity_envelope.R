#'###################################################################
#' Function to extract reference climatic envelope during maturation:
#'
#' Started 20th January 2017
#' Ignacio Morales-Castilla & Elizabeth M. Wolkovich
#'###################################################################

## remove objects (should not be needed if the function is sourced)
#rm(list=ls())
#options(stringsAsFactors = FALSE)



#' Function to forecast phenological models 
#' 
#' @param dir.pheno - local folder containing simulated phenology for each variety at each site 
#' @param dir.climate -  local folder containing .nc daily climate data from BEST
#' @param dir.prec - local folder containing .nc precipitation data from 
#'                  (https://www.esrl.noaa.gov/psd/data/gridded/data.gpcc.html)          
#' @param varieties - identifier for a given variety (1-11 in our study) - this
#'        scalar is used to subset the object "list.var.coords" containing a list    
#'        with the coordinates for each variety from Anderson & Aryal (2015) Which 
#'         winegrape varieties are grown where? A global empirical picture.      
#' @return storing.envelope.each.gcm - list of dataframes for each GCM and the 
#'        selected variety, with coordinates where that variety is currently planted
#'        and 8 bioclimatic variables characterizing its climatic envelope. For a 
#'        description of the variables see Modelling Maturity section in Supplementary Materials 

extract.envelope<-function(dir.pheno,dir.climate,dir.prec,varieties){
  periods<-"1950-1980."
  past_pheno.dir<-dir(dir.pheno)
  #order.corr<-c(2,3,6,11,8,1,4,5,7,9,10) ## if needed to make sure varieties match order in varieties_parameters_best.csv  
  
  storing.envelope.each.gcm<-list()
  gcms.good<-strsplit(past_pheno.dir[which(grepl(periods,past_pheno.dir,fixed=T) )],".",fixed=T)
  
  for(gcm in gcms.good){
    ## get name ready for phenology
    GCM.target<-paste('GCM.',gcm,'.',sep='')
    name.past.pheno<-past_pheno.dir[which(grepl(GCM.target,past_pheno.dir,fixed=T) &
                                            grepl(periods,past_pheno.dir,fixed=T) )]
    
    
    ## load phenology data
    load(paste(dir.pheno,name.past.pheno,sep='/')) # year/variety/stage
    
    
    # defining variety data
    var.i.data<-list.var.coords[[varieties]]
    var.i.coords<-var.i.data[which(!is.na(var.i.data$latitude)),c(14,13)]
    
    coord.shp<-as.data.frame(var.i.data[which(!is.na(var.i.data$latitude)),c(14,13)])
    coord.shp$longitude<-as.numeric(coord.shp$longitude)
    coord.shp$latitude<-as.numeric(coord.shp$latitude)
    north.locs<-which(coord.shp$latitude>0)
    var.ids<-as.numeric(rownames(coord.shp))
    
    ## obtaining starting date (we will be counting since veraison occurs)
    veraisonlist<-list()
    for (yea in 10:length(GCM.future.proj)){ #yea=1
      
      ## average veraison dates
      ################
      # extract values
      veraisonlist[[yea]]<-extract(GCM.future.proj[[yea]][[order.corr[varieties]]][[3]],coord.shp[,1:2])
    }  
    
    ## sumarize veraison dates
    veraison.allyears<-list()
    for(i in 1:nrow(coord.shp)){
      print(i)
      for(j in 10:length(GCM.future.proj)){
        if(j==10){
          veraison.allyears[[i]]<-veraisonlist[[j]][[i]]
        }
        if(j>10){
          veraison.allyears[[i]]<-c(veraison.allyears[[i]],veraisonlist[[j]][[i]])
        }
        
      }
      
    }
    
    veraison.date<-unlist(lapply(veraison.allyears,mean,na.rm=T))
    coord.shp$veraison.date.corr<-round(ifelse(seq(1,nrow(coord.shp),1)%in%north.locs,
                                               veraison.date-365,veraison.date),0)
    coord.shp$veraison.date.corr<-ifelse(coord.shp$veraison.date.corr>365,
                                         coord.shp$veraison.date.corr-365,coord.shp$veraison.date.corr)
    
    coord.shp<-coord.shp[which(!is.na(coord.shp$veraison.date.corr)),]
    
    # transform into date (use of 1970 as a dummy year)
    ver.julian<-as.Date(paste(coord.shp$veraison.date.corr,"-1970",sep=""),format="%j-%Y")
    
    
    ## precipitation -- loading data 
    # call names for GPCC
    name.prec<-"precip.mon.total.v7.nc"
    GPCC.set<-paste(dir.prec,name.prec,sep='')
    #nc.GPCC <- nc_open(GPCC.set)
    #variabs <- names(nc.GPCC$var)
    # open as brick
    precip.GPCC <- brick(GPCC.set, varname="precip", layer="time")
    names.nc<-names(precip.GPCC)
    
    ## temperatures -- loading data
    ## load avg BEST climatologies
    BEST.mean.clim.min<-brick("BEST.mean.clim_min.nc")
    BEST.mean.clim.max<-brick("BEST.mean.clim_max.nc")
    
    storing.temp.site.i<-list()
    # dates of starting and end of data extraction
    for(prec.i in 1:nrow(coord.shp)){ #prec.i=1
      subs.var.i<-coord.shp[prec.i,]
      
      ## defining dates from veraison
      start.date<-ver.julian[prec.i]
      start.month<-as.numeric(format(start.date,"%m"))
      ndays<-ifelse(start.month%in%c(1,3,5,7,8,10,12),31,30)
      weight.start.month<-(ndays-(as.numeric(format(start.date,"%d")))+1)/ndays
      end.date<-start.date+44
      end.month<-as.numeric(format(end.date,"%m"))
      ndays.end<-ifelse(end.month%in%c(1,3,5,7,8,10,12),31,30)
      weight.end.month<-as.numeric(format(end.date,"%d"))/ndays.end
      
      if(end.month-start.month>1 | end.month-start.month<0){
        if(end.month==1){
          middle.month<-11
          weight.middle<-1
        } 
        if(end.month==2) {
          middle.month<-12
          weight.middle<-1
        }
        if(!end.month%in%1:2) {
          middle.month<-end.month-1
          weight.middle<-1
        }
      } else {
        weight.middle<-0
        middle.month<-NA
      }
      
      
      ## loop for precipitation
      storing.prec.site.i<-list()
      cummul.precs<-list()
      for(date in seq(1950,1980,1)){ 
        print(paste("GCM-",gcm, prec.i,date))
        layer.month.a<-which(grepl(date,names.nc))[start.month]
        layer.month.b<-which(grepl(date,names.nc))[end.month]
        
        # extract precipitation
        vals.month.a<-extract(rotate(precip.GPCC[[layer.month.a]]),subs.var.i[,1:2])
        vals.month.b<-extract(rotate(precip.GPCC[[layer.month.b]]),subs.var.i[,1:2])
        
        if(!is.na(middle.month)){
          layer.month.c<-which(grepl(date,names.nc))[middle.month]
          vals.month.c<-extract(rotate(precip.GPCC[[layer.month.c]]),subs.var.i[,1:2])
        } else {
          vals.month.c<-rep(NA,length(vals.month.a))
        }
        
        storing.prec.site.i[[which(seq(1950,1980,1)==date)]]<-c(vals.month.a,vals.month.b,vals.month.c)
        cummul.precs[[which(seq(1950,1980,1)==date)]]<-sum(c(mean(vals.month.a,na.rm=T),
                                                             mean(vals.month.b,na.rm=T),
                                                             mean(vals.month.c,na.rm=T))*c(weight.start.month,
                                                                                           weight.end.month,
                                                                                           weight.middle),na.rm=T)
        
      }
      
      
      #completing dataset with extracted values for prec
      weights.i=c(rep(weight.start.month,length(vals.month.a)),
                  rep(weight.end.month,length(vals.month.b)),
                  rep(weight.middle,length(vals.month.c)))
      coord.shp$lim.low.prec[prec.i]<-quantile(unlist(storing.prec.site.i),c(0.05,0.95),na.rm=T)[1]
      coord.shp$lim.up.prec[prec.i]<-quantile(unlist(storing.prec.site.i),c(0.05,0.95),na.rm=T)[2]
      coord.shp$w.mean.prec[prec.i]<-weighted.mean(unlist(storing.prec.site.i),w=rep(weights.i,31),na.rm=T)
      coord.shp$w.sum.prec[prec.i]<-mean(unlist(cummul.precs),na.rm=T)
      coord.shp$w.sum.prec.low[prec.i]<-quantile(unlist(cummul.precs),c(0.05,0.95),na.rm=T)[1]
      coord.shp$w.sum.prec.up[prec.i]<-quantile(unlist(cummul.precs),c(0.05,0.95),na.rm=T)[2]
      
      ## extracting and storing temperature data
      ncells<-length(extract(BEST.mean.clim.min[[1]],subs.var.i[,1:2]))
      store.temps<-array(NA,dim=c(45,3)) #days/ncells/min-max-mean
      
      dates.extract<-as.numeric(format(seq(start.date,end.date,1),"%j"))
      store.temps[,1]<-extract(subset(BEST.mean.clim.min,dates.extract),subs.var.i[,1:2])
      store.temps[,2]<-extract(subset(BEST.mean.clim.max,dates.extract),subs.var.i[,1:2])
      
      if(ncells>1){
        store.temps[,3]<-apply(store.temps[,,1:2],1:2,mean,na.rm=T)
      } else {
        store.temps[,3]<-rowMeans(store.temps[,1:2],na.rm=T)
      }
      
      #completing dataset with extracted values for prec
      storing.temp.site.i[[prec.i]]<-store.temps
      
      
    }
    
    
    # process temp data to get the rest of variables:
    coord.shp$gdd=coord.shp$above40=coord.shp$below10=rep(NA,nrow(coord.shp))
    coord.shp$mintemp=coord.shp$maxtemp=coord.shp$gddover25=rep(NA,nrow(coord.shp))
    for (i in 1:nrow(coord.shp)){#i=2
      
      # GDD
      if(!is.null(dim(storing.temp.site.i[[i]]))){  
        daysabove<-which(storing.temp.site.i[[i]][,3]>10)
        coord.shp$gdd[i]<-sum(storing.temp.site.i[[i]][daysabove,3]-10,na.rm=T)
      } else {
        daysabove<-which(storing.temp.site.i[[i]][,3]>10)
        coord.shp$gdd[i]<-sum(storing.temp.site.i[[i]][daysabove,3]-10,na.rm=T)
      }
      
      #  number of days above 40
      if(!is.null(dim(storing.temp.site.i[[i]]))){  
        coord.shp$above40[i]<-sum(storing.temp.site.i[[i]][,2]>40,na.rm=T)
      } else {
        coord.shp$above40[i]<-sum(storing.temp.site.i[[i]][,2]>40,na.rm=T)
      }
      
      #  number of days below 10
      if(!is.null(dim(storing.temp.site.i[[i]]))){  
        coord.shp$below10[i]<-sum(storing.temp.site.i[[i]][,1]<10,na.rm=T)
      } else {
        coord.shp$below10[i]<-sum(storing.temp.site.i[[i]][,1]<10,na.rm=T)
      }
      
      #  number of gdd above 25 (optimum temperature in the models)
      if(!is.null(dim(storing.temp.site.i[[i]]))){  
        daysabove<-which(storing.temp.site.i[[i]][,3]>25)
        coord.shp$gddover25[i]<-sum(storing.temp.site.i[[i]][daysabove,3]-25,na.rm=T)
      } else {
        daysabove<-which(storing.temp.site.i[[i]][,3]>25)
        coord.shp$gddover25[i]<-sum(storing.temp.site.i[[i]][daysabove,3]-25,na.rm=T)
      }
      
      #  minimum temps 
      if(!is.null(dim(storing.temp.site.i[[i]]))){  
        coord.shp$mintemp[i]<-min(storing.temp.site.i[[i]][,1],na.rm=T)
      } else {
        coord.shp$mintemp[i]<-min(storing.temp.site.i[[i]][,1],na.rm=T)
      }
      
      #  maximum temps 
      if(!is.null(dim(storing.temp.site.i[[i]]))){  
        coord.shp$maxtemp[i]<-max(storing.temp.site.i[[i]][,2],na.rm=T)
      } else {
        coord.shp$maxtemp[i]<-max(storing.temp.site.i[[i]][,2],na.rm=T)
      }
      
      
    }
    
    
    ## saving results fore each variety
    #storing.envelope.each.var[[vars]]<-coord.shp
    storing.envelope.each.gcm[[gcm]]<-coord.shp
    
  }
  
  return(storing.envelope.each.gcm)
  
}
