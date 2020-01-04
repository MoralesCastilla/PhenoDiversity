##############################################################################################################
# Function to:
#' * Extract  BEST climate data for each location given a list of coordinates and dates
#'
#' Data is first downscaled and interpolated to match GCM data
#'
#'  Ignacio Morales-Castilla & Elizabeth M. Wolkovich
#'  Started 02 Nov 2016
##############################################################################################################

## remove objects
rm(list=ls())

## load packages
packs.to.extract <- list('raster','ncdf4','maptools','sp','foreach')
lapply(packs.to.extract,require, character.only=T)


#' function to extract data from BEST (download from http://berkeleyearth.org/data/)
#'
#' @param dir.BEST - directory with BEST avg daily data
#' @param dir.BEST.max - directory with BEST max daily data
#' @param dir.BEST.min - directory with BEST min daily data
#' @param GCM.day.i - raster file from GCM to interpolate data
#' @param coastal.rast - raster file from BEST with coastal cells (interpolated)      
#' @param coords.site - data.frame with sitename and lon lat coords for each site format:
#' station_code  Lon  Lat   
#' 1
#' ...
#' n
#' 
#' @param site.year - data.frames with site, lat, year and doy columns format:
#'  SITE  Year
#'  1
#'  ...
#'  n
#'
#' @return climate.data - data frame with Tavg, Tmax and Tmin for each day in each location
##### Parallel version of the function: coded to work as a function of ndec (number of decade)


Extract.BEST.fx.parallel <- function(ndec,dir.BEST,dir.BEST.max,dir.BEST.min,
                                   coords.site,site.year,GCM.day.i,coast.rast,land.rast){
  
  
  # Define where BEST data is downloaded
  dfiles <- dir(dir.BEST)
  dfiles.max <- dir(dir.BEST.max)
  dfiles.min <- dir(dir.BEST.min)
  years.list <- sort(unique(site.year$Year))
  
  # generate a data frame to store data
  init.dims.to.store <- expand.grid(coords.site$station_code,site.year$Year,seq(1,365,1))
  to.remove <- which(!apply(init.dims.to.store[,1:2],1,paste, collapse=" ")%in%apply(site.year[,1:2],1,paste, collapse=" "))
  
  if(length(to.remove)!=0){
    init.dims.to.store <- init.dims.to.store[-to.remove,]}
  lat <- coords.site[match(init.dims.to.store[,1],coords.site[,1]),3]
  
  climate.data.store <- as.data.frame(array(NA,dim=c(nrow(init.dims.to.store),7)))
  colnames(climate.data.store) <- c("Site","Lat","Year","DOY","Tmean","Tmin","Tmax")
  climate.data.store[,c("Site","Year","DOY")] <- init.dims.to.store
  climate.data.store[,c("Lat")] <- lat
  climate.data.store <- climate.data.store[!duplicated(climate.data.store),]
  
  # how many decades to process?
  ndecades <- unique((years.list %/% 10) * 10)
  
  BEST.set <- paste(dir.BEST,dfiles[which(grepl(ndecades[ndec],dfiles))],sep='/')
  BEST.set.max <- paste(dir.BEST.max,dfiles.max[which(grepl(ndecades[ndec],dfiles.max))],sep='/')
  BEST.set.min <- paste(dir.BEST.min,dfiles.min[which(grepl(ndecades[ndec],dfiles.min))],sep='/')
  nc.BEST  <-  nc_open(BEST.set)
  nc.BEST.max  <-  nc_open(BEST.set.max)
  var  <-  names(nc.BEST$var)
  
  # open as brick
  anomalies.BEST  <-  brick(BEST.set, varname=var[7], layer="time")
  climatology.BEST <- brick(BEST.set,varname=var[8],layer="time")
  
  anomalies.BEST.max  <-  brick(BEST.set.max, varname=var[7], layer="time")
  climatology.BEST.max <- brick(BEST.set.max,varname=var[8],layer="time")
  
  anomalies.BEST.min  <-  brick(BEST.set.min, varname=var[7], layer="time")
  climatology.BEST.min <- brick(BEST.set.min,varname=var[8],layer="time")
  
  
  # naming layers in bricks
  doy <-  ncvar_get(nc.BEST,"day_of_year")
  year <-  ncvar_get(nc.BEST,"year")
  month <-  ncvar_get(nc.BEST,"month")
  day <- ncvar_get(nc.BEST,"day")
  timevec <- paste(year,month,day,sep=" ")
  
  doy.max <-  ncvar_get(nc.BEST.max,"day_of_year")
  year.max <-  ncvar_get(nc.BEST.max,"year")
  month.max <-  ncvar_get(nc.BEST.max,"month")
  day.max <- ncvar_get(nc.BEST.max,"day")
  timevec.max <- paste(year.max,month.max,day.max,sep=" ")
  years.vec <- as.character(strptime(timevec,format="%Y %m %d"))
  years.vec.max <- as.character(strptime(timevec.max,format="%Y %m %d"))
  
  names(climatology.BEST)[60:365] <- paste("X",seq(61,366,1),sep="") ## correcting for leap years
  climdates <- format(as.Date(names(climatology.BEST), "X%j"),"%m.%d")
  names(anomalies.BEST)=as.vector(years.vec)
  names(anomalies.BEST.max)=names(anomalies.BEST.min) <- as.vector(years.vec.max)
  names(climatology.BEST)=names(climatology.BEST.max)=names(climatology.BEST.min) <- as.vector(climdates)
  
  # years to run analyses for
  unique.years <- sort(unique(year))
  years.decade <- unique.years[which(unique.years%in%years.list)]
  
  for(j in length(years.decade):length(years.decade)){
    #for(j in 1:length(years.decade)){
    year.i <- years.decade[j]
    anom.best.year <- anomalies.BEST[[which(year%in%year.i)]]
    anom.best.year.max <- anomalies.BEST.max[[which(year.max%in%year.i)]]
    anom.best.year.min <- anomalies.BEST.min[[which(year.max%in%year.i)]]
    sites.that.year <- site.year[which(site.year[,2] == year.i),]
    coords.sites.that.year <- coords.site[which(coords.site[,1]%in%sites.that.year$SITE),]
    is.leap <- nlayers(anom.best.year)
    
    if(is.leap == 366){
      anom.best.year <- anom.best.year[[-60]]
      anom.best.year.max <- anom.best.year.max[[-60]]
      anom.best.year.min <- anom.best.year.min[[-60]]
    }
    
    for(k in 325:nlayers(anom.best.year.max)){
      
      #for(k in 1:nlayers(anom.best.year.max)){
      
      print(paste("decade",ndec,"year",year.i,"day",k))
      ## comparing date
      dayi <- names(anom.best.year)[k]
      anom.day.i <- anom.best.year[[dayi]] #raster with temp anomalies.BEST for day i
      anom.day.i.max <- anom.best.year.max[[dayi]] #raster with temp anomalies.BEST for day i
      anom.day.i.min <- anom.best.year.min[[dayi]] #raster with temp anomalies.BEST for day i
      
      anom.day.compare <- paste("X",format(as.Date(dayi, "X%Y.%m.%d"),"%m.%d"),sep="")
      clim.day.i <- climatology.BEST[[anom.day.compare]]
      clim.day.i.max <- climatology.BEST.max[[anom.day.compare]]
      clim.day.i.min <- climatology.BEST.min[[anom.day.compare]]
      doy.i <- as.numeric(format(as.Date(dayi, "X%Y.%m.%d"),"%j"))
      
      if(is.leap == 366 & doy.i>60){
        doy.i <- doy.i-1
      }
      
      ## computing that date's temperature (climatology + anomaly)
      temp.day.i <- sum(clim.day.i,anom.day.i,na.rm=T)
      temp.day.i.max <- sum(clim.day.i.max,anom.day.i.max,na.rm=T)
      temp.day.i.min <- sum(clim.day.i.min,anom.day.i.min,na.rm=T)
      
      # interpolating to same resolution as GCM 
      temp.interpolated <- projectRaster(from=temp.day.i,to=GCM.day.i,method="bilinear")
      temp.interpolated.max <- projectRaster(from=temp.day.i.max,to=GCM.day.i,method="bilinear")
      temp.interpolated.min <- projectRaster(from=temp.day.i.min,to=GCM.day.i,method="bilinear")
      
      #which location is coastal based on land and on GCM
      #is.coastal <- extract(coast.rast,coords.sites.that.year[,2:3])
      is.coastal <- extract(land.rast,coords.sites.that.year[,2:3])
      
      #only for coastal cells 
      if(sum(is.na(is.coastal))<length(is.coastal)){
        coast.coords <- coords.sites.that.year[which(is.na(is.coastal)),2:3]
        
        for(costa in 1:nrow(coast.coords)){
          dist.coords <- distanceFromPoints(land.rast, coast.coords[costa,])
          #LonLat(coast.coords[costa,])
          dist.coords[is.na(land.rast)] <- NA
          coast.neib.cell <- which.min(dist.coords)
          
          ## extract climate for that day and locations represented
          position.sites <- which(climate.data.store$Site%in%coords.sites.that.year[which(is.na(is.coastal)),1] & climate.data.store$Year == year.i & 
                                  climate.data.store$DOY == doy.i)[costa]
          
          ### needs correction here
          #sites.in.matrix <-  climate.data.BEST.PMP[position.sites,1]
          climate.data.store[position.sites,5] <- temp.interpolated[rowColFromCell(dist.coords, coast.neib.cell)]
          climate.data.store[position.sites,6] <- temp.interpolated.min[rowColFromCell(dist.coords, coast.neib.cell)]
          climate.data.store[position.sites,7] <- temp.interpolated.max[rowColFromCell(dist.coords, coast.neib.cell)]
          
        }
        
        
      }
      
      
      ## extract climate for that day and locations represented
      position.sites <- which(climate.data.store$Site%in%coords.sites.that.year[which(!is.na(is.coastal)),1] & climate.data.store$Year == year.i & 
                              climate.data.store$DOY == doy.i)
      #sites.in.matrix <-  climate.data.BEST.PMP[position.sites,1]
      climate.data.store[position.sites,5] <- extract(temp.interpolated,coords.sites.that.year[which(!is.na(is.coastal)),2:3])
      climate.data.store[position.sites,6] <- extract(temp.interpolated.min,coords.sites.that.year[which(!is.na(is.coastal)),2:3])
      climate.data.store[position.sites,7] <- extract(temp.interpolated.max,coords.sites.that.year[which(!is.na(is.coastal)),2:3])
      
    }
    
    
  }
  
  
  
  return(climate.data.store[which(!is.na(climate.data.store[,5])),])
  
  
}





