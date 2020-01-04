#'#############################################################################################################
#' Script and functions to fit phenological models:
#' 
#' functions include: 
#' Smoothed Utah model - for chilling units 
#' Chuine model - for chilling units as in Chuine et al. 2000
#' Wang & Engel - for forcing 
#'
#'  Ignacio Morales-Castilla, Iñaki García de Cortazar-Atauri, Elizabeth M. Wolkovich et al.
#'  19 Oct 2016
#'#############################################################################################################


#### housekeeping ####
#rm(list=ls())

#### load auxiliary functions ####
is.leapyear=function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

## Model performance test - Root Mean Square error
#' @param Xa - observed vector
#' @param Xb - predicted vector
#' @return RMSE 

RMSE  <-  function(Xa,Xb){
  N <- length(Xa)
  rmse <- (sum((Xa-Xb)^2)/N)^0.5
  return(rmse)
}

#### functions ####

#' Chuine model Unichill. Unimodal and symetrical function. 2 parts - compute the function and get date
#'   After Cstar date (t1) is obtained, 
#' Part 1a - Fxchuine
#' @param a - width window funciton != null
#' @param b - sharpness of the response curve and asymetry
#' @param c - mid-response value    
#' @param Tavg - vector with average daily temperature
#' @return Chill.units - number of chill units produced each day
#' @return Accummulated.chill - accumulated chill units
#' 
#' Part 1b - SmoothUtah # alternative to Fxchuine 
#' assumes chilling to occur only under certain temperatures (there's negative chilling on warm days)
#' @param Tm1 - sharptness of decrease in cold efficiency for but endodormacy
#' @param Topt - optimal daily temp for which 1 chilling unit is accumulated every day
#' @param Tn2 - mid-response (temperature above Topt wich has half the efficciency to induce endodormancy)    
#' @param minn - how negative is the impact of too high temperatures?
#' @param Year.previous - year previous to that modelled
#' @param Year.target - year modelled
#' @return dataframe binding Chill.units and Accummulated.chill, daily
#' 
#' Part 2 - GetCstardate 
#' @param Cstar - sum of chilling units from day 0 till day of Budbrake (BBdate)
#' @param Fstar - sum of forcing units from day Budbrake
#' @param Tavg - vector with average daily temperature
#' @param day1 - day when Cstar is reached
#' @return Cstardate - date of bud brake (Cstar reached)  
#' @return Cstardate.julian - BBdate in julian dates


Fxchuine <- function(a,b,c,Tavg,Year.previous,Year.target,Hemisphere){
  
  if(Hemisphere=="NONE"|Hemisphere=="NORTH"){
    if(is.leapyear(Year.previous)){
      day1 <- 213
    } else {
      day1 <- 214
    }
    } else { ## for Southern Hemisphere locations
    if(is.leapyear(Year.previous)){
      day1 <- 213-182
    } else {
      day1 <- 214-182
    }
    }
  Chill.units=1/(1+exp((a*((Tavg[day1:nrow(Tavg),5]-c)^2))+(b*(Tavg[day1:nrow(Tavg),5]-c))))
  Chill.accum=cumsum(Chill.units)
  
  return(data.frame(cbind(Chuill.units=Chill.units,Accummulated.chill=Chill.accum)))
}


SmoothUtah <- function(Tm1,Topt,Tn2,minn,Tavg,Year.previous,Year.target,Hemisphere){
  
  if(Hemisphere=="NONE"|Hemisphere=="NORTH"){
    if(is.leapyear(Year.previous)){
      day1 <- 213
    } else {
      day1 <- 214
    }
  } 
  if(Hemisphere=="SOUTH"){
  if(is.leapyear(Year.previous)){
      day1 <- 213-182
    } else {
      day1 <- 214-182
    }
  }
  
  Chill.units <- as.vector(array(0,dim=c(nrow(Tavg),1)))
  for(i in day1:nrow(Tavg)){
  if(Tavg[i,5]<Tm1){ # note this is incorrect in original paper
    Chill.units[i]=1/(1+exp(-4*((Tavg[i,5]-Tm1)/(Topt-Tm1))))
  
    } else {
  
    if(Tavg[i,5]>=Tm1 & Tavg[i,5]<Topt){
      Chill.units[i]=1+(-0.5/((Tm1-Topt)^2))*(Tavg[i,5]-Topt)^2
  
    } else {
  
    if(Tavg[i,5]>Topt & Tavg[i,5]<Tn2){
    Chill.units[i]=1-(1-minn)*(((Tavg[i,5]-Topt)^2)/(2*((Tn2-Topt)^2)))
  
    } else {
  
    if(Tavg[i,5]>Tn2){
    Chill.units[i]=minn+((1-minn)/(1+exp(-4*((Tn2-Tavg[i,5])/(Tn2-Topt)))))
  }
  }
  }
  }
}
  Chill.accum=cumsum(Chill.units)
  
  return(data.frame(cbind(Chuill.units=Chill.units,Accummulated.chill=Chill.accum)))
}


GetCstardate <- function(Cstar,Tavg,result.chilling){ ## interpolated version
  
  #day1 <- which(result.chilling[,2]>=Cstar)[1]
  day1=approx(result.chilling[,2],seq(1,nrow(result.chilling)),xout=Cstar)$y
  
  return(list(Cstardate=Tavg[round(day1,0),3:4],Cstardate=day1))
  
}



#' Wang and Engel model. This model computes the rate of phenological 
#' development in response to temperature using Topt bounded by Tmin and Tmax
#' (below and above which there is no Temp action on the plant) in 3 parts
#'
#' Part 1- Alpha - based on fixed parameters
#' @param Tmin - Minimum temperature
#' @param Tmax - Maximum temperature
#' @param Topt - Optimum temperature
#' @return Alpha

Alphafx  <-  function(Tmin, Tmax, Topt){
  Alpha  <-  log10(2)/(log10((Tmax-Tmin)/(Topt-Tmin)))
  return(Alpha)
}

#' Part 2- WangEngelfx - computes temperature action across all days of the year
#' @param Tmin - Minimum temperature
#' @param Tmax - Maximum temperature
#' @param Topt - Optimum temperature
#' @param Alpha - alpha parameter 
#' @param Tavg - daily average temperature
#' @return cTt - data frame with Temperature action on developement and its daily accumulated values - 
#' 
WangEngelfx  <-  function(Tmin, Tmax, Topt, Alpha, Tavg){
  cTt <- array(NA, dim=c(length(Tavg),2))
  colnames(cTt) <- c("Temp.action","Accum.Temp.action")
  for(i in 1:length(Tavg)){
    if (Tmin<Tavg[i] & Tavg[i]<Tmax){
      cTt[i,1]  <-  (2*(Tavg[i]-Tmin)^Alpha*(Topt-Tmin)^Alpha-(Tavg[i]-Tmin)^(2*Alpha))/((Topt-Tmin)^(2*Alpha))
    } 
    if (Tavg[i]<=Tmin | Tavg[i]>=Tmax){
      cTt[i,1]  <- 0
    }
  }
  cTt[,2] <- cumsum(cTt[,1])
  
  return(cTt)
}


#' Part 3- GetBBdate, GetFlowdate and GetVeraisdate
#' @param Tavg - daily average temperatures 
#' @param result.WangEng - data frame with temperature action daily and accumulated
#' @param dateBB - day when budbrake happens 
#' @param t0 - day of when temperature action starts for flowering #same as Ndb?
#' @param Fstar - Accummulated temperature action from start date till flow/verais date  
#' @return Flowdate - date of flowering
#' @return Flowdate.julian - 
#' @return Veraisdate - date of veraison
#' @return Veraisdate.julian


GetBBdate <- function(FstarBB,Tavg,result.WangEng,Cstardate){
  Cstardate <- Cstardate[[2]]
  
  accum.Funits <- cumsum(result.WangEng[Cstardate:nrow(result.WangEng),1])
  
  day.BB <- which(accum.Funits>=FstarBB)[1]+Cstardate
  
  return(list(BBdate=Tavg[day.BB,3:4],BBdate.julian=day.BB))
  
}


GetFlowdate <- function(FstarFl,Tavg,result.WangEng,dateBB){
  dateBB <- dateBB[[2]]
  
  accum.Funits <- cumsum(result.WangEng[dateBB:nrow(result.WangEng),1])
  
  day.Flow <- which(accum.Funits>=FstarFl)[1]+dateBB
  
  return(list(Flowdate=Tavg[day.Flow,3:4],Flowdate.julian=day.Flow))
  
}


GetVeraisdate <- function(FstarVer,Tavg,result.WangEng,dateFlow){
  dateFlow <- dateFlow[[2]]
  
  accum.Funits <- cumsum(result.WangEng[dateFlow:nrow(result.WangEng),1])
  
  day.Verais <- which(accum.Funits>=FstarVer)[1]+dateFlow
  
  return(list(Veraisdate=Tavg[day.Verais,3:4],Veraisdate.julian=day.Verais))
  
}




  
## Function to fit models sequentially and estimate yearly phenology for 
## varieties and sites
#' @param climate.data - data frame containing climate data with 7 columns:
#'        Site | Latitude | Year | DOY | TAV | TMIN | TMAX
#' @param variety.data - data frame with 7 columns to be used as model to store results:
#'        Variety | Site | Provenance | Year | BB | FLOR | VER
#' @param variety.parameters - fixed parameters for each variety in:
#'        ../data/phenology_parameterization/varieties_parameters_best.csv
#' @param UseBestparam - whether or not using BEST phenological parameters
#' @param chillingfx - what chilling function to use to simulate dormancy
#' @return  Results.phenology - data frame (variety, site, year, BBdate, Fldate, Verdate)
#' 
yearly.phenology.dates <- function(climate.data,variety.data,variety.parameters,
                                 UseBestparam,chillingfx=c("Chuine","Utah")){

  ## format/handle input data
  list.varieties <- unique(variety.data$Variety)
  site.var <- variety.data[,1:2]
  list.site.var <- site.var[!duplicated(site.var),]
  list.site.year <- climate.data[!duplicated(climate.data[,c(1,3)]),c(1,3)]
  #list.varieties[6] <- "Cabernet-Sauvignon"
  list.sites <- as.character(unique(variety.data$Site))
  list.years <- sort(unique(variety.data$Year))
  list.variety_parameter <- variety.parameters$Variety[variety.parameters$Variety%in%list.varieties]
  
  ## array to store results
  Results.phenology <- as.data.frame(array(NA,dim=dim(variety.data)))
  Results.phenology[,1:4] <- variety.data[,1:4]
  colnames(Results.phenology) <- colnames(variety.data)
  
  
  ## loop computation for each year, variety and site 
  for(j in 1:length(list.sites)){
    Site <- list.sites[j]
    hemisphere <- ifelse("Hemisphere"%in%colnames(variety.data),
                       variety.data[which(variety.data$Site==Site)[1],"Hemisphere"],"NONE")
    for(k in 1:length(list.variety_parameter)){
      Variety <- list.variety_parameter[k]
      
      if(length(which(list.site.var[,1]==Variety & list.site.var[,2]==Site))!=0){
        
        ## temperature vector for the subsetted year, site and variety
        temps.avg <- climate.data[which(climate.data$SITE==Site),]
        head(climate.data)
        years.vector.met <- list.site.year[which(list.site.year[,1]==Site),2]
        years.vector <- variety.data[which(variety.data[,1]==Variety & variety.data[,2]==Site),4]
        
        for(i in 1:length(years.vector)){
          #print(paste("site: ",j,"; variety: ", k,"; year: ",i,sep=""))
          Year.target <- years.vector[i]
          if(which(years.vector.met==Year.target)>1){
            Year.previous <- years.vector.met[which(years.vector.met==Year.target)-1]
            temps.avgs <- temps.avg[which(temps.avg$Year==Year.previous | temps.avg$Year==Year.target),]
            parameters <- variety.parameters[which(variety.parameters[,"Variety"]==Variety),2:23]
            
            if(UseBestparam==F & sum(is.na(parameters))==0 & !-9999.0%in%temps.avgs[,5]){
              ## 2) get Budburst date (based on Chuine model)
              
              if(chillingfx=="Chuine"){
              Chillunits <- Fxchuine(parameters$a,parameters$b,parameters$c,temps.avgs,Year.previous,Year.target,hemisphere)
              }
              
              if(chillingfx=="Utah"){
              Chillunits <- SmoothUtah(parameters$Tm1,parameters$Topt,parameters$Tn2,parameters$minn,temps.avgs,Year.previous,Year.target,hemisphere)
              }
              
              DateCstar <- GetCstardate(parameters$Cstar_BB,temps.avgs,Chillunits)
              
              ## 3) compute temperature action (forcing) units based on WE model
              alpha.BB <- Alphafx(parameters$Tmin,parameters$Tmax,parameters$Topt.1)
              ActionunitsBB <- WangEngelfx(parameters$Tmin,parameters$Tmax,parameters$Topt.1,alpha.BB,temps.avgs$TAV)
              
              alpha.Fl <- Alphafx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL)
              ActionunitsFl <- WangEngelfx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL,alpha.Fl,temps.avgs$TAV)
              
              alpha.Ver <- Alphafx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER)
              ActionunitsVer <- WangEngelfx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER,alpha.Ver,temps.avgs$TAV)
              
              ## 4) get BB, Flowering and Veraison dates
              if(!is.na(DateCstar[[2]])){
              DateBB <- GetBBdate(parameters$Fstar_BB,temps.avgs,ActionunitsBB,DateCstar)
              Flow.date <- GetFlowdate(parameters$Fstar_BB_Fl_obs,temps.avgs,ActionunitsFl,DateBB)
              Verais.date <- GetVeraisdate(parameters$Fstar_Fl_Ver_obs,temps.avgs,ActionunitsVer,Flow.date)
              
              row.to.store <- which(Results.phenology[,1]==Variety & Results.phenology[,2]==Site & Results.phenology[,4]==Year.target)
              Results.phenology[row.to.store,5] <- DateBB[[1]][1,2]
              Results.phenology[row.to.store,6] <- Flow.date[[1]][1,2]
              Results.phenology[row.to.store,7] <- Verais.date[[1]][1,2]
              }
            }
            
            ## if climate data comes from BEST
            if(UseBestparam==T & sum(is.na(parameters))==0 & !-9999.0%in%temps.avgs[,5]){
              ## 2) get Budburst date (based on Chuine model) 
              if(chillingfx=="Chuine"){
                Chillunits <- Fxchuine(parameters$a,parameters$b,parameters$c,temps.avgs,Year.previous,Year.target,hemisphere)
              }
              
              if(chillingfx=="Utah"){
                Chillunits <- SmoothUtah(parameters$Tm1,parameters$Topt,parameters$Tn2,parameters$minn,temps.avgs,Year.previous,Year.target,hemisphere)
              }
              
              DateCstar <- GetCstardate(parameters$Cstar_BB,temps.avgs,Chillunits)
              
              ## 3) compute temperature action (forcing) units based on WE model
              alpha.BB <- Alphafx(parameters$Tmin,parameters$Tmax,parameters$Topt.1)
              ActionunitsBB <- WangEngelfx(parameters$Tmin,parameters$Tmax,parameters$Topt.1,alpha.BB,temps.avgs$TAV)
              
              alpha.Fl <- Alphafx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL)
              ActionunitsFl <- WangEngelfx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL,alpha.Fl,temps.avgs$TAV)
              
              alpha.Ver <- Alphafx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER)
              ActionunitsVer <- WangEngelfx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER,alpha.Ver,temps.avgs$TAV)
              
              ## 4) get BB, Flowering and Veraison dates 
              if(!is.na(DateCstar[[2]])){
              DateBB <- GetBBdate(parameters$Fstar_BB,temps.avgs,ActionunitsBB,DateCstar)
              Flow.date <- GetFlowdate(parameters$Fstar_BB_Fl_obs,temps.avgs,ActionunitsFl,DateBB)
              if(!is.na(Flow.date$Flowdate.julian)){
                Verais.date <- GetVeraisdate(parameters$Fstar_Fl_Ver_obs,temps.avgs,ActionunitsVer,Flow.date)
                
                row.to.store <- which(Results.phenology[,1]==Variety & Results.phenology[,2]==Site & Results.phenology[,4]==Year.target)
                Results.phenology[row.to.store,5] <- DateBB[[2]][1]
                Results.phenology[row.to.store,6] <- Flow.date[[2]][1]
                Results.phenology[row.to.store,7] <- Verais.date[[2]][1]
              }
              }
            }
            
          }
        }
        
      }
    }
    
  }
  
  return(Results.phenology)
}



## Combined function to be applied to very large matrices (from raster datasets)
## To be utilized with BEST data only where NAs are found in certain dates
#' @param x - vector with average daily temperatures   
#' @param variety - which variety to use parameters from
#' @return  vector with DateBB (date of budbreak),
#'          dateFlow (date of flowering), dateVER (date of veraison)
pheno.seq.fx.best <- function(x,variety){
  
  Tm1=variety$Tm1
  Topt=variety$Topt
  Tn2=variety$Tn2
  minn=variety$minn
  Tavg=x[3:length(x)]
  index.year <- which(period==Year.pheno)
  Year.previous <- period[index.year-1]
  Year.target <- period[index.year]
  
  if(sum(is.na(x))==0){
    ## start with Smoothed-UTAH
    if(x[2]>0){Hemisphere="NORTH"}else{Hemisphere="SOUTH"}
    
    
    if(Hemisphere=="NONE"|Hemisphere=="NORTH"){
      if(is.leapyear(Year.previous)){
        day1 <- 213
      } else {
        day1 <- 214
      }
    } 
    if(Hemisphere=="SOUTH"){
      if(is.leapyear(Year.previous)){
        day1 <- 213-182
      } else {
        day1 <- 214-182
      }
    }
    
    Chill.units <- as.vector(array(0,dim=c(length(Tavg),1)))
    for(h in day1:length(Tavg)){
      if(Tavg[h]<Tm1){ # note this is incorrect in original paper
        Chill.units[h]=1/(1+exp(-4*((Tavg[h]-Tm1)/(Topt-Tm1))))
        
      } else {
        
        if(Tavg[h]>=Tm1 & Tavg[h]<Topt){
          Chill.units[h]=1+(-0.5/((Tm1-Topt)^2))*(Tavg[h]-Topt)^2
          
        } else {
          
          if(Tavg[h]>Topt & Tavg[h]<Tn2){
            Chill.units[h]=1-(1-minn)*(((Tavg[h]-Topt)^2)/(2*((Tn2-Topt)^2)))
            
          } else {
            
            if(Tavg[h]>Tn2){
              Chill.units[h]=minn+((1-minn)/(1+exp(-4*((Tn2-Tavg[h])/(Tn2-Topt)))))
            }
          }
        }
      }
    }
    Chill.accum=cumsum(Chill.units)
    
    result.chilling <- data.frame(cbind(Chuill.units=Chill.units,Accummulated.chill=Chill.accum))
    #return(data.frame(cbind(Chuill.units=Chill.units,Accummulated.chill=Chill.accum)))
    Cstar=variety$Cstar_BB
    day1=approx(result.chilling[,2],seq(1,nrow(result.chilling)),xout=Cstar)$y
    
    DateCstars <- list(Cstardate=Tavg[round(day1,0)],Cstardate=day1)
    
    alpha.BB <- Alphafx(variety$Tmin,variety$Tmax,variety$Topt.1)
    ActionunitsBB <- WangEngelfx(variety$Tmin,variety$Tmax,variety$Topt.1,alpha.BB,Tavg)
    
    alpha.Fl <- Alphafx(variety$Tmin_BB_FL,variety$Tmax_BB_FL,variety$Topt_BB_FL)
    ActionunitsFl <- WangEngelfx(variety$Tmin_BB_FL,variety$Tmax_BB_FL,variety$Topt_BB_FL,alpha.Fl,Tavg)
    
    alpha.Ver <- Alphafx(variety$Tmin_FL_VER,variety$Tmax_FL_VER,variety$Topt_FL_VER)
    ActionunitsVer <- WangEngelfx(variety$Tmin_FL_VER,variety$Tmax_FL_VER,variety$Topt_FL_VER,alpha.Ver,Tavg)
    
    
    ## 4) get BB, Flowering and Veraison dates
    if(!is.na(DateCstars[[2]])){
      
      ##BB date
      Cstardate <- DateCstars[[2]]
      accum.Funits <- cumsum(ActionunitsBB[Cstardate:nrow(ActionunitsBB),1])
      day.BB <- which(accum.Funits>=variety$Fstar_BB)[1]+Cstardate
      DateBB=list(BBdate=Tavg[day.BB],BBdate.julian=day.BB)
      
      ##FL date
      if(!is.na(DateBB[[2]])){
        dateBB <- DateBB[[2]]
        accum.Funits <- cumsum(ActionunitsFl[dateBB:nrow(ActionunitsFl),1])
        day.Flow <- which(accum.Funits>=variety$Fstar_BB_Fl_obs)[1]+dateBB
        dateFlow <- list(Flowdate=Tavg[day.Flow],Flowdate.julian=day.Flow)
        
        ##VER date
        if(!is.na(dateFlow[[2]])){
          dateFlows <- dateFlow[[2]]
          accum.Funits <- cumsum(ActionunitsVer[dateFlows:nrow(ActionunitsVer),1])
          day.Verais <- which(accum.Funits>=variety$Fstar_Fl_Ver_obs)[1]+dateFlows
          dateVER <- list(Veraisdate=Tavg[day.Verais],Veraisdate.julian=day.Verais)
          
        } 
      }
    } else {
      DateBB <- list(NA,NA)
    }
    if(is.na(DateBB[[2]])){
      dateFlow <- list(NA,NA)
      dateVER <- list(NA,NA)
    }
    if(is.na(dateFlow[[2]])){
      dateVER <- list(NA,NA)
    }
  } else {
    DateBB <- list(NA,NA)
    dateFlow <- list(NA,NA)
    dateVER <- list(NA,NA)
  }
  
  
  return(c(DateBB[[2]],dateFlow[[2]],dateVER[[2]]))
}


