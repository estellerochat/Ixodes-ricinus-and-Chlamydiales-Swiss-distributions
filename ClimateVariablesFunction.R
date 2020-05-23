window <- function(x,n,f){ 
  nc <- ncol(x)
  out<- x[,1,drop=F]
  for(c in 1:(nc-n)){
    out<-cbind(out,apply(x[,c:(c+n-1),drop=F],1,f,na.rm=T))
  }
  out<-out[,-1,drop=F]
  return(out)
}

climvar<-function(tmean,tmax,tmin,prec){
  # tmean is a dataframe containing n rows for n occurrence points 
  # and one column for each monthly mean temperature of the period of interest
  # (e.g. if the period of interest is 3 months, tmean contains 3 columns with for each occurrence the monthly mean temperature of the 3 months preceding the sampling date)
  # tmin, tmax and prec have a similar strucutre but contain values of minimum temperature, maximum temperature and sum of precipitation, respectively.
  
  # 3 wettest, driest, warmest and coldest consecutive months if the time considered is at least 6 months
  if(ncol(tmean)>=6){
    wettest=apply(window(prec,3,sum),1,which.max)
    driest=apply(window(prec,3,sum),1,which.min)
    warmest=apply(window(tmean,3,mean),1,which.max)
    coldest=apply(window(tmean,3,mean),1,which.min)
  }
  
  if(ncol(tmean)>=6){
    c<-data.frame(matrix(0,nrow(tmean),31))
  }else{
    c<-data.frame(matrix(0,nrow(tmean),21))
  }
  
  # 1. meantmean
  c[,1] <- apply(tmean,1,mean,na.rm=T)
  colnames(c)[1]="meantmean"
  # 2. maxtmax
  c[,2] <- apply(tmax,1,max,na.rm=T)
  colnames(c)[2]="maxtmax"
  # 3. mintmin
  c[,3] <- apply(tmin,1,min,na.rm=T)
  colnames(c)[3]="mintmin"
  # 4. tRge
  c[,4] <- c[,2] - c[,3]
  colnames(c)[4]="tRge"
  # 5. tmeanMoRge
  c[,5] <- apply(tmax-tmin, 1, mean,na.rm=T)
  colnames(c)[5]="tmeanMoRge"
  # 6. isotherm
  c[,6] <- 100 * c[,5] / c[,4]
  colnames(c)[6]="isotherm"
  # 7. tseason
  c[,7] <- 100 * apply(tmean, 1, sd,na.rm=T)
  colnames(c)[7]="tseason"
  # 8. mintmean
  c[,8] <- apply(tmean,1,min,na.rm=T)
  colnames(c)[8]="mintmean"
  # 9. maxtmean
  c[,9] <- apply(tmean,1,max,na.rm=T)
  colnames(c)[9]="maxtmean"
  # 10. mintmax
  c[,10] <- apply(tmax,1,min,na.rm=T)
  colnames(c)[10]="mintmax"
  # 11. maxtmin
  c[,11] <- apply(tmin,1,max,na.rm=T)
  colnames(c)[11]="maxtmin"
  
  if(ncol(tmean)>=6){
    # 12. meantmean3cold
    c[,12] <- apply(window(tmean,3,mean),1,min,na.rm=T)
    colnames(c)[12]="meantmean3cold"
    # 13. meantmin3cold
    c[,13] <- window(tmin,3,mean)[cbind(1:length(coldest),coldest)]
    colnames(c)[13]="meantmin3cold"
    # 14. meantmax3cold
    c[,14] <- window(tmax,3,mean)[cbind(1:length(coldest),coldest)]
    colnames(c)[14]="meantmax3cold"
    # 15. meantmean3warm
    c[,15] <- apply(window(tmean,3,mean),1,max,na.rm=T)
    colnames(c)[15]="meantmean3warm"
    # 16. meantmin3warm
    c[,16] <- window(tmin,3,mean)[cbind(1:length(warmest),warmest)]
    colnames(c)[16]="meantmin3warm"
    # 17. meantmax3warm
    c[,17] <- window(tmax,3,mean)[cbind(1:length(warmest),warmest)]
    colnames(c)[17]="meantmax3warm"
  }
  
  # 18. sumprec
  c[,18] <- apply(prec, 1, sum,na.rm=T)
  colnames(c)[18]="sumprec"
  # 19. maxprec
  c[,19] <-  apply(prec, 1, max,na.rm=T)
  colnames(c)[19]="maxprec"
  # 20. minprec
  c[,20] <-  apply(prec, 1, min,na.rm=T)
  colnames(c)[20]="minprec"
  # 21. pseason
  c[,21] <- apply(prec+1, 1, cv)
  colnames(c)[21]="pseason"
  
  if(ncol(tmean)>=6){
    # 22. prec3wet
    c[,22] <- apply(window(prec,3,sum),1,max,na.rm=T)
    colnames(c)[22]="prec3wet"
    # 23. prec3dry
    c[,23] <- apply(window(prec,3,sum),1,min,na.rm=T)
    colnames(c)[23]="prec3dry"
    
    # 24. meantmean3wet
    c[,24] <- window(tmean,3,mean)[cbind(1:length(wettest),wettest)]
    colnames(c)[24]="meantmean3wet"
    # 25. meantmax3wet
    c[,25] <- window(tmax,3,mean)[cbind(1:length(wettest),wettest)]
    colnames(c)[25]="meantmax3wet"
    # 26. meantmin3wet
    c[,26] <- window(tmin,3,mean)[cbind(1:length(wettest),wettest)]
    colnames(c)[26]="meantmin3wet"
    # 27. meantmean3dry
    c[,27] <- window(tmean,3,mean)[cbind(1:length(driest),driest)]
    colnames(c)[27]="meantmean3dry"
    # 28. meantmax3dry
    c[,28] <- window(tmin,3,mean)[cbind(1:length(driest),driest)]
    colnames(c)[28]="meantmax3dry"
    # 29. meantmin3dry
    c[,29] <- window(tmax,3,mean)[cbind(1:length(driest),driest)]
    colnames(c)[29]="meantmin3dry"
    # 30. prec3warm
    c[,30] <- window(prec,3,sum)[cbind(1:length(warmest),warmest)]
    colnames(c)[30]="prec3warm"
    # 31. prec3cold
    c[,31] <- window(prec,3,sum)[cbind(1:length(coldest),coldest)]
    colnames(c)[31]="prec3cold"
  }
  
  if(ncol(tmean)<6){c<-c[,-(12:17)]}
  
  return(c)
}

RHvar<-function(RHmean,RHq050,RHq025,RHq075){
  # RHmean is a dataframe containing n rows for n occurrence points 
  # and one column for each monthly mean relative humidity of the period of interest
  # (e.g. if the period of interest is 3 months, RHmean contains 3 columns with for each occurrence the monthly relative humidity values of the 3 months preceding the sampling date)
  # RHq050, RHq025, RHq075 have a similar strucutre but contain values of monthly quantile 0.5, 0.25 and 0.75 of relative humidity values
  
  rh=data.frame(matrix(0,nrow(RHmean),9))
  
  # 1. meanRHmean
  rh[,1]=apply(RHmean,1,mean,na.rm=T)
  colnames(rh)[1]="meanRHmean"
  
  # 2. meanRHq050
  rh[,2]=apply(RHq050,1,mean,na.rm=T)
  colnames(rh)[2]="meanRHq050"
  
  # 3. minRHmean
  rh[,3]=apply(RHmean,1,min,na.rm=T)
  colnames(rh)[3]="minRHmean"
  
  # 4. maxRHmean
  rh[,4]=apply(RHmean,1,max,na.rm=T)
  colnames(rh)[4]="maxRHmean"
  
  # 5. minRHq025
  rh[,5]=apply(RHq025,1,min,na.rm=T)
  colnames(rh)[5]="minRHq025"
  
  # 6. minRHq075
  rh[,6]=apply(RHq075,1,min,na.rm=T)
  colnames(rh)[6]="minRHq075"
  
  # 7. maxRHq075
  rh[,7]=apply(RHq075,1,max,na.rm=T)
  colnames(rh)[7]="maxRHq075"
  
  # 8. RHrge
  rh[,8]=rh[,7]-rh[,5]
  colnames(rh)[8]="RHrge"
  
  # 9. RHMoRge
  rh[,9]=apply(RHq075-RHq025,1,mean,na.rm=T)
  colnames(rh)[9]="RHMoRge"
  
  return(rh)
}

