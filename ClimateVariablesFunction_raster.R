windowStack <- function(x,n,f)  {
  nc <- length(names(x))
  out<- x[[1]]
  for(c in 1:(nc-n)){
    out<-stack(out,calc(x[[c:(c+n-1)]],f))
  }
  out<-out[[-1]]
  return(out)
}

climvarRast<-function(tmean,tmax,tmin,prec){
  
  beginCluster()
  
  # 3 wettest, driest, warmest and coldest consecutive months if time period considered is at least 6 months
  if(length(names(tmean))>=6){
    print("windows")
    precWindow<-windowStack(prec,3,sum)
    tmeanWindow<-windowStack(tmean,3,mean)
    
    print("which.max")
    wettest=which.max(precWindow)
    driest=which.min(precWindow)
    warmest=which.max(tmeanWindow)
    coldest=which.min(tmeanWindow)
  }
  
  c<-raster()
  
  # 1. meantmean
  c <- stack(c,mean(tmean,na.rm=T))
  names(c)[1]="meantmean"
  # 2. maxtmax
  c <- stack(c,max(tmax,na.rm=T))
  names(c)[2]="maxtmax"
  # 3. mintmin
  c <- stack(c,min(tmin,na.rm=T))
  names(c)[3]="mintmin"
  # 4. tRge
  c <- stack(c,c[[2]] - c[[3]])
  names(c)[4]="tRge"
  # 5. tmeanMoRge
  c <- stack(c,mean(tmax-tmin,na.rm=T))
  names(c)[5]="tmeanMoRge"
  # 6. isotherm
  c <- stack(c, 100 * c[[5]] / c[[4]])
  names(c)[6]="isotherm"
  # 7. tseason
  c <- stack(c, 100 * calc(tmean,sd,na.rm=T))
  names(c)[7]="tseason"
  # 8. mintmean
  c <- stack(c, min(tmean,na.rm=T))
  names(c)[8]="mintmean"
  # 9. maxtmean
  c <- stack(c,max(tmean,na.rm=T))
  names(c)[9]="maxtmean"
  # 10. mintmax
  c <- stack(c,min(tmax,na.rm=T))
  names(c)[10]="mintmax"
  # 11. maxtmin
  c <- stack(c,max(tmin,na.rm=T))
  names(c)[11]="maxtmin"
  
  if(length(names(tmean))>=6){
    # 12. meantmean3cold
    c <- stack(c,min(tmeanWindow,na.rm=T))
    names(c)[12]="meantmean3cold"
    # 13. meantmin3cold
    c <- stack(c,stackSelect(windowStack(tmin,3,mean),coldest))
    names(c)[13]="meantmin3cold"
    # 14. meantmax3cold
    c <- stack(c,stackSelect(windowStack(tmax,3,mean),coldest))
    names(c)[14]="meantmax3cold"
    # 15. meantmean3warm
    c <- stack(c,max(tmeanWindow,na.rm=T))
    names(c)[15]="meantmean3warm"
    # 16. meantmin3warm
    c <- stack(c,stackSelect(windowStack(tmin,3,mean),warmest))
    names(c)[16]="meantmin3warm"
    # 17. meantmax3warm
    c <- stack(c,stackSelect(windowStack(tmax,3,mean),warmest))
    names(c)[17]="meantmax3warm"
  }
  
  nc=length(names(c))
  
  # 18. sumprec
  c <- stack(c,sum(prec,na.rm=T))
  names(c)[nc+1]="sumprec"
  # 19. maxprec
  c <- stack(c,max(prec,na.rm=T))
  names(c)[nc+2]="maxprec"
  # 20. minprec
  c <-  stack(c,min(prec,na.rm=T))
  names(c)[nc+3]="minprec"
  # 21. pseason
  c <- stack(c,calc(prec+1,fun=cv,na.rm=T))
  names(c)[nc+4]="pseason"
  
  if(length(names(prec))>=6){
    # 22. prec3wet
    c <- stack(c,max(precWindow,na.rm=T))
    names(c)[nc+5]="prec3wet"
    # 23. prec3dry
    c <- stack(c,min(precWindow,na.rm=T))
    names(c)[nc+6]="prec3dry"
    # 24. meantmean3wet
    c <- stack(c,stackSelect(tmeanWindow,wettest))
    names(c)[nc+7]="meantmean3wet"
    # 25. meantmax3wet
    c <- stack(c,stackSelect(windowStack(tmax,3,mean),wettest))
    names(c)[nc+8]="meantmax3wet"
    # 26. meantmin3wet
    c <- stack(c,stackSelect(windowStack(tmin,3,mean),wettest))
    names(c)[nc+9]="meantmin3wet"
    # 27. meantmean3dry
    c <- stack(c,stackSelect(windowStack(tmean,3,mean),driest))
    names(c)[nc+10]="meantmean3dry"
    # 28. meantmax3dry
    c <- stack(c,stackSelect(windowStack(tmin,3,mean),driest))
    names(c)[nc+11]="meantmax3dry"
    # 29. meantmin3dry
    c <- stack(c,stackSelect(windowStack(tmax,3,mean),driest))
    names(c)[nc+12]="meantmin3dry"
    # 30. prec3warm
    c <- stack(c,stackSelect(precWindow,warmest))
    names(c)[nc+13]="prec3warm"
    # 31. prec3cold
    c <- stack(c,stackSelect(precWindow,coldest))
    names(c)[nc+14]="prec3cold"
  }
  
  endCluster()
  
  return(c)
}
