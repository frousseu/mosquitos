### Packages #########################################
library(sp)
library(rgdal)
library(scales)
library(INLA)
library(rgeos)
library(FRutils)
library(readxl)
library(raster) 
library(rasterVis)
library(data.table)
library(alphahull)
library(concaveman)
library(mapview)
library(fasterize)
library(sf)
library(velox) 
library(viridis)
library(foreach)
library(doParallel)
library(rmapshaper)
library(ncdf4)
library(daymetr)
library(gdalUtils)
library(rgdal)
library(exactextractr)
library(foreach)
library(doParallel)
library(DHARMa)
library(corrplot)
library(abind)
library(mgcv)
library(dplyr)
library(abind)
library(rasterVis)
library(terra)
library(doSNOW)
library(magick)
library(splines)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(car)


path<-"C:/Users/God/Documents/mosquitos/data"
setwd(path)


### Daymet download #########################################

# The Daymet data needs to be downloaded first with the following code

### code to download daymet data
###rasterOptions(chunksize=1e+09,maxmemory=5e+10)
download_daymet_tiles(location = c(46.75,-77.5,45,-70),
##  location = st_bbox(st_transform(st_as_sf(ds),crs=4326))[c(4,1,2,3)],
##  location=c(ymax=45.9022729397129,xmin=-74.4357158151628,ymin=45.1889865783487,xmax=-72.8946830482002), # topleft and bottom right coordinates of ds object
###  tiles=c("12472", "12473", "12474", "12475"),
  start = 2003,
  end = 2016,
  param = c("tmin","tmax","prcp"),
  path = file.path(path,"daymet"))

# faster to stack rasters, summarize them to weekly values and then merge them, not merge and then summarize

### Build Tmean #########################################

# Computes Tmean from Tmax and Tmin

lf<-list.files(file.path(path,"daymet"),full=TRUE,pattern="tmax")
lf<-lf[order(lf)]
ids<-do.call("rbind",strsplit(lf,"_"))
o<-order(ids[,3],ids[,2])
ids<-ids[o,]
lf<-lf[o]
lf<-split(lf,paste(ids[,2],gsub(".nc","",ids[,3]),sep="_"))
lf<-lapply(seq_along(lf),function(i){
  tmax<-rast(lf[[i]])
  tmin<-rast(gsub("tmax_","tmin_",lf[[i]]))
  tmean<-((tmax-tmin)/2)+tmin
  #writeRaster(tmean,
  writeCDF(tmean,filename=gsub("tmax_","tmean_",lf[[i]]),overwrite=TRUE)
  print(i)
}) 


### Anomalies ################################################

lf<-list.files(file.path(path,"daymet"),full=TRUE,pattern="tmean")
lf<-lf[order(lf)]
ids<-do.call("rbind",strsplit(lf,"_"))
o<-order(ids[,3],ids[,2])
ids<-ids[o,]
lf<-lf[o]
lf<-split(lf,paste(ids[,2],gsub(".nc","",ids[,3]),sep="_"))
lf<-lapply(lf,function(i){
  stack(i)
})

#bbb<-lapply(lf,function(i){
#  st_as_sfc(st_bbox(i),crs=st_crs(i))
#})
#bbb<-Reduce(st_union,bbb)
#plot(st_geometry(bbb))
#plot(st_geometry(st_transform(st_as_sf(Q),crs=st_crs(lr[[1]]))),add=TRUE)
#for(i in seq_along(lf)){
#  bb<-st_as_sfc(st_bbox(lf[[i]]),crs=st_crs(lf[[i]]))
#  plot(bb,add=TRUE)
#  text(st_coordinates(st_centroid(bb))[,1],st_coordinates(st_centroid(bb))[,2],label=sapply(strsplit(names(lf)[i],"_"),"[",2#))
#}

tiles<-paste0("_",gsub(".nc","",unique(sapply(strsplit(list.files(file.path(path,"daymet"),full=TRUE,pattern="tmean"),"_"),"[",3))))
#tiles<-tiles[1]

# reste 4 and 6+

foreach(k=seq_along(tiles)[1:13],.packages=c("abind","mgcv","raster")) %do% {
  tile<-tiles[k]
  g<-grep(tile,names(lf))
  lr<-lf[g]
  lr<-lapply(lr,aggregate,2)
  a<-abind(lapply(lr,function(i){
    a<-as.array(i)
    dimnames(a)[[3]]<-names(i)
    a
  }))
  ii<-expand.grid(i=seq(dim(a)[1]),j=seq(dim(a)[2]))
  ts<-lapply(1:nrow(ii),function(x){
    a[ii[x,"i"] ,ii[x,"j"],]
  })
  
  empty<-rep(NA,365)
  t1<-Sys.time()
  #cl <- makeCluster(detectCores ()-6)
  #registerDoSNOW(cl)
  an<-foreach(j=seq_along(ts),.packages=c("mgcv")) %do% {
    #an<-lapply(seq_along(ts),function(j){
    i<-ts[[j]]
    x<-data.frame(tmean=i)
    x$date<-as.Date(substr(names(i),1,11),format="X%Y.%m.%d")
    x$jul<-as.integer(format(x$date,"%j"))
    if(all(is.nan(x[,"tmean"]))){
      p<-empty
    }else{
      #plot(x[,"jul"],x[,"tmean"],ylim=c(-25,33))
      m<-gam(tmean~s(jul,bs="cc"),data=x)
      p<-predict(m,data.frame(jul=1:365))
      #lines(1:365,p,lwd=5,col=alpha("red",0.5))
      #x
    }
    #closeAllConnections()
    x$normal<-rep(p,nrow(x)/365)
    x$anomaly<-x$tmean-x$normal
    cat("\r",paste(k,length(tiles)," - ",j,length(ts),sep=" / "))
    #plot(x$anomaly)
    x$anomaly
    
    #matrix(c(x$normal,x$anomaly),ncol=2)
  }#)
  t2<-Sys.time()
  t2-t1
  rm(ts);gc()
  
  for(i in 1:nrow(ii)){
    a[ii[i,1],ii[i,2],]<-an[[i]]
  }
  rm(an);gc()
  
  a<-lapply(split(1:dim(a)[3], ceiling(1:dim(a)[3]/365)),function(i){
    a[,,i]
  })
  
  temp<-lapply(seq_along(lr),function(i){
    r<-setValues(lr[[i]],a[[i]])
    names(r)<-names(lr[[i]])
    r
  })
  names(temp)<-names(lr)
  
  
  #levelplot(temp[[1]][[ceiling(seq(1,365,length.out=20))]],zlim=c(-8,0),col.regions=viridis(200),cuts=199,margin=FALSE) +
  #  layer(sp.polygons(as(st_transform(st_as_sf(Q),crs=st_crs(lr[[1]])),"Spatial"),col=gray(0,0.25)))
  
  for(kk in seq_along(temp)){
    writeRaster(temp[[kk]],filename=file.path(path,"daymet",paste0("anom_",names(temp)[kk])),format="CDF",overwrite=TRUE)
  }
  rm(temp,a);gc()
  cat("\r",paste(k,length(tiles)))
}

