
library(raster)
library(abind)
library(mgcv)
library(data.table)
library(sf)
library(viridis)
library(rasterVis)

path<-"C:/Users/God/Documents/mosquitos/data"

lf<-list.files(file.path(path,"daymet"),full=TRUE,pattern="tmax")
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

tiles<-paste0("_",gsub(".nc","",unique(sapply(strsplit(list.files(file.path(path,"daymet"),full=TRUE,pattern="tmax"),"_"),"[",3))))
#tiles<-tiles[1]

for(k in seq_along(tiles)){

tile<-tiles[k]
g<-grep(tile,names(lf))
lr<-lf[g]
lr<-lapply(lr,aggregate,10)
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
an<-lapply(seq_along(ts),function(j){
  i<-ts[[j]]
  x<-data.frame(tmax=i)
  x$date<-as.Date(names(i),format="X%Y.%m.%d")
  x$jul<-as.integer(format(x$date,"%j"))
  if(all(is.nan(x[,"tmax"]))){
    p<-empty
  }else{
    #plot(x[,"jul"],x[,"tmax"],ylim=c(-25,33))
    m<-gam(tmax~s(jul,bs="cc"),data=x)
    p<-predict(m,data.frame(jul=1:365))
    #lines(1:365,p,lwd=5,col=alpha("red",0.5))
    #x
  }
  x$normal<-rep(p,nrow(x)/365)
  x$anomaly<-x$tmax-x$normal
  cat("\r",paste(k,length(tiles)," - ",j,length(ts),sep=" / "))
  x
  #matrix(c(x$normal,x$anomaly),ncol=2)
})
t2<-Sys.time()
t2-t1

aa<-a
for(i in 1:nrow(ii)){
  aa[ii[i,1],ii[i,2],]<-an[[i]]$anomaly
  #print(i)
}

aa<-lapply(split(1:dim(aa)[3], ceiling(1:dim(aa)[3]/365)),function(i){
  aa[,,i]
})

temp<-lapply(seq_along(lr),function(i){
  r<-setValues(lr[[i]],aa[[i]])
  names(r)<-names(lr[[i]])
  r
})
names(temp)<-names(lr)


#levelplot(temp[[1]][[ceiling(seq(1,365,length.out=20))]],zlim=c(-8,0),col.regions=viridis(200),cuts=199,margin=FALSE) +
#  layer(sp.polygons(as(st_transform(st_as_sf(Q),crs=st_crs(lr[[1]])),"Spatial"),col=gray(0,0.25)))

for(kk in seq_along(temp)){
  writeRaster(temp[[kk]],filename=file.path(path,"daymet",paste0("anom_",names(temp)[kk])),format="CDF")
}

}



r<-stack(paste0(file.path(path,"daymet",paste0("anom_",names(temp)[1])),".nc"))

#a<-rbindlist(an)
#a[,tile:=tile]
#a[,year:=as.integer(stri_sub(date,1,4))]
#a<-split(a,paste(a$year,a$tile))

r<-lf[[paste0("2013",tile)]][[1]]
plot(aggregate(r,2),add=FALSE,col=plasma(200))
plot(r,add=TRUE)
plot(as(st_transform(st_as_sf(Q),crs=st_crs(lr[[1]])),"Spatial"),border=gray(0,0.25),add=TRUE)

r[30000:30100]<-NA
plot(r)
plot(aggregate(r,2),add=FALSE,col=plasma(200))


r<-stack(paste0(file.path(path,"daymet",paste0("anom_",names(temp)[1])),".nc"))



tiles<-paste0("_",gsub(".nc","",unique(sapply(strsplit(list.files(file.path(path,"daymet"),full=TRUE,pattern="tmax"),"_"),"[",3))))
lr<-lapply(2003:2016,function(i){
  stack(paste0(file.path(path,"daymet",paste0("anom_",i,tiles[4])),".nc"))
})
par(mfrow=n2mfrow(length(lr)))
lapply(lr,function(i){
  plot(i[[1]])
})


