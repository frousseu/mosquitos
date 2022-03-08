
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

#cat("\014")
options(device = "X11")
grDevices::windows.options(record=TRUE)

Sys.setlocale("LC_ALL","English")

# first set working directory
# all files should be in this folder
path<-"C:/Users/God/Documents/mosquitos/data"
setwd(path)
load("mosquitos.RData")

### Daymet download #########################################

# The Daymet data needs to be downloaded first with the following code

### code to download daymet data
###rasterOptions(chunksize=1e+09,maxmemory=5e+10)
#download_daymet_tiles(location = c(46.75,-77.5,45,-70),
##  location = st_bbox(st_transform(st_as_sf(ds),crs=4326))[c(4,1,2,3)],
##  location=c(ymax=45.9022729397129,xmin=-74.4357158151628,ymin=45.1889865783487,xmax=-72.8946830482002), # topleft and bottom right coordinates of ds object
###  tiles=c("12472", "12473", "12474", "12475"),
#  start = 2003,
#  end = 2016,
#  param = c("tmin","tmax","prcp"),
#  path = file.path(path,"daymet"))

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



### Functions #########################################

### useful functions for plotting and predictions
source("https://raw.githubusercontent.com/frousseu/FRutils/master/R/colo.scale.R")
source("https://raw.githubusercontent.com/frousseu/UdeS/master/GStecher/newdata.R")


### function for turning mesh to sp polygons
inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
      calling inla.mesh2sp()."))
  }
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],1:2,drop = FALSE])),ID = x)
      }
    ),proj4string = crs),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  list(triangles = triangles, vertices = vertices)
  }

### projection that will be used and a contour map of Quebec
prj<-"+proj=utm +zone=18 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"
can<-raster::getData("GADM", country = "CAN", level = 2)
que<-can[can$NAME_1=="Québec",]
Q<-as(ms_simplify(st_as_sf(spTransform(que,CRS(prj))),keep=0.03),"Spatial")


### Read mosquito data #######################################

# Do we have trap nights when counts were 0 ???

### GDG #############################################
gdg<-as.data.frame(read_excel("GDG.xls",sheet="Total"))
ngdg<-names(gdg)
gdg$id<-gdg$SiteP
gdg$lon<-gdg$Longitude
gdg$lat<-gdg$Latitude
gdg$date<-as.character(gdg$Collection_date) # make sure the date is OK cause it is read as POSIX
gdg$nights<-gdg$NbNuitOp_Calc
gdg$code<-gdg$Species
gdg$species<-gdg$Taxa_Name
gdg$count<-gdg$Total_MS
gdg$type<-gdg$Typ
gdg$method<-gdg$MethP

gdg<-gdg[,setdiff(names(gdg),ngdg)]
gdg$db<-"gdg"

table(gdg$count>0)

            
### INSPQ ###########################################

# A single trap id can have multiple method...
# What is H and E dates?
# What are the different coordinates?
# Is this the total count taking into account the subsampling?

### check count values in column Pools_Entomo__NbMS in every sheet
### All sheets seem to have values over 50, except the two last
### I assume we have to sum every line
sheets<-excel_sheets("INSPQ.xlsx")
par(mfrow=n2mfrow(length(sheets)),mar=c(5,5,1,1),oma=c(0,0,3,1))
lapply(sheets,function(i){
  x<-as.data.frame(read_excel("INSPQ.xlsx",sheet=i))
  h<-hist(x$Pools_Entomo__NbMS,main=i,breaks=20,plot=FALSE)
  b<-barplot(ifelse(h$counts<1,NA,h$counts),log="y",yaxt="n",ylab="Freq",xlab="",border="grey70",plot=FALSE)[,1]
  barplot(ifelse(h$counts<1,NA,h$counts),log="y",yaxt="n",ylab="Freq",xlab="",border="grey70")
  axis(1,at=b,labels=paste0("[",h$breaks[-length(h$breaks)]," - ",h$breaks[-1],"]"),las=2,cex.axis=0.75)
  axis(2,las=2)
  mtext(side=3,line=-2,text=paste("Sheet",i))
})
mtext(side=3,line=1,text="Column Pools_Entomo__NbMS",outer=TRUE)


### Read data from inspq first sheet
inspq<-as.data.frame(read_excel("INSPQ.xlsx",sheet="Sheet1"))
ninspq<-names(inspq)
inspq$id<-inspq$Pools_Entomo__Site
inspq$lon<-inspq$Stations_Entomo__Longdec
inspq$lat<-inspq$Stations_Entomo__Latdec
inspq$method<-inspq$Pools_Entomo__MethP
inspq$species<-inspq$Pools_Entomo__TaxaName
inspq$code<-inspq[,7]
inspq$count<-inspq$Pools_Entomo__NbMS
inspq$type<-inspq$Pools_Entomo__Typ
inspq$date<-as.character(inspq$Echantillons_Entomo__ColDate) # make sure the date is OK cause it is read as POSIX
inspq<-inspq[!is.na(inspq$date),] # drop weird lines without dates
inspq$nights<-inspq$Information_Supp_Entomo__NB_Nuits_Operation
inspq<-inspq[,setdiff(names(inspq),ninspq)]
inspq$db<-"inspq"

table(inspq$count>0)


### Pool GDG and INSPQ ############################

### all names are common to both database
setdiff(names(gdg),names(inspq))
setdiff(names(inspq),names(gdg))

### bind data
d<-as.data.table(rbind(gdg,inspq[,names(gdg)]))

### some cleaning
d<-d[!is.na(d$species),] # remove NA species
d<-d[d$method%in%c("LT"),] # keep only LTs
table(d$type) # no male listed
d<-d[,type:=NULL] # drop type column
d$species<-gsub("à","a",gsub("é","e",gsub("/| |-","_",gsub("\\.","",d$species)))) # remove accents and special characters

### check number of lines to sum (or duplicate lines based on id, lon, lat, date, species)
d[,nlines:=.N,by=.(id,lon,lat,date,code,species,db)] # counts the nb of lines per by
d[nlines>1,][order(-nlines,id,date)] # shows the highest nb off lines
d[,nlines:=NULL] # removes le the nb of lines

### sum values and keep varying nights
# some traps have more than one nb of nights of effort
d<-d[,.(count=sum(count)),by=.(id,lon,lat,date,nights,code,species,db)] # sums values per by
nrow(d)

### show varying night (3 cases)
d[,nlines:=.N,by=.(id,lon,lat,date,code,species,db)]
d[nlines>1,][order(-nlines,id,date)]
nrow(d)

### read in species / code correspondance
taxa<-as.data.frame(read_excel("GDG.xls",sheet="Taxa"))
taxa$species<-gsub("/| |-","_",gsub("\\.","",taxa$TaxaName))
taxa$code<-taxa$TaxaCode

### show code species that don't have matches in taxa names and drop them
# some species seem to have wrong codes
# those species not used in this analysis
table(paste(d$code,d$species)[!paste(d$code,d$species)%in%paste(taxa$code,taxa$species)])
nrow(d)
d<-d[paste(d$code,d$species)%in%paste(taxa$code,taxa$species),]
nrow(d)

### check if several lines per species
# seems to be a problem with different night values, will fix later
d[,nlines:=.N,by=.(id,lon,lat,date,species,db)]
d[nlines>1,]

### drop QUW code
d<-d[d$code!="QUW",]
d[,nlines:=NULL]

### verify if single line per species, per trap, per date
table(d[,.(nlines=.N),by=.(id,lon,lat,date,code,species,db)]$nlines)
table(d[,.(nlines=.N),by=.(id,lon,lat,date,species,db)]$nlines)

### paste code and species to ID species and remove old columns
d[,species:=paste(code,species,sep="_"),]
d[,code:=NULL]

### cast to get species columns
# there should not be any waring here, otherwise it means there are more than one values
d<-dcast(d,...~species,value.var="count",fill=0) 

### check if there duplicate id/dates 
d[,nlines:=.N,by=.(id,date)]
d[nlines>1,]
d[,nlines:=NULL]

### order by night and drop the longest number of nights (don't know what to do with this...)
nrow(d)
setorder(d,id,date,nights)
d<-d[!duplicated(d[,c("id","lon","lat","date")]),]
nrow(d)

### are there any duplicates leftovers for id,lon,lat,date?
nrow(d)
setorder(d,id,date)
table(duplicated(d[,c("id","date")]))

### order species per reverse total abundance
notspecies<-c("id","lon","lat","date","nights","db")
species<-setdiff(names(d),notspecies)
total<-colSums(d[,..species])
species<-species[rev(order(total))]
d<-d[,c(notspecies,species),with=FALSE]
rev(sort(total))

### random checks
d[d$date=="2011-07-12" & d$id=="TER001",]


### add info
d$jul<-as.integer(format(as.Date(d$date),"%j"))
#d$jul<-d$jul/100
#d$julsquare<-d$jul^2
d$year<-substr(d$date,1,4)
d$week<-format(as.Date(d$date),"%Y_W%U")
#d$year_week<-paste(d$year,d$week,sep="_")

setorder(d,id,week)
unique(d[,c("id","lon","lat","year")])[,lapply(.SD,length),by=year,.SDcols="id"][order(year),]
d[,lapply(.SD,sum),by=year,.SDcols=species][order(year),1:7]

l<-split(d,d$year)
ds<-d
coordinates(ds)<-~lon+lat
ds$longitude<-d$lon
ds$latitude<-d$lat
proj4string(ds)<-"+init=epsg:4326"
ds<-spTransform(ds,CRS(prj))


### Region ##############################################

reg<-st_read("C:/Users/God/Documents/mosquitos/data/SHP/SHP",layer="regio_s")
reg$region[grep("Outaouais|Laurentides|Lanaudière|Abitibi|Côte|Centre|Capitale|Nord|Mauricie|Saguenay",reg$RES_NM_REG)]<-"Nord"
reg$region[grep("Montréal",reg$RES_NM_REG)]<-"Montréal"
reg$region[grep("Laval",reg$RES_NM_REG)]<-"Laval"
reg$region[grep("Estrie|Montérégie|Centre|Chaudière|Gaspésie|Bas",reg$RES_NM_REG)]<-"Sud"
reg <- reg %>% group_by(region) %>% summarise(geometry = st_union(geometry))
plot(st_geometry(reg))
reg<-st_transform(reg,st_crs(ds))

ds$region<-st_join(st_as_sf(ds),reg,join=st_intersects)$region

plot(ds,col=factor(ds$region))

### Show trap data ######################################

#### By year ###########################################
l<-split(ds[!ds$db%in%"map",],ds$year[!ds$db%in%"map"])
par(mfrow=n2mfrow(length(l)),mar=c(0.25,0.25,0.25,0.25))
lapply(l,function(i){
  x<-i[!duplicated(i$longitude),]
  print(nrow(x))
  plot(x,col="white")
  plot(Q,add=TRUE,col="grey90",border="white")
  plot(x,add=TRUE,pch=16,cex=1,col="forestgreen")
  plot(st_geometry(reg),add=TRUE,border=gray(0,0.15),lwd=3)
  mtext(i$year[1],side=3,adj=c(0,0),line=-1.25)
})
par(mfrow=c(1,1))

#### By weeks ###########################################
l<-split(ds[!ds$db%in%"map",],ds$week[!ds$db%in%"map"])
names(l)<-sapply(l,function(i){i$week[1]})
ee<-expand.grid(year=sort(unique(substr(names(l),1,4))),week=sort(unique(substr(names(l),6,8))))
m<-match(apply(ee,1,function(i){paste(i[1],i[2],sep="_")}),names(l))
ee$nbtraps<-sapply(m,function(i){if(is.na(i)){0}else{nrow(l[[i]])}})
ee<-ee[order(ee$year,ee$week),]


par(mfrow=c(length(unique(ds$year)),nrow(ee)/length(unique(ds$year))),mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1))
lapply(1:nrow(ee),function(j){
  m<-match(paste(ee[j,"year"],ee[j,"week"],sep="_"),names(l))
  if(is.na(m)){
    plot(1,1,type="n",axes=FALSE)
  }else{
    i<-l[[m]]
    x<-i[!duplicated(i$longitude) & !duplicated(i$longitude),]
    print(nrow(x))
    plot(x,col="white",xlim=bbox(ds[ds$db!="map",])[1,],ylim=bbox(ds[ds$db!="map",])[2,])
    plot(Q,add=TRUE,col="grey90",border="white")
    plot(x,add=TRUE,pch=16,cex=0.5,col=gray(0,0.75))
    mtext(i$week[1],side=3,adj=c(0,0),line=-0.5,cex=0.4)
  }
})
par(mfrow=c(1,1))

#### Plot traps/weeks #####################################
nbtraps<-split(ee,ee$year)
par(mfrow=c(length(nbtraps),1),mar=c(0.5,3,0,0),oma=c(3,0,0,0))
lapply(nbtraps,function(i){
  b<-barplot(i$nbtraps,names.arg=if(identical(i,nbtraps[[length(nbtraps)]])){i$week}else{NULL},ylim=c(0,max(ee$nbtraps)),border=NA,col="forestgreen",xpd=FALSE,yaxt="n")
  text(b[,1],rep(-15,nrow(i)),i$nbtraps,cex=0.75,xpd=TRUE)
  mtext(side=3,line=-1.5,text=i$year[1],adj=0.005)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = gray(0,0.05),border=NA)
  axis(2,las=TRUE)
  #grid()
})
par(mfrow=c(1,1))


#### Keep what's in the zone #######################################
plot(ds)
plot(Q,add=TRUE)
#l<-locator()
l<-list(x = c(537.794779946664, 537.794779946664, 558.632058968996, 
           596.566592573754, 638.77543982412, 653.201248378042, 666.558478520563, 
           671.901370577571, 665.489900109161, 647.324067115333, 629.158234121505, 
           579.469337991328, 558.097769763295), y = c(5063.9912118604, 5036.20817316396, 
                                                      4999.8765071763, 4999.8765071763, 5003.08224241051, 5013.23373731882, 
                                                      5029.26241348985, 5049.03111410078, 5067.19694709461, 5083.22562326563, 
                                                      5088.56851532264, 5081.62275564853, 5079.48559882573))

h<-gConvexHull(SpatialPoints(cbind(l$x,l$y),proj4string=CRS(proj4string(ds))))
plot(h,add=TRUE)
o<-over(ds,h)
d<-d[!is.na(o),]
ds<-ds[!is.na(o),]


#### Remove weeks with few traps #################################
# currently ignored, better to do this before running the models when subsetting
d[,nbtrapweek:=.N,by=.(week)]
k<-d$nbtrapweek>=0
d<-d[k,]
ds<-ds[k,]


### Show trap weeks ###################################

d$pch<-15
x<-dcast(unique(d[,c("id","week","pch")]),id~week,value.var="pch")
x<-melt(x,id=c("id","week"))
x$variable<-factor(x$variable,levels=unique(x$variable))
x$id<-factor(x$id,levels=unique(x$id))
par(mar=c(5,4,1,1))
plot(as.integer(x$variable),as.integer(x$id),pch=x$value,xaxt="n",yaxt="n",xlab="weeks",ylab="traps")
axis(1,at=1:nlevels(x$variable),labels=levels(x$variable),las=2,cex.axis=0.5,srt=45)
axis(2,at=1:nlevels(x$id),labels=levels(x$id),las=2,cex.axis=0.35)


### Build mapping zone #############################

mappingzone<-concaveman(coordinates(ds[ds$db!="map",]),2)
mappingzone<-gBuffer(spPolygons(mappingzone,crs=CRS(proj4string(ds))),width=10)
plot(mappingzone)
plot(Q,add=TRUE,border="grey80")
plot(ds,add=TRUE,pch=1)

### Build prediction grid ###################################

# add data for prediction maps to ds with NAs for numbers of mosquitos
# this allows to make the lcc/daymet extractions in a single step
# and INLA does predictions for NA response values

pgrid<-raster(ext=extent(mappingzone),res=c(1,1),crs=CRS(proj4string(mappingzone)))
pgrid<-setValues(pgrid,1)
g<-xyFromCell(pgrid,1:ncell(pgrid),spatial=TRUE)
plot(pgrid)
plot(g,add=TRUE,cex=0.3)

ee<-expand.grid(week=sort(unique(ds$week)))
ee$id<-"map"
ee$year<-as.integer(substr(ee$week,1,4))
# need to specify day of week when converting week to date see https://stackoverflow.com/questions/9380435/how-to-parse-year-week-number-in-r
ee$date<-as.Date(paste0(ee$week,3),format="%Y_W%U%u")
ee$jul<-as.integer(format(ee$date,"%j"))
ee$date<-as.character(ee$date)
ee$nights<-1
ee$db<-"map"
eesp<-as.data.frame(matrix(rep(NA,nrow(ee)*length(species)),ncol=length(species)))
names(eesp)<-species
ee<-cbind(ee,eesp)
coo<-as.data.frame(coordinates(g))
names(coo)<-c("longitude","latitude")
coo$lon<-coo$longitude
coo$lat<-coo$latitude

info<-ee[rep(1:nrow(ee),each=nrow(coo)),]
coords<-coo[rep(1:nrow(coo),times=nrow(ee)),]

map<-cbind(info,coords)
coordinates(map)<-~lon+lat
proj4string(map)<-CRS(proj4string(ds))

# subset predictions only for a given week (reduces extractions times/memory use)
#map<-map[substr(map$week,6,8)%in%"W32",]
map<-map[map$week%in%"2003_W32",]

# id region in case of spatial validation
map$region<-st_join(st_as_sf(map),reg,join=st_intersects)$region
#plot(map,col=factor(map$region))

setdiff(names(ds),names(map))
setdiff(names(map),names(ds))

map<-map[,names(ds)]


ds<-rbind(ds,map)

# remove what is not in the mapping zone
o<-over(ds,mappingzone)
ds<-ds[!is.na(o),]


### Weather lags #####################################

weathervars<-c("tmean","prcp","anom")
#weathervars<-c("tmean","tmax","tmin","prcp","anom")
#weathervars<-c("anom")
#weathervars<-c("tmax","anom")
#weathervars<-c("anom")

for(v in seq_along(weathervars)){

  l<-list.files(file.path(path,"daymet"),full=TRUE,pattern=paste0("tmax","_2015_"))
  g<-lapply(l,function(i){
    r<-stack(i)[[1]]
    r<-setValues(r,ifelse(is.na(values(r)),NA,as.integer(gsub(".nc","",tail(strsplit(i,"_")[[1]],1)))))
    r
  })
  bb<-lapply(g,function(i){
    st_as_sfc(st_bbox(i))
  })
  g<-Reduce(merge,g)

  #plot(g[[1]])
  #invisible(lapply(bb,plot,add=TRUE))

  dsbuffer<-st_transform(st_as_sf(ds),crs=st_crs(g[[1]]))
  # the fun is to ensure that a tile is returned in case values are missing in the raster
  # assumes that the most common value is the tile number
  # it is possible the in some case the value returned will be in a neigbouring tile
  e<-extract(g,dsbuffer,buffer=2000,fun=function(i){names(rev(sort(table(i)))[1])})
  dsbuffer$tile<-e
  dsbuffer$year_tile<-paste(dsbuffer$year,dsbuffer$tile,sep="_")

  lf<-list.files(file.path(path,"daymet"),full=TRUE,pattern=weathervars[v])
  lf<-lf[order(lf)]
  ids<-do.call("rbind",strsplit(lf,"_"))
  o<-order(ids[,3],ids[,2])
  ids<-ids[o,]
  lf<-lf[o]
  lf<-split(lf,paste(ids[,2],gsub(".nc","",ids[,3]),sep="_"))
  lf<-lapply(lf,function(i){
    stack(i)
  })
  nameslf<-names(lf)
  
  if(weathervars[v]%in%c("anom","tmean")){ # gives correct names to these generated vars
    lf<-lapply(seq_along(lf),function(i){
      #r<-aggregate(i,2)
      r<-lf[[i]]
      names(r)<-paste0("X",substr(names(lf)[i],1,4),gsub("-",".",substr(seq.Date(as.Date("2010-01-01"),as.Date("2010-12-31"),1),5,10)))
      r
    })
    names(lf)<-nameslf
  }
  
  
  
  #lf2<-lapply(lf[1],function(i){
  #  # this fills missing values with the mean of neighbouring cells
  #  lras<-lapply(1:nlayers(i),function(j){
  #    focal(i[[j]],w=matrix(rep(1,3^2),ncol=3),fun=mean,na.rm=TRUE,NAonly=TRUE)
  #  })
  #  stack(lras)
  #})

  # not working cause needs too much ram
  registerDoParallel(detectCores()-4) 
  getDoParWorkers()

  ptemps<-c(2,7,30,90) # number of days over which to average weather values

  #l<-lapply(seq_along(lf)[1:5],function(i){
  l<-foreach(i=seq_along(lf),.packages=c("raster","sf")) %do% { # parallel version needs too much ram 
    w<-which(dsbuffer$year_tile==names(lf)[i])
    if(any(w)){
      e<-extract(lf[[i]],dsbuffer[w,],method="bilinear") # some cells with NA so this takes neighbouring cells  
      ll<-lapply(seq_along(w),function(j){
        m<-match(dsbuffer$date[w[j]],gsub("\\.","-",gsub("X","",dimnames(e)[[2]])))
        a<-sapply(ptemps,function(k){
          mean(e[j,(m-k):m])
        })
        a
      })
      ll<-do.call("rbind",ll)  
      #plot(lf[[i]][[1]])
      #plot(st_geometry(dsbuffer[w,]),add=TRUE)
      dimnames(ll)[[1]]<-w
      cat("\r",paste(weathervars[v],i,length(lf),sep=" / "))
      ll
    }
  }

  vals<-do.call("rbind",l)
  vals<-vals[order(as.integer(dimnames(vals)[[1]])),]
  vals<-as.data.frame(vals)
  names(vals)<-paste0(weathervars[v],ptemps)
  ds<-cbind(ds,vals)
  
}

#plot(st_geometry(st_transform(st_as_sf(ds),crs=crs(lf[[1]]))))
#plot(st_geometry(st_transform(st_as_sf(Q),crs=crs(lf[[1]]))),axes=TRUE,add=TRUE)
#lapply(1:length(lf),function(i){
#  plot(lf[[i]][[1]],legend=FALSE,add=F,zlim=c(-15,5))
#  text(colMeans(coordinates(lf[[i]]))[1],colMeans(coordinates(lf[[i]]))[2],i,cex=5)
#  plot(st_geometry(st_transform(st_as_sf(ds[!1:nrow(ds)%in%as.integer(dimnames(vals)[[1]]),]),crs=crs(lf[[1]]))),add=TRUE)
#  Sys.sleep(1)
#})
#plot(lf[[5]][[6]],add=TRUE,legend=FALSE)
#plot(st_geometry(st_transform(st_as_sf(ds),crs=crs(lf[[1]]))),add=TRUE)
#par(mfrow=c(1,2))
#test<-lf[[8]]
#plot(test,zlim=c(-10,0))
#test<-disaggregate(test,fact=2,method="bilinear")
#plot(test,zlim=c(-10,0))
#par(mfrow=c(1,1))

#plot(ds)
#plot(ds[!1:nrow(ds)%in%as.integer(dimnames(vals)[[1]]),],add=TRUE,col="red",pch=1,lwd=3)
#plot(st_geometry(st_transform(st_as_sf(ds[!1:nrow(ds)%in%as.integer(dimnames(vals)[[1]]),]),crs=crs(lf[[1]]))),add=TRUE)

#rr<-lf[[8]]#[[1]]
#par(mfrow=c(1,2))
#plot(rr)
#rr<-focal(rr,w=matrix(rep(1,3^2),ncol=3),fun=mean,na.rm=TRUE,NAonly=TRUE)
#plot(rr)
#par(mfrow=c(1,1))

#ll<-lapply(1:nlayers(rr),function(i){
#  focal(rr[[i]],w=matrix(rep(1,3^2),ncol=3),fun=mean,na.rm=TRUE,NAonly=TRUE)
#})


#e<-exact_extract(weather,buff)
#temp<-sapply(e,function(i){
#  w<-sample(1:(ncol(i)-1),1)
#  sum(i[,w]*i[,"coverage_fraction"])/sum(i[,"coverage_fraction"])
#})

#x<-1:5
#p<-c(0.99,0.99,0.005,0.005,0.005)
#sum(x*p)/sum(p)


plot(ds@data[ds$db!="map",unique(c("jul",names(ds)[grep("tmean|prcp",names(ds))]))])


### LULC data #########################################

l<-list.files(file.path(path,"LULC/LULC/"),pattern=".tif",full.names=TRUE)
l<-l[substr(l,nchar(l)-3,nchar(l))==".tif"]
lulc<-stack(l)
water<-lulc[["LULC2011"]]
water<-rast(water)
water<-crop(water,st_transform(st_as_sf(mappingzone),crs=crs(water)))
water[water!=20]<-NA
water<-ms_simplify(st_union(st_transform(st_as_sf(as.polygons(water)),crs=crs(ds))),0.05)
#NAvalue(lulc)<-0

# build ids to get unique location/year
ds$idlulc<-paste(ds$longitude,ds$latitude)

# make buffers of different widths for unique ids
bwidth<-c(50,1000) # list any number of widths
lbuffer<-lapply(bwidth,function(i){
  st_buffer(st_as_sf(spTransform(ds[!duplicated(ds$idlulc),],CRS(proj4string(lulc)))),i)  
})
dsbuffer<-do.call("rbind",lbuffer)

### visualize buffers on a satellite map
buffers<-st_transform(dsbuffer,4326)
#mapview(list(st_centroid(buffers),buffers),alpha=0.25)


# crop raster for faster computations
lulc<-crop(lulc,extent(st_bbox(dsbuffer)))

e<-exact_extract(lulc[["LULC2011"]],dsbuffer)

l<-lapply(seq_along(e),function(i){
  cov<-data.table(classn=e[[i]][,"value"],area=e[[i]][,"coverage_fraction"]*res(lulc)[1]*res(lulc)[2])
  a<-cov[,.(area=sum(area)),by=.(classn)]
  a[,area:=area/sum(area)]
  a[,idlulc:=dsbuffer$idlulc[i]]
  a[,bwidth:=rep(bwidth,each=length(e)/length(bwidth))[i]]
  cat("\r",paste(i,length(e),sep=" / "))
  a
})
l<-do.call("rbind",l)
ld<-setDT(l)
pcov<-dcast(ld,idlulc+bwidth~classn,value.var="area",fill=0)
table(rowSums(pcov[,-(1:2)])) # should expect only 1's

lccnames<-as.data.frame(read_excel(file.path(path,"0_ListeReclassification.xlsx"),sheet="Reclassification",range="C52:D64"))
names(lccnames)<-c("classn","orig")
lccnames<-cbind(lccnames,class=c("clouds","water","barren","urban","wet","pond","marsh","swamp","crop","pasture","shrub","forest"))
names(pcov)[-1]<-lccnames$class[match(names(pcov)[-1],lccnames$classn)]

# combine some lcc
pcov[,wetland:=wet+swamp+pond+marsh]
pcov[,agriculture:=crop+pasture]
pcov[,natural:=shrub+forest]

# set names for each buffer size
pcov<-lapply(split(pcov,pcov$bwidth),function(i){
  setorder(i,idlulc)
  names(i)[3:ncol(i)]<-paste0(names(i)[3:ncol(i)],i$bwidth[1])
  i[,bwidth:=NULL]
  i
})
# merge them to get a column by lcc/width
pcov<-Reduce(function(...) merge(..., all=T),pcov)

### merge with obs
ds<-merge(ds,pcov,by="idlulc",all.x=TRUE)

#ds$wetland<-ds$wet+ds$swamp+ds$marsh+ds$pond
#ds$agriculture<-ds$crop+ds$pasture
#ds$natural<-ds$shrub+ds$forest


#### Map LULC data #########################################

xlim<-bbox(mappingzone)[1,]*1000
ylim<-bbox(mappingzone)[2,]*1000
#arg <- list(at=rat$ID, labels=rat$class)
cols<-c("gray","skyblue","grey20","grey50","brown","skyblue","brown","brown","lightgoldenrod","lightgoldenrod","green","forestgreen")
par(mar=c(0,0,0,8))
plot(lulc[[1]],col=cols,breaks=c(0,lccnames$classn),axis.args=arg,xlim=xlim,ylim=ylim,zlim=c(0,220),legend=FALSE,legend.width=1,legend.shrink=1.25,axes=TRUE,bty="n")
legend(x=par("usr")[2],y=par("usr")[4],legend=paste(lccnames$class,lccnames$classn),fill=cols,bty="n",border=NA,cex=1.2,xpd=TRUE)
plot(st_geometry(st_buffer(st_transform(st_as_sf(ds[ds$id!="map",]),proj4string(lulc)),1000)),add=TRUE)
#loc<-locator()
#plot(lulc[[1]],col=cols,breaks=c(0,lccnames$classn),axis.args=arg,xlim=range(loc$x),ylim=range(loc$y),zlim=c(0,220),legend=FALSE,legend.width=1,legend.shrink=1.25,axes=TRUE)
#plot(st_geometry(st_buffer(st_transform(st_as_sf(ds),proj4string(lulc)),1000)),add=TRUE)
#legend(x=par("usr")[2],y=par("usr")[4],legend=paste(lccnames$class,lccnames$classn),fill=cols,bty="n",border=NA,cex=1.2,xpd=TRUE)

ds$jul<-ds$jul/100
ds$julsquare<-ds$jul^2
ds$db<-gsub("pred","map",ds$db)

options(scipen=20)
rev(sort(sapply(ls(),function(i){object.size(get(i))})))/1024^2

rm(e,lf,dsbuffer,can,map,info,coords,lbuffer,buffers,que,inspq,gdg);gc();gc()

##### Check if lcc OK ############################################

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(ds[ds$db=="map",],pch=21,cex=0.5,col=NA,bg=gray(scales::rescale(ds$urban1000[ds$db=="map"],c(0.05,0.99))))

##### Show locations and % #######################################
#x<-ds[!duplicated(ds$idlulc),]
#x<-x[x$db!="map",]
#x<-st_as_sf(x)
#x50<-st_transform(st_buffer(x,0.050),4326)
#x1000<-st_transform(st_buffer(x,1),4326)
#x<-st_transform(x,4326)
#vars<-names(ds)[grep("1000|50",names(ds))]
#mapviewOptions(basemaps = c("Esri.WorldImagery"))
#mapview(list(x50[,vars],x[,vars]),layer.name=c("barren50 buffer50m","trap"),zcol="barren50",alpha.regions=0.5)
#mapview(list(x1000[,vars],x[,vars]),layer.name=c("barren1000 buffer1000m","trap"),zcol="barren1000",alpha.regions=0.5)


### save the loaded data in a session
#save.image("mosquitos2.RData")

## Correlations ###############################################

#### Pairwise #################################################
vars<-names(ds)[grep("jul|tmean|tmin|tmax|prcp|anom|forest|agriculture|water|urban|pond|swamp|pasture|crop|wet|barren",names(ds))]
vars<-vars[-grep("CQ",vars)]
corrplot(cor(ds@data[ds@data$db!="map",vars]),method="number",number.cex=0.5)

par(mar=c(0,0,0,0),oma=c(4,4,1,1))
plot(ds@data[ds@data$db!="map",vars[grep("jul|tmean|anom|prcp",vars)]],pch=16,cex=0.5,col=gray(0,0.1))

par(mar=c(0,0,0,0),oma=c(4,4,1,1))
plot(ds@data[ds@data$db!="map",vars[grep("tmean|tmin|tmax",vars)]],pch=16,cex=0.5,col=gray(0,0.1))

ws<-vars[grep("jul|tmean|prcp|anom",vars)]
#ws<-vars[grep("jul|anom",vars)]
ws<-ws[-grep("square",ws)]
par(mar=c(0,0,0,0),oma=c(4,4,1,1))
plot(ds@data[ds@data$db!="map",ws],pch=16,cex=0.5,col=gray(0,0.075))

ws<-vars[grep("agriculture|forest|shrub|urban|wetland",vars)]
par(mar=c(0,0,0,0),oma=c(4,4,1,1))
plot(ds@data[ds@data$db!="map",ws],pch=16,cex=0.5,col=gray(0,0.075))
cor(ds@data[ds@data$db!="map",ws])

#### VIF ######################################################

### check vif for specific model
mo<-lm(VEX_Aedes_vexans~jul+urban1000+forest1000+wetland1000+tmax7+prcp7,data=ds@data[ds@data$db!="map",])
vif(mo)
hist(unique(apply(ds@data[ds@data$db!="map",names(model.frame(mo)[,-1])],1,sum)),xlim=0:1)

### check vif for all combinations of lcc and climate
ws1<-vars[grep("agriculture|forest|shrub|urban|wetland",vars)]
mo1<-combn(ws1,3,simplify=FALSE)
ws2<-vars[grep("anom|prcp",vars)]
mo2<-combn(ws2,3,simplify=FALSE)
ee<-expand.grid(1:length(mo1),1:length(mo2))
mo<-lapply(1:nrow(ee),function(i){c(mo1[[ee[i,1]]],mo2[[ee[i,2]]])})
#ws<-vars[grep("agriculture|forest|shrub|urban|wetland|anom|prcp",vars)]
#mo<-combn(ws,4,simplify=FALSE)
registerDoParallel(detectCores()-1) 
getDoParWorkers()
v<-foreach(i=seq_along(mo),.packages=c("car","sp")) %dopar% {
  f<-formula(paste("VEX_Aedes_vexans~jul+",paste(mo[[i]],collapse="+")))
  mod<-lm(f,data=ds@data[ds@data$db!="map",])
  #print(i)
  vif(mod)  
}
v<-as.data.frame(do.call("rbind",v))
names(v)<-paste0("vif",1:length(mo[[1]]))
v<-cbind(model=sapply(mo,paste,collapse=" + "),v)
v[apply(v,1,function(i){any(i[-1]>3)}),]


vif(lm(VEX_Aedes_vexans~urban1000+forest1000+agriculture1000,data=ds@data[ds@data$db!="map",]))


library(MASS)
library(glmnet)
d<-ds@data[ds@data$db!="map",]
fit <- cv.glmnet(scale(as.matrix(d[,vars])), d$VEX_Aedes_vexans, family = negative.binomial(theta = 3))


## MODELS #####################################################

#load("mosquitos.RData")

# summary of most abundant species per year
#d[,lapply(.SD,sum,na.rm=TRUE),by=year,.SD=species][order(year),][,1:6]

#### Model list ##########################################

climate<-list(
  ~ anom2 + prcp2,
  ~ anom7 + prcp7,
  ~ anom30 + prcp30,
  ~ anom90 + prcp90,
  ~ anom2 + prcp2 + anom7 + prcp7,
  ~ anom2 + prcp2 + anom30 + prcp30,
  ~ anom2 + prcp2 + anom90 + prcp90
)

lcc<-list(
  VEX=list(
    ~ agriculture50 + forest50,
    ~ agriculture1000+ forest1000,
    ~ agriculture50 + forest50 + agriculture1000+ forest1000
  ),
  CPR=list(
    ~ urban50 + forest50,
    ~ urban1000+ forest1000,
    ~ urban50 + forest50 + urban1000+ forest1000
  ),
  CQP=list(
    ~ agriculture50 + forest50,
    ~ agriculture1000+ forest1000,
    ~ agriculture50 + forest50 + agriculture1000+ forest1000
  ),
  STM=list(
    ~ wetland50 + forest50,
    ~ wetland1000+ forest1000,
    ~ wetland50 + forest50 + wetland1000+ forest1000
  )
)

models<-lapply(names(lcc),function(n){
  i<-lcc[[n]]
  ex<-expand.grid(seq_along(i),seq_along(climate))
  res<-lapply(1:nrow(ex),function(j){
    formula(paste(
      "y ~ -1 + ns(jul,knots=knots)",
      paste(paste(all.vars(i[[ex[j,1]]]),collapse=" + "),paste(all.vars(climate[[ex[j,2]]]),collapse= " + "),sep=" + "),
      paste("offset(lognights) + f(spatial, model=spde)"),
      sep=" + "
    ))  
  })
  names(res)<-paste(n,1:nrow(ex),sep="_")
  res
})
names(models)<-names(lcc)


#### VIFS #################################################

vifs<-lapply(unlist(models),function(i){
  v<-all.vars(i)
  v<-v[!v%in%c("y","knots","intercept","lognights","spatial","spde")]
  f<-formula(paste("VEX_Aedes_vexans~jul+",paste(v,collapse="+")))
  mod<-lm(f,data=ds@data[ds@data$db!="map",])
  #print(i)
  vif(mod)  
})
vifs[sapply(vifs,function(i){any(i>3)})]


#### Subset data #############################################
cat("\014")
inla.setOption(inla.mode="experimental")
year<-c(2003:2016);
weeks<-10:50
spcode<-"VEX_"
lweeks<-lapply(year,function(i){list(i,weeks)})
#lweeks<-list(list(2014,weeks),list(2015,29:32))
weeks<-apply(do.call("rbind",lapply(lweeks,function(i){expand.grid(year=i[[1]],week=i[[2]])})),1,function(i){paste(i[1],i[2],sep="_W")})
xs<-ds[ds$year%in%year,]
sp<-names(xs)[grep(spcode,names(xs))]
xs$sp<-xs@data[,sp]
xs<-xs[xs$week%in%weeks,]
xs<-xs[order(xs$week),]
a<-aggregate(sp~week,data=xs@data,FUN=function(i){length(i)})
names(a)[2]<-"nbtraps"
xs$nbtraps<-a$nbtraps[match(xs$week,a$week)]
xs<-xs[xs$nbtraps>=0,]

xs$jul<-as.integer(format(as.Date(xs$date),"%j"))

#meanjul<-mean(xs$jul[xs$db!="map"])
#sdjul<-sd(xs$jul[xs$db!="map"])

#### Scale variables ########################

vs<-unique(unlist(lapply(unlist(models),all.vars)))
vs<-vs[!vs%in%c("lognights","spde","knots","y","spatial")]
vscale<-lapply(xs@data[xs$db!="map",vs],function(i){c(mean=mean(i),sd=sd(i))})
xs@data[vs]<-lapply(vs,function(i){(xs@data[[i]]-vscale[[i]][1])/vscale[[i]][2]})
bscale<-function(x,v=NULL){
  (x*vscale[[v]]["sd"])+vscale[[v]]["mean"]
}

#xs$jul<-(xs$jul-meanjul)/sdjul
xs$julsquare<-xs$jul^2
xs$lognights<-log(xs$nights) # for offset
#xs$wetland50<-log(xs$wetland50+0.5)
#xs$wetland1000<-log(xs$wetland1000+0.5)

l<-split(xs,xs$week)
par(mfrow=n2mfrow(length(l)),mar=c(0,0,0,0))
lapply(l,function(i){
  plot(mappingzone)
  plot(i,cex=i$sp/200,pch=1,add=TRUE)
  mtext(side=3,line=-1,text=i$week[1])
})
par(mfrow=c(1,1),mar=c(0,0,0,0))


#w<-which(xs$sp>5000000)
#if(any(w)){
#  xs<-xs[-w,]
#}

#xs<-xs[substr(xs$week,7,8)%in%c(10:23,25:38),]

#xs$week[xs$db!="map"]<-sample(xs$week[xs$db!="map"])


#### Mesh #####################################################

if(TRUE){
  xs2<-st_as_sf(xs)
  xs2map<-xs2[xs2$db=="map",]
  xs2<-xs2[xs2$db!="map",]
  
  predmap<-st_buffer(concaveman(xs2,2.5),0)
  plot(st_geometry(predmap))
  plot(st_geometry(xs2),add=TRUE)
  xs2map<-xs2map[predmap,]
  plot(st_geometry(xs2map),add=TRUE)
  
  #set.seed(1234)
  xs3pts<-as(st_sample(predmap,3000),"Spatial")
  xs2pts<-as(st_cast(predmap,"MULTIPOINT"),"Spatial")
  xs2<-as(xs2,"Spatial")
  xs2map<-as(xs2map,"Spatial")
  
  plot(st_geometry(predmap),add=TRUE,border="red")
  
  plot(xs2pts,add=TRUE)
  plot(xs3pts,add=TRUE)
  
}


edge<-2
#domain <- inla.nonconvex.hull(coordinates(ds),convex=-0.015, resolution = c(100, 100))
#mesh<-inla.mesh.2d(loc.domain=coordinates(ds),max.edge=c(edge,3*edge),offset=c(edge,1*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
#domain <- inla.nonconvex.hull(coordinates(xs2pts),convex = -0.15, concave = 0.5, resolution = c(340,340))
domain <- inla.nonconvex.hull(rbind(coordinates(xs2pts),coordinates(xs3pts)),convex = -0.05)
#domain <- inla.nonconvex.hull(rbind(coordinates(xs2pts)),convex = -0.12)
#ims<-inla.mesh.segment(loc=coordinates(xs2pts))
mesh<-inla.mesh.2d(loc.domain=NULL,max.edge=c(edge,3*edge),offset=c(edge,3*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
#mesh<-inla.mesh.2d(loc=coordinates(xs),max.edge=c(2,8),cutoff=2)
#mesh<-inla.mesh.2d(boundary=domain,max.edge=c(edge,2*edge),offset=NULL,cutoff=0.5*edge,crs=CRS(proj4string(xs)))
plot(mesh,asp=1)
plot(Q,col="grey90",border="white",add=TRUE,lwd=2)
plot(mesh,asp=1,add=TRUE)
plot(xs[!xs$db%in%"map",],add=TRUE,pch=16,col="forestgreen")
#plot(mappingzone,add=TRUE)
mtext(side=3,line=-2,text=paste("edge =",edge,"km"),font=2,cex=1.5)


#### Restrict predictions ######################################
xsmap<-xs[xs$db=="map",]
xs<-xs[xs$db!="map",]

#### SPDE #################################################
spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  constr = FALSE, # not exactly sure what this does
  prior.range=c(10, 0.01), ### P(practic.range<0.05)=0.01
  prior.sigma=c(0.5,0.01)) ### P(sigma>1)=0.01

#### Priors on hyperpar ##################################
#h.spec <- list(theta=list(prior="pc.prec", param=c(0.5,0.5)), rho=list(prior="pc.cor1", param=c(0.9,0.9)))
#h.spec <- list(theta = list(prior="pc.prec", param=c(1, NA)),
#               rho = list(prior="pc.cor0", param=c(0.1, NA)))

# priors on zero inflation probability
vals<-rnorm(10000,-2,1.5)
hist(inla.link.invlogit(vals))
hist(exp(vals)/(1+exp(vals)))


# priors on nb size parameter (not working check inla docs...)
theta<-7
vals<-rgamma(10000,exp(-theta),exp(-theta))
par(mfrow=c(1,2))
hist(vals,main="Size parameter for Negative Binomial")
mu<-5
hist(rnbinom(10000,mu=mu,size=exp(vals)),main=paste("Simulated counts with mu =",mu))
par(mfrow=c(1,1))
inla.pc.dgamma(x, lambda = 1, log = FALSE)

#co<-seq(-0.999,0.999,by=0.001)
#u<-seq(0.00001,0.2,length.out=10)
#alpha<-sqrt((1-u)/2)+0.1
#alpha<-rep(0.75,length(u))

#u<-c(0.0000001,0.01,0.1,0.5,0.9,0.95,0.99,0.99)
#sqrt((1-u)/2)
#alpha<-c(0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.999)
#l<-lapply(seq_along(u),function(i){
#  print(i)
#  scales::rescale(inla.pc.dcor1(co,u=u[i],alpha=alpha[i]),to=c(0,1))
#})
#plot(0:1,0:1,type="n",xlim=c(-1,1),ylim=range(unlist(l)))
#lapply(l,function(i){
#  lines(co,i)  
#})
#abline(h=0,lty=3)


#u<-c(0.001,0.01,0.1,0.5,0.7,0.9,0.95,0.99)
#alpha<-c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
#l<-lapply(seq_along(u),function(i){
#  print(i)
#  scales::rescale(inla.pc.dcor0(co,u=u[i],alpha=alpha[i]),to=c(0,1))
#})
#plot(0:1,0:1,type="n",xlim=c(-1,1),ylim=range(unlist(l)))
#lapply(l,function(i){
#  lines(co,i)  
#})
#abline(h=0,lty=3)

#jul<-c(1,25,75,runif(100,1,100),100)
#sp<-rnorm(length(jul))
#x<-data.frame(sp=sp,jul=jul)
#knots<-seq(min(xs$jul),max(xs$jul),length.out=4)
#unname(head(model.matrix(sp~ns(jul,df=NULL,knots=knots,Boundary.knots=range(x$jul)),data=x)))


#### Model formula ########################################
#model <- y ~ -1 + intercept + jul + julsquare + forest50 + urban50 + urban1000 + agriculture1000  + tmax7 + tmax2 + prcp30 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))

# using knots allow to fix the values for each jul across newdata
knots<-seq(min(xs$jul)+0.5,max(xs$jul)-0.5,length.out=7)

model<-models[["VEX"]][[21]]
#model<- y ~ -1 + ns(jul, knots = knots) + offset(lognights) + f(spatial, model = spde)
#model<- y ~ -1 + ns(jul,knots=knots) + urban50 + forest50 + urban1000 + forest1000 + anom2 + prcp2 + anom90 + prcp90 + offset(lognights) + f(spatial, model = spde, group = spatial.group, control.group = list(model = "ar1",hyper = h.spec))

#model <- y ~ -1 + ns(jul,knots=knots) + urban50 + urban1000+ forest50 + forest1000 + anom2 + prcp2 + anom30 + prcp30 + offset(lognights) + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + bs(jul) + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + jul + julsquare + forest1000 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + jul + julsquare + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#model <- y ~ -1 + intercept + jul + julsquare + forest100 + urban100 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#model <- y ~ -1 + intercept + jul + julsquare + forest + urban + tmax15 + tmax1
#model <- y ~ -1 + intercept + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#formulae <- y ~ 0 + w + f(spatial, model=spde) + f(week,model="rw1")
#formulae <- y ~ 0 + w + f(spatial, model=spde, group=spatial.group,control.group=list(model='exchangeable')) 





v<-setdiff(all.vars(model),c("y","spatial","intercept","spde","year","knots","lognights")) 

#### Priors on fixed effects ##############################
vals<-list(intercept=1/35^2,default=1/35^2) #5-30
control.fixed<-list(prec=vals,mean=list(intercept=0,default=0),expand.factor.strategy = "inla")

#### Newdata ########################################
n<-100
v2<-v[grep("square",v)]
v1<-setdiff(v,v2)
lp<-newdata(x=xs@data[,v1,drop=FALSE],v=v1,n=n,fun=mean,list=FALSE)
#lp<-newdata(x=rbind(xs@data[,v1,drop=FALSE],xsmap@data[,v1,drop=FALSE]),v=v1,n=n,fun=mean,list=FALSE) # this "partially?) takes the range of observed and predicted, need to better this and make it general to all variables)
if(length(v2)){
  lp<-lapply(lp,function(i){
    a<-as.data.frame(lapply(i[,gsub("square","",v2),drop=FALSE],"^",2))
    names(a)<-v2
    res<-cbind(i,a)
    res[,order(names(res))]
  })
}
if(TRUE){ # adds basis columns to new data
  lp<-lapply(names(lp),function(j){
    i<-lp[[j]]
    if(j=="jul"){
      mm<-model.matrix(jul~0+ns(jul,knots=knots),data=i)
    }else{
      mm<-model.matrix(jul~0+ns(jul,knots=knots),data=lp[["jul"]])
      mm<-mm[which.min(abs(lp[["jul"]]$jul-mean(xs$jul))),,drop=FALSE][rep(1,nrow(i)),] # finds the smallest difference with the jul date sequence to get basis jul values
    }
    res<-cbind(i,mm)
    names(res)<-gsub("ns\\(jul\\, knots = knots\\)","X",names(res))
    res
  })
}
lp<-lapply(lp,function(i){cbind(intercept=1,i)}) # add intercept to lp
names(lp)<-v1
#head(lp[[1]])
#model<-formula(paste("y~-1+",paste(v,collapse="+")))
#lp<-lapply(lp,function(i){mmatrix(model,i)})
#lpmed<-mmatrix(model,newdata(x=xs[,v,drop=FALSE],v=v,n=1,fun=median,list=FALSE)[[1]][rep(1,length(g)),])[[1]]


#### Make index ##########################################
iset<-inla.spde.make.index("spatial",n.spde=spde$n.spde)

#### A matrix ##############################################
Aest<-inla.spde.make.A(mesh=mesh,loc=coordinates(xs)) 
Amap<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap)) 

#### Stacks ################################################
stackest<-inla.stack(tag='est',data=list(y=xs$sp),A=list(Aest,1),effects=list(c(iset,list(intercept=1)),xs@data)) 
stackmap<-inla.stack(tag='map',data=list(y=xsmap$sp),A=list(Amap,1),effects=list(c(iset,list(intercept=1)),xsmap@data)) 
stackfull<-inla.stack(stackest,stackmap)


#### Stack vars ############################################
for(i in seq_along(v1)){
  le<-nrow(lp[[v1[i]]])
  if(le!=n){
    #AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,le),,drop=FALSE]) # for categorical variables
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(580,5045),ncol=2)[rep(1,n),,drop=FALSE])
  }else{
    #AA<-Apn # for numerical variables
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(580,5045),ncol=2)[rep(1,n),,drop=FALSE])
  }
  stack<-inla.stack(tag=v1[i],data=list(y=NA),A=list(AA,1),effects=list(c(iset,list(intercept=1)),lp[[v1[i]]][-1])) # -1 removes intercept in lp not ure if essential    
  stackfull<-inla.stack(stackfull,stack)
}

#### Full stack ############################################
#stackfull<-inla.stack(stackest)


#### Index #################################################
index.est<-inla.stack.index(stackfull,tag="est")$data
index.map<-inla.stack.index(stackfull,tag="map")$data
index<-list(est=index.est,map=index.map)
for(i in seq_along(v1)){
  index<-c(index,list(inla.stack.index(stackfull,tag=v1[i])$data))
}  
names(index)[3:length(index)]<-v1


#### Model ##################################################
m <- inla(model,data=inla.stack.data(stackfull), 
          control.predictor=list(compute=TRUE, A=inla.stack.A(stackfull),link=1), 
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=control.fixed,
          control.inla = list(strategy='gaussian',int.strategy = "eb"),
          num.threads="1:1",
          verbose=FALSE,
          control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
          #control.mode = list(result = m, restart = TRUE)), # to rerun the model with NA predictions according to https://06373067248184934733.googlegroups.com/attach/2662ebf61b581/sub.R?part=0.1&view=1&vt=ANaJVrHTFUnDqSbj6WTkDo-b_TftcP-dVVwK9SxPo9jmPvDiK58BmG7DpDdb0Ek6xypsqmCSTLDV1rczoY6Acg_Zb0VRPn1w2vRj3vzHYaHT8JMCEihVLbY
          family="nbinomial")#"zeroinflatednbinomial1"


#### Posterior samples ####################################

# from haakon bakk a, BTopic112
nsims<-500
samples<-inla.posterior.sample(nsims,m,num.threads="1:1")
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor" # not sure if should use APredictor or Predictor
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-do.call("cbind",samples.effect)
nbsize<-sapply(samples, function(x) x$hyperpar[grep("size",names(x$hyperpar))])
zeroprob<-sapply(samples, function(x) x$hyperpar[grep("probability",names(x$hyperpar))])

#### Hyperpar samples #####################################
# I don't think this works, hyperpars cannot be sampled independently of fixed effects and everything should be sampled jointly, this will overestimate variability
sampleshyper<-inla.hyperpar.sample(nsims,m)
nbsize<-sampleshyper[,grep("size",dimnames(sampleshyper)[[2]])]
zeroprob<-sampleshyper[,grep("probability",dimnames(sampleshyper)[[2]])]


#save.image("mosquitos_model.RData")
#load("mosquitos_model.RData")

#load("mosquito_models_do.RData")

## RESULTS ################################################

#### Show posteriors ########################################

posteriors<-c(m$marginals.fixed,m$marginals.hyper)
par(mfrow=n2mfrow(length(posteriors),asp=1.49),mar=c(3,2.5,1,1))
for (j in 1:length(posteriors)) {
  k<-posteriors[[j]][,2]>=1e-5*max(posteriors[[j]][,2])
  posteriors[[j]]<-posteriors[[j]][k,]
}
for (j in 1:length(posteriors)) {
  plot(posteriors[[j]][,1],posteriors[[j]][,2],type='l',xlab=names(posteriors)[j],ylab='Density',tcl=-0.2,mgp=c(1.5,0.45,0),xlim=NULL)#range(do.call("rbind",posteriors)[,1]))
}
par(mfrow=c(1,1))



#### Spatial fields #########################################

#### Projection grid 

stepsize <- 0.5*1/1
coords<-st_coordinates(st_cast(st_buffer(st_as_sf(xs),10),"MULTIPOINT"))
nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
projgrid <- inla.mesh.projector(mesh, xlim=range(coords[,1]),ylim=range(coords[,2]), dims=nxy,crs=CRS(proj4string(xs)))

#### Extract mean and sd 
vfield<-c("mean","sd")
field<-list()
for(i in seq_along(vfield)){
  xmean <- inla.mesh.project(projgrid,m$summary.random$spatial[[vfield[i]]])
  #### Set NAs 
  b<-gBuffer(gConvexHull(SpatialPoints(domain$loc,p=CRS(proj4string(ds)))),width=0.1,byid=FALSE)
  o <- over(SpatialPoints(projgrid$lattice$loc,p=CRS(proj4string(ds))),b)
  xmean[is.na(o)] <- NA
  r<-raster(nrows=nxy[2], ncols=nxy[1], xmn=min(projgrid$x), xmx=max(projgrid$x), ymn=min(projgrid$y), ymx=max(projgrid$y),crs=CRS(proj4string(xs)),vals=as.vector(xmean[,ncol(xmean):1])) ## some crazy ordering in INLA output be careful
  names(r)<-vfield[i]
  field[[i]]<-r
}
names(field)<-vfield

#### Mask 
xsbuff<-st_coordinates(st_cast(st_buffer(st_as_sf(xs),7),"MULTIPOINT"))[,1:2]
buf<-concaveman(xsbuff,10)
buf<-spPolygons(buf,crs=CRS(proj4string(xs)))
buf<-gBuffer(buf,width=1)
#field<-stack(field)
r<-mask(field[["mean"]],buf)

#### Plot fields 
cols<-colo.scale(seq(range(values(r),na.rm=TRUE)[1],range(values(r),na.rm=TRUE)[2],length.out=200),c("darkblue","dodgerblue","ivory2","tomato2","firebrick4"),center=TRUE)#,"grey20"))
p.strip<-list(cex=0.65,lines=1,col="black")
levelplot(r,col.regions=cols,cuts=199,margin=FALSE,par.strip.text=p.strip,par.settings = list(axis.line = list(col = "grey90"),strip.background = list(col = 'transparent'),strip.border = list(col = 'grey90')),scales = list(col = "black")) +
  layer(sp.points(xs,col="black",pch=1,cex=scales:::rescale(c(max(xs$sp),identity(xs$sp)),to=c(0.5,15))[-1]))+
  layer(sp.points(xs,col="black",pch=3,cex=0.5))+
  layer(sp.polygons(as(water,"Spatial"),col=gray(0,0.2)))
par(mfrow=c(1,1))


#### Marginal effects ########################################

# page 263 in Zuur
table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  #grep(paste0(i,":"),row.names(samples[[1]]$latent))  # old version
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
#table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
nweights<-grep("spatial",row.names(samples[[1]]$latent))

par(mfrow=n2mfrow(length(v1),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
for(k in seq_along(v1)){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
    fixed<-as.matrix(lp[[v1[k]]][,names(betas)]) %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
       #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
      spatial<-as.matrix(AA) %*% wk # stack was Apn in fire
    #}
    #p<-fixed+spatial
    p<-fixed # ignores spatial part
    p
  })
  p<-do.call("cbind",p)
  p<-t(apply(p,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  p<-exp(p)
  if(nrow(lp[[v1[k]]])==n){
    vals<-bscale(lp[[v1[k]]][,v1[k]],v=v1[k])
    #if(v1[k]%in%c("wetland50","wetland1000")){vals<-exp(vals)-0.5}
    plot(vals,p[,2],type="l",ylim=c(0,min(c(max(p[,3]),max(xs$sp))))*1.3,xlab=v1[k],font=2,ylab="",lty=1,yaxt="n",mgp=c(2,0.45,0),tcl=-0.3)
    points(bscale(xs@data[,v1[k]],v=v1[k]),xs$sp,pch=1,col=gray(0,0.1))
    lines(vals,p[,2],lwd=3,col=gray(0,0.8))
    #lines(vals,p[,1],lty=3)
    #lines(vals,p[,3],lty=3)
    polygon(c(vals,rev(vals),vals[1]),c(p[,1],rev(p[,3]),p[,1][1]),col=alpha("black",0.1),border=NA)
  }else{
    #plot(unique(sort(size[,v[k]])),p[,2],type="l",ylim=c(0,100),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    #points(jitter(as.integer(size[,v[k]])),size$tTotal,pch=16,col=gray(0,0.1))
    #segments(x0=as.integer(unique(sort(size[,v[k]]))),x1=as.integer(unique(sort(size[,v[k]]))),y0=p[,1],y1=p[,3],lty=3)
  }
  axis(2,las=2)
}
mtext(paste("Mosquitos per trap"),outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)


#### Marginal effects with spatial uncertainty ##########################

# this section is not that useful because it is a prediction for a given location, hence it includes uncertainty in the spatial field

par(mfrow=n2mfrow(length(v1),asp=1.5),mar=c(3,3,1,1),oma=c(0,10,0,0))
for(i in seq_along(v1)){
  p<-m$summary.fitted.values[index[[v1[i]]],c("0.025quant","0.5quant","0.975quant")]
  #p[]<-lapply(p,transI)
  dat<-data.frame(lp[[v1[i]]])
  if(nrow(p)==n){
    vals<-bscale(dat[[v1[i]]],v=v1[i])
    plot(vals,p[,2],type="l",ylim=c(0,min(c(max(p[,3]),max(xs$sp)))),xlab=v1[i],font=2,ylab="",lty=1,yaxt="n",mgp=c(2,0.45,0),tcl=-0.3)
    #plot(dat[[v1[i]]],p[,2],type="l",ylim=c(0,300),xlab=v1[i],font=2,ylab="",lty=1,yaxt="n",mgp=c(2,0.45,0),tcl=-0.3)
    #lines(vals,p[,1],lty=3,lwd=1)
    #lines(vals,p[,3],lty=3,lwd=1)
    points(bscale(xs@data[,v1[i]],v=v1[i]),xs$sp,pch=16,col=gray(0,0.07))
    polygon(c(vals,rev(vals),vals[1]),c(p[,1],rev(p[,3]),p[,1][1]),col=alpha("black",0.1),border=NA)
  }else{
    #plot(unique(sort(size[,v1[i]])),p[,2],type="l",ylim=c(0,100),xlab=v1[i],font=2,ylab="",lty=1,yaxt="n")
    #segments(x0=as.integer(unique(sort(size[,v[i]]))),x1=as.integer(unique(sort(size[,v[i]]))),y0=p[,1],y1=p[,3],lty=3,lwd=2)
    #points(jitter(as.integer(size[,v[i]]),fac=2),transI(size$tTotal),pch=16,col=gray(0,0.07))
  }
  axis(2,las=2)
}
mtext("Weekly number of mosquitos",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)


#### Map predictions ####################################

# The code for the graph below really sucks... should make it better...

frange<-function(x,n=10){seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length.out=n)}

quantities<-c("mean","X0.025quant","X0.975quant","sd")
xsmappred<-cbind(xsmap[,"id"],data.frame(m$summary.linear.predictor[index.map,gsub("X","",quantities)]))
pred<-lapply(seq_along(quantities),function(i){
  rasterize(xsmappred,pgrid,field=quantities[i])
})
pred<-stack(pred)
names(pred)<-quantities
pred<-exp(pred) # transform back to counts after all calculations

#cw<-unique(xsmap$temporal) # watch out if there is more than one week predicted
#mcw<-match(paste0("X",cw),names(field[["mean"]]))
meansd<-stack(field)
meansd<-resample(meansd,pred[[1]])
names(meansd)<-c("mean.spatial.field","sd.spatial.field")
pred<-stack(pred,meansd)
#pred<-disaggregate(pred,fact=5,method="bilinear") # hack to make the map smoother

### use tighter mapping zone instead of mappingzone
xsbuff<-st_coordinates(st_cast(st_buffer(st_as_sf(xs),7),"MULTIPOINT"))[,1:2]
buf<-concaveman(xsbuff,10)
buf<-spPolygons(buf,crs=CRS(proj4string(xs)))
buf<-gBuffer(buf,width=1)
pred<-mask(pred,buf)

cols<-alpha(colo.scale(200,c("steelblue3","lightgoldenrod","orange","red3","darkred","grey10")),0.80)
colssd<-rev(cividis(200))
colsfield<-colo.scale(seq(range(values(pred[["mean.spatial.field"]]),na.rm=TRUE)[1],range(values(pred[["mean.spatial.field"]]),na.rm=TRUE)[2],length.out=200),c("navyblue","steelblue","ivory2","firebrick3","firebrick4"),center=TRUE)#,"grey20"))

par(mfrow=n2mfrow(nlayers(pred),asp=1.5),mar=c(1,0.5,1,5),oma=c(0,0,0,0),bty="n")
lapply(names(pred)[c(1,2,5,4,3,6)],function(i){
  print(i)
  if(i%in%names(meansd)){
    pred2<-pred
    f<-function(x){identity(x)}# exp vs identity
  }else{
    pred2<-log(pred)
    f<-function(x){exp(x)}# exp vs identity
  } # the .local warning NaNs produced comes from logging 0-negative values in mean/sd, not a problem
  col<-if(i%in%c("sd","sd.spatial.field")){
    colssd
  }else{
    if(i%in%"mean.spatial.field"){
      colsfield
    }else{
      cols
    }
  }
  #if(i%in%c("sd","mean")){
  if(i%in%quantities){
    zlim<-NULL
    axis.args=list(at=frange(values(pred2[[i]])),labels=round(f(frange(values(pred2[[i]]))),0),cex.axis=0.8,lwd=0,tck=-0.2,mgp=c(3,0.3,0),lwd.ticks=1)
    legend.args=list(text='Nb of mosquitos / trap', side=4, font=2, line=-2.5, cex=0.8)
  }else{
    if(i%in%names(meansd)){
      zlim<-NULL
      axis.args=list(at=frange(values(pred2[[i]])),labels=round(f(frange(values(pred2[[i]]))),2),cex.axis=0.8,lwd=0,tck=-0.2,mgp=c(3,0.3,0),lwd.ticks=1)
      legend.args=list(text='Nb of mosquitos / trap', side=4, font=2, line=-2.5, cex=0.8)
    }else{
      zlim<-range(values(pred2[[quantities[1:3]]]),na.rm=TRUE)
      axis.args=list(at=c(frange(zlim),range(values(pred2[[i]]),na.rm=TRUE)),labels=round(f(c(frange(zlim),range(values(pred2[[i]]),na.rm=TRUE))),0),cex.axis=0.8,lwd=0,tck=-0.2,mgp=c(3,0.3,0),lwd.ticks=1)
      legend.args=list(text='Nb of mosquitos / trap', side=4, font=2, line=-2.5, cex=0.8)
   }
  }
  plot(pred2[[i]],col=col,zlim=zlim,legend.width=2.5, legend.shrink=1,axis.args=axis.args,legend.args=legend.args,axes=FALSE,box=FALSE,tcl=0.2,mgp=c(1.5,0.0,0),cex.axis=0.7)
  plot(st_geometry(water),border=NA,col="white",add=TRUE)
  obs<-xs$sp[xs$week%in%unique(xsmap$week)]
  cexminmax<-c(1,10)  
  ocex<-scales::rescale(c(ifelse(obs==0,1,identity(obs))),to=cexminmax)
  opch<-ifelse(obs==0,4,1)
  plot(xs[xs$week%in%unique(xsmap$week),],add=TRUE,cex=ocex,pch=opch,lwd=1,col=gray(0.2,1))
  #plot(Q,add=TRUE,border=gray(0,0.25))
  #plot(mappingzone,add=TRUE)
  #plot(mesh,add=TRUE)
  mtext(side=3,line=-3,text=i,adj=0.05,font=2,cex=2)
  if(i %in% names(meansd)){
    plot(xs,add=TRUE,cex=0.5,pch=3,lwd=1,col=gray(0.2,1))
  }
  if(any(obs<1)){
    vobs<-frange(identity(obs[obs>0]),n=7)
    vcex<-c(1,scales::rescale(c(vobs),to=cexminmax))
    vleg<-round(c(0,identity(vobs)),0)
    vpch<-c(4,rep(1,length(vobs)))
  }else{
    vobs<-frange(identity(obs),n=8)
    vcex<-scales::rescale(c(vobs),to=cexminmax)
    vleg<-round(identity(vobs),0)
    vpch<-rep(1,length(vobs))
  }
  posx<-seq(par("usr")[1],par("usr")[2],length.out=length(vpch)+2)
  posx<-posx[-c(1,length(posx))]-diff(posx)[1]*0.25
  posy<-rep(diff(par("usr")[c(3,4)])*0.075+par("usr")[3],length(vpch))
  points(posx,posy,cex=vcex,pch=vpch,col=gray(0.2,1))
  text(posx,posy,label=vleg,adj=c(0.5,3))
})


#### Map predictions with posteriors samples #############################

# not done yet and not general enough
# this is mostly to make sure that the posterior sample approach gives the same results as the NA approach in the map stack

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  #grep(paste0(i,":"),row.names(samples[[1]]$latent))  
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
#table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
nweights<-grep("spatial",row.names(samples[[1]]$latent))
Amapmatrix<-as.matrix(Amap)

#par(mfrow=n2mfrow(length(v1),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
#for(k in seq_along(v1)){
p<-lapply(1:nsims,function(i){
  dat<-xsmap@data
  juls<-lp[["jul"]][which.min(abs(lp[["jul"]]$jul-dat$jul[1])),,drop=FALSE] # finds the closest jul value to get the corresponding sline basis
  dat<-cbind(dat,juls[,names(juls)%in%paste0("X",1:50)][rep(1,nrow(dat)),])
  betas<-samples[[i]]$latent[nparams]
  names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
  dat<-dat[,names(betas)]
  fixed<-as.matrix(dat) %*% betas # make sure betas and vars are in the same order
  # this if we want a spatial part
  wk<-samples[[i]]$latent[nweights]
  #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
  #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
  #}else{
  spatial<-Amapmatrix %*% wk
  #}
  p<-fixed+spatial
  #p<-fixed # ignores spatial part
  #p<-spatial
  print(i)
  p
})
p<-do.call("cbind",p)
p<-t(apply(p,1,function(i){c(quantile(i,0.0275,na.rm=TRUE),mean(i),quantile(i,0.975,na.rm=TRUE))}))
#p<-exp(p)

xsmap$preds<-p[,2]

pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
pr<-mask(pr,buf)
pr<-exp(pr)
#pr<-disaggregate(pr,fact=2,method="bilinear")

par(mfrow=c(1,2),oma=c(0,0,0,4))
f<-function(i){log(i)}
zlim<-range(f(c(values(pred[[1]]),values(pr))),na.rm=TRUE)
zlim<-range(f(c(values(pr))),na.rm=TRUE)
xxs<-xs[xs$week%in%xsmap$week[1],]
plot(f(pred[[1]]),zlim=zlim)
plot(st_geometry(water),border=NA,col="white",add=TRUE)
plot(xxs,add=TRUE,pch=1,cex=scales::rescale(xxs$sp,c(0.5,10)))
#plot(resample(pred[[1]],pr))
plot(f(pr),zlim=zlim)
plot(st_geometry(water),border=NA,col="white",add=TRUE)
plot(xxs,add=TRUE,pch=1,cex=scales::rescale(xxs$sp,c(0.5,10)))

#### Map predictions across season #############################

# not done yet and not general enough

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  #grep(paste0(i,":"),row.names(samples[[1]]$latent))  
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
#table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
yearpred<-"2003"
nweights<-grep("spatial",row.names(samples[[1]]$latent))
Amapp<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap)) 
Amapmatrix<-as.matrix(Amapp)
#Amapmatrix<-as.matrix(Amap) # previously
#Amapp1<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap),group=gss) 
#Amap1<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap),group=1)
#Amap3<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap),group=1)
#a1<-unique(as.vector(as.matrix(Amap3)[,1:212]))
#a2<-unique(as.vector(as.matrix(Amap3)[,213:424]))

#par(mfrow=n2mfrow(length(v1),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
#for(k in seq_along(v1)){

days<-seq(min(xs$jul[xs$year==yearpred]),max(xs$jul[xs$year==yearpred]),length.out=10)
lpr<-foreach(j=seq_along(days),.packages=c("raster")) %do% {
  juls<-lp[["jul"]][which.min(abs(lp[["jul"]]$jul-days[j])),,drop=FALSE]
  standardv<-names(nparams)[!names(nparams)%in%c("intercept","jul","julsquare",1:50)]
  dat<-as.matrix(xsmap@data[,standardv])
  dat<-cbind(dat,data.frame(jul=days[j],julsquare=days[j]^2)[rep(1,nrow(dat)),])
  dat<-cbind(dat,juls[,names(juls)%in%paste0("X",1:50)][rep(1,nrow(dat)),])
  dat<-cbind(intercept=1,dat)
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
    fixed<-as.matrix(dat[,names(betas)]) %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
    #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    spatial<-Amapmatrix %*% wk
    #}
    p<-fixed+spatial
    #p<-fixed # ignores spatial part
    #p<-spatial
    #print(i)
    p
  })
  p<-do.call("cbind",p)
  p<-t(apply(p,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  p<-exp(p)
  xsmap$preds<-p[,2]
  pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
  pr<-mask(pr,buf)
  #pr<-disaggregate(pr,fact=1,method="bilinear")
  print(j)
  pr
}

lpr<-lapply(lpr,rast)



#par(mfrow=c(1,2))
#plot(pred[[1]])

jul<-round(days*vscale[["jul"]]["sd"]+vscale[["jul"]]["mean"],0)
datelim<-range(as.Date(format(as.Date(jul,origin=paste0(yearpred,"-01-01")),"%Y-%m-%d")))+c(-5,5)

zlim1<-range(sapply(lpr,function(i){range(values(i),na.rm=TRUE)}))
zlim2<-range(c(sapply(lpr,function(i){range(values(i),na.rm=TRUE)}),xs$sp[xs$year==yearpred]))
zlim<-c(zlim1[1],zlim2[2])

img <- image_graph(1500, 1000, res = 150)
lapply(seq_along(lpr),function(i){
  j<-round(days[i]*vscale[["jul"]]["sd"]+vscale[["jul"]]["mean"],0)
  xdate<-as.Date(format(as.Date(j,origin=paste0(yearpred,"-01-01")),"%Y-%m-%d"))
  gw<-layout(matrix(c(1,rep(2,20)),ncol=1))
  par(mar=c(1,0,0,3),oma=c(0,0,0,5))
  plot(xdate,1,pch=25,xlim=datelim,yaxt="n",cex=2,bty="n",col=1,bg=1)
  par(mar=c(0,0,0,0))
  at<-seq(min(values(log(lpr[[i]])),na.rm=TRUE),max(values(log(lpr[[i]])),na.rm=TRUE),length.out=5)
  #lab<-ifelse(round(exp(at),0)==0,round(exp(at),2),round(exp(at),0))
  lab<-sapply(at,function(a){
    e<-exp(a)
    if(e<=0.001){return(round(e,4))}
    if(e<=0.01){return(round(e,3))}
    if(e<=0.1){return(round(e,2))}
    if(e<=2){return(round(e,1))}
    if(e>2){return(round(e,0))}
  })
  labels<-paste(lab,c("min pred. > 0",rep("",length(at)-2),"max pred."))
  plot(log(lpr[[i]]),range=log(zlim),col=cols,asp=1,axes=FALSE,bty="n",plg=list(at=at,labels=labels,cex=1.5))
  plot(st_geometry(water),border=NA,col="white",add=TRUE)
  rd<-as.character(seq.Date(xdate-3,xdate+3,by=1))
  xxs<-st_transform(st_as_sf(xs[xs$date%in%rd,]),crs=crs(lpr[[1]]))
  colobs<-c(log(zlim),log(xxs$sp))
  colobs<-ifelse(is.infinite(colobs),NA,colobs)
  colobs<-colo.scale(colobs,cols)[-(1:2)]
  colobs<-ifelse(is.na(colobs),"#FFFFFF",colobs)
  #plot(st_geometry(xxs),pch=1,cex=scales::rescale(identity(c(0,max(xs$sp[xs$year==yearpred],na.rm=TRUE),xxs$sp)+0.01),to=c(0.25,10))[-(1:2)],add=TRUE)
  plot(st_geometry(xxs),pch=21,cex=2,add=TRUE,bg=colobs,col="grey10",lwd=0.4)
  text(st_coordinates(st_geometry(xxs)),label=xxs$sp,cex=0.7,col="grey10",adj=c(0.5,-1))
  #plot(log(lpr[[i]]),zlim=log(zlim),col=cols,asp=1,legend.only=TRUE)
  mtext(side=3,line=-2,text=paste(gsub("_","",spcode),yearpred,"  observations:",paste(format(as.Date(range(rd)),"%b-%d"),collapse=" to "),sep="  "),adj=0.15)
  mtext(side=4,line=-1,text="Number of mosquitos per trap (observed and predicted)",adj=0.5)
  xp<-xmin(lpr[[1]])+((xmax(lpr[[1]])-xmin(lpr[[1]]))*c(0.58,0.65))
  yp<-rep(ymin(lpr[[1]])+((ymax(lpr[[1]])-ymin(lpr[[1]]))*0.97),2)  
  points(xp,yp,pch=21,cex=2,bg=colobs[c(which.min(xxs$sp),which.max(xxs$sp))],col="grey10",lwd=0.4)
  text(xp,yp,label=xxs$sp[c(which.min(xxs$sp),which.max(xxs$sp))],cex=0.7,col="grey10",adj=c(0.5,-1))
  text(xp,yp,label=c("min obs.","max obs."),cex=0.7,col="grey10",adj=c(1.2,0.5))
  ### hist of fit optional
  #print(i)
  #xxsb<-st_buffer(xxs,2)
  #e<-extract(lpr[[i]],vect(xxsb))
  #par(mar=c(2,0,0,0))
  #by<-50
  #maxn<-max(c(e[,2],xxs$sp))*1.05
  ##brks<-unique(c(seq(0,100,by=by),seq(100,200,by=by),seq(200,max(maxn,400)*1.05,by=by)))
  #brks<-seq(0,maxn+by,length.out=20)
  #h1<-hist(e[,2],breaks=brks,plot=FALSE)
  #h2<-hist(xxs$sp,breaks=brks,plot=FALSE)
  #ylim<-range(c(h1$density,h2$density))
  #h1<-hist(e[,2],breaks=brks,xlim=c(0,maxn),freq=FALSE,border=NA,ylim=ylim)
  #par(new=TRUE,lwd=3)
  #h2<-hist(xxs$sp,breaks=brks,xlim=c(0,maxn),freq=FALSE,col=NA,border=alpha("darkred",0.5),lwd=3,ylim=ylim)
  
})
dev.off()
animation <- image_animate(img, fps = 1, optimize = TRUE)
#print(animation)

image_write(animation,file.path("C:/Users/God/Downloads",paste0(paste0(spcode,yearpred),".gif")))



#### Model checks ##############################################

# make sure this is the right way to do it and check if paramete rs are ok. I'm sampling from the sampled hyperpar for each sims, but this is hacky and not correctly jointly sampled

# figure out if the predicted response used already includes the zero-inflation part. If so, the simulation of observations below is likely wrong. If not, meaning the predictor is the nbinom mean before zero-inflation, then the simulated value are probably ok.

#### Simulated
prob<-m$summary.fitted.values[index[["est"]],"mean"]
matprob<-apply(s.eff,2,function(i){
  #rbinom(length(i),size=1,prob=1-sample(zeroprob,1))*rnbinom(n=length(i),mu=exp(i),size=sample(nbsize,1))
  1*rnbinom(n=length(i),mu=exp(i),size=sample(nbsize,1)) # no zeroinflation
})

##### DHARMa plots
o<-createDHARMa(simulatedResponse=matprob,observedResponse=xs$sp,fittedPredictedResponse=prob,integerResponse=TRUE)
par(mfrow=c(2,2))
plotQQunif(o)
plotResiduals(o)
testZeroInflation(o)
testDispersion(o)
#hist(o$scaledResiduals)

##### Histograms 
par(mfrow=c(1,1))
brks<-unique(c(seq(0,100,by=5),seq(100,200,by=10),seq(200,max(matprob)*1.05,by=25)))
h1<-hist(matprob,breaks=brks,xlim=c(0,1000),freq=FALSE,border=NA)
ylim<-range(h1$density)
par(new=TRUE,lwd=3)
h2<-hist(xs$sp,breaks=brks,xlim=c(0,1000),freq=FALSE,col=NA,border=alpha("darkred",0.5),lwd=3,ylim=range(ylim))
par(lwd=1)

##### Plot simulated and observed 
plot(0,0,xlim=c(0,1000),ylim=c(0,max(h1$density)*1.07),type="n")
invisible(lapply(1:ncol(matprob[,1:100]),function(i){
  h<-hist(matprob[,i],breaks=brks,xlim=c(0,1000),freq=FALSE,plot=FALSE)
  points(h$mids,h$density,pch=16,cex=2,col=gray(0,0.02))
}))
points(h2$mids,h2$density,pch=16,cex=1.25,col=alpha("red",0.7))

#### Explanatory/Predictive power #############################

plot(m$summary.fitted.values[index[["est"]],1],xs$sp,asp=1)
cor(m$summary.fitted.values[index[["est"]],1],xs$sp)^2


#### SPDE posteriors ##########################################

res <- inla.spde.result(m, "spatial", spde)
par(mfrow=c(2,1))
plot(res$marginals.range.nominal[[1]],
     type="l", main="Posterior density for range")
plot(inla.tmarginal(sqrt, res$marginals.variance.nominal[[1]]),
     type="l", main="Posterior density for std.dev.")
par(mfrow=c(1,1))


### Trap variograms ###########################################

# Check structure (~variogram) in traps in relation to distance 
# to determine if there is a microhabitat/trap placement effect

l<-split(dsbuffer,dsbuffer$date)
l<-lapply(seq_along(l),function(i){
  print(i)
  x<-l[[i]]
  s<-st_distance(st_centroid(x))
  e<-expand.grid(trap1=1:nrow(s),trap2=1:nrow(s))
  sp<-"Aedes_vexans"#Culex_pipiens_restuans_gr"#"Coquillettidia_perturbans"#"Aedes_vexans"
  e$diff<-abs(x[[sp]][e[,"trap1"]]-x[[sp]][e[,"trap2"]])
  e$dist<-sapply(1:nrow(e),function(j){s[e[j,"trap1"],e[j,"trap2"]]})
  e
})

x<-rbindlist(l)
x<-x[x$trap1!=x$trap2,]
blocks<-50
x$cut<-cut(x$dist,breaks=seq(-0.01,max(x$dist)+blocks,by=blocks))
xx<-x[,.(n=.N,mean=mean(diff)),by="cut"][order(cut)][1:(5000/blocks),]
xx<-droplevels(xx)
plot(as.integer(xx$cut)+0.5,xx$mean,xaxt="n",type="b",cex=scales::rescale(xx$n,to=c(0.1,5)))
axis(1,at=1:nlevels(xx$cut),label=as.integer(sapply(strsplit(gsub("\\(|\\]","",levels(xx$cut)),","),"[",1)),las=2,cex.axis=0.5)



## Visualize spatial fields ######################################

xlim<-bbox(mappingzone)[1,]
ylim<-bbox(mappingzone)[1,]

proj<-inla.mesh.projector(mesh,xlim=xlim,ylim=ylim,dims=c(300,300))

mfield<-inla.mesh.project(projector=proj,field=m$summary.random[["spatial"]][['mean']])
sdfield<-inla.mesh.project(projector=proj,field=m$summary.random[["spatial"]][['sd']])

par(mfrow=c(1,2),mar=c(3,3,2,5))

image.plot(list(x=proj$x,y=proj$y,z=mfield),col=viridis(100),asp=1,main="Spatial field (log scale)") 
axis(1)
axis(2)
plot(swediv,add=TRUE,border=gray(0,0.5))
# not sure how to represent the q quantile form observations...
#plot(sizesdiv,pch=1,cex=0.1*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.13),add=TRUE)
#brks<-c(0.01,0.25,0.50,0.75,0.99)
#legend("topleft",pch=1,pt.cex=3*brks,col=gray(0,0.3),legend=brks,bty="n",title="Probability of location\nbeing an actual fire",inset=c(0.02,0.05))

image.plot(list(x=proj$x,y=proj$y,z=sdfield),col=viridis(100),asp=1,main="sd of spatial field (log scale)") 
axis(1)
axis(2)
plot(swediv,add=TRUE,border=gray(0,0.5))
#plot(sizesdiv,pch=1,cex=0.1*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.13),add=TRUE)



### Graphical predictions with spatial uncertainty ############################

mm<-glmmTMB(sp~jul+julsquare+forest+urban+tmax1+tmax15,data=xs@data[!is.na(xs$sp),],family=nbinom2())

plot(mesh,asp=1)
plot(mappingzone,add=TRUE)
plot(ds,add=TRUE)
plot(xsmap,add=TRUE)


library(patchwork)
library(ggplot2)
library(ggeffects)


dat<-xs@data[xs$db!="map",]
dat$we<-log(dat$wetland1000+0.001)


mm<-glmmTMB(sp ~ -1 + ns(jul, df=9) + agriculture50 + forest50 + agriculture1000 + forest1000 + anom2 + prcp2 + anom90 + prcp90 + offset(lognights), ziformula=~0, data=dat,family=nbinom2())

vs<-all.vars(formula(mm))[-1]
gl<-lapply(vs,function(i){
  g<-ggeffect(mm,terms=paste(i,"[n=100]"))
  plot(g,add=TRUE,jitter=FALSE)+coord_cartesian(ylim = c(0, 500)) 
})
wrap_plots(gl,nrow=2)


lp<-newdata(x=dat[,vs,drop=FALSE],v=vs,n=100,fun=mean,list=FALSE)
we<-seq(min(dat$we),0,length.out=100)
nd<-cbind(we=we,lp[["we"]][rep(1,length(we)),][,-2])
p<-predict(mm,newdata=nd,type="response")
plot(exp(nd$we),p)


### The following id for checking the stability of results

par(mfrow=c(1,2),oma=c(0,0,0,4))
f<-function(i){log(i)}
zlim<-range(f(c(values(pred[[1]]),values(pr))),na.rm=TRUE)
zlim<-range(f(c(values(pr))),na.rm=TRUE)
xxs<-xs[xs$week%in%xsmap$week[1],]
plot(f(pred[[1]]),zlim=zlim)
plot(st_geometry(water),border=NA,col="white",add=TRUE)
plot(xxs,add=TRUE,pch=1,cex=scales::rescale(xxs$sp,c(0.5,10)))
#plot(resample(pred[[1]],pr))
plot(f(pr),zlim=zlim)
plot(st_geometry(water),border=NA,col="white",add=TRUE)
plot(xxs,add=TRUE,pch=1,cex=scales::rescale(xxs$sp,c(0.5,10)))

#preds<-list()
preds[[length(preds)+1]]<-pred[[1]]

m<-readRDS("C:/Users/God/Downloads/m.rds")
load("C:/Users/God/Downloads/results.RData")
lf<-list.files("C:/Users/God/Downloads",pattern=".rds",full=TRUE)
lf<-sort(lf[grep("pred",lf)])
preds<-lapply(lf,function(i){
  readRDS(i)[[5]]
})


#pred<-readRDS("C:/Users/God/Documents/mosquitos/data/pred.rds")
#preds<-pred[[1]]

#pmesh<-inla.mesh2sp(mesh)$triangles
spreds<-stack(preds)
spreds<-disaggregate(spreds,fact=5,method="bilinear")
#spreds<-log(spreds)
#spreds[spreds>600]<-NA
par(mar=c(0,0,0,0),oma=c(1,1,1,2),mfrow=n2mfrow(nlayers(spreds),asp=2))
lapply(1:nlayers(spreds),function(i){
  plot(spreds[[i]],zlim=range(values(spreds),na.rm=TRUE),axes=FALSE,col=cols)
  plot(st_geometry(water),border=NA,col="white",add=TRUE)
  #plot(xs,add=TRUE,pch=1,cex=scales::rescale(xs$sp,c(0.5,10)))
})
par(mfrow=c(1,1))


par(mar=c(0,0,0,0),oma=c(1,1,1,2))
plot(spreds[[nlayers(spreds)]],zlim=range(values(spreds),na.rm=TRUE))
plot(pmesh,col=NA,border=gray(1),lwd=1,add=TRUE)

# scp C:/Users/God/Documents/mosquitos/models.R root@do:/root

# nohup Rscript models.R --no-save > verbose.out 2>&1 &

# top -o %MEM
