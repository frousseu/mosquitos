
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
load("data.RData")


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

#save("data.RData")
