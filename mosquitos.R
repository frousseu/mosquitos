
library(sp)
library(rgdal)
library(scales)
library(gstat)
library(geoR)
library(geostatsp)
library(INLA)
library(rgeos)
library(FRutils)
library(RColorBrewer)
library(visreg)
library(readxl)
library(mgcv)
library(qgam)
library(raster) 
library(rasterVis)
library(data.table)
library(plyr)
library(alphahull)
library(concaveman)
library(mapview)
library(fasterize)
library(sf)
library(velox)
library(glmmTMB)
library(viridis)
library(foreach)
library(doParallel)
library(rmapshaper)

source("https://raw.githubusercontent.com/frousseu/FRutils/master/R/colo.scale.R")
source("https://raw.githubusercontent.com/frousseu/UdeS/master/GStecher/newdata.R")
### this is to tunr the formula using dummy variables for factors
mmatrix<-function(i,dat){
  if(!is.list(i)){
    i<-list(i)  
  }
  lapply(i,function(j){
    mod<-as.character(j)
    mod[3]<-gsub(" \\+ f\\(spatial, model = spde\\)","",mod[3]) # remove spatial effect
    mod[3]<-gsub("-1 \\+ intercept \\+ ","",mod[3]) # remove spatial effect and intercept notation
    mod<-as.formula(paste0(mod[c(1,3)],collapse=" "))
    model.matrix(mod,dat)[,-1]
  })
}
mm<-mmatrix(modell,size)


### this is the model with the dummy variable from the model.matrix as suggested by the pdf E Krainski on using factors
modellmm<-lapply(mm,function(i){
  as.formula(paste("y ~ -1 + intercept +",paste(dimnames(i)[[2]],collapse=" + "),"+ f(spatial, model = spde)"))
})

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

d<-as.data.frame(read_excel("C:/Users/God/Documents/UdeS/Documents/UdeS/Consultation/JAllostry/Doc/BD.xlsx"))
d$Long<-d$LongdecSatScan
d$Lat<-d$LatdecSatScan
d$Long<-d$Long+0.3
d$Lat<-d$Lat-0.2
d$Site<-ifelse(nchar(d$Site)==6,sapply(strsplit(d$Site,"(?<=.{3})", perl = TRUE),paste,collapse=" "),d$Site)
#unique(d$Site)

d$date<-as.Date(d$Day)
d$year<-d$Annee
d$jul<-as.integer(format(d$date,"%j"))


d<-d[order(d$Site,d$Annee_Week),]

l<-split(d,d$year)

ds<-d
coordinates(ds)<-~Long+Lat
proj4string(ds)<-"+init=epsg:4326"

prj<-"+proj=utm +zone=18 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"
ds<-spTransform(ds,CRS(prj))

#mapview(ds,zcol="Site")

l<-list.files("C:/Users/God/Documents/UdeS/Documents/UdeS/Consultation/JAllostry/Doc/lcc2000",pattern=".shp",full.names=TRUE)
s<-lapply(l,st_read)
s<-do.call("rbind",s)
#s<-st_transform(s,crs=26918)
ss<-st_transform(s,crs=st_crs(st_as_sf(ds)))

#https://open.canada.ca/data/en/dataset/97126362-5a85-4fe0-9dc2-915464cfdbb7#
cla<-as.data.frame(read_excel("C:/Users/God/Documents/UdeS/Documents/UdeS/Consultation/JAllostry/Doc/lcc2000class.xlsx"))

clan<-sapply(strsplit(cla[,1],"-"),"[",1)
cla<-strsplit(cla[,1],"-") %>% 
  sapply("[",-1) %>% 
  sapply(paste0,collapse="") %>% 
  sapply(function(i){gsub('[[:digit:]]+', '', i)}) %>% 
  sapply(function(i){gsub("  "," ",i)}) %>% 
  sapply(function(i){gsub("  ","",i)}) %>% 
  unname


cla[grep("Conifères|Coniférien|Mixte|Feuillu|Forêt",cla)]<-"Natural"
cla[grep("Cultures|Prairies|agricoles",cla)]<-"Agriculture"
cla[grep("humide",cla)]<-"Wet"
cla[grep("Arbustes|arbustes|herbacées",cla)]<-"Natural"
cla[grep("découvert|Stérile|Bryo",cla)]<-"Exposed"
cla[grep("développées",cla)]<-"Urban"
cla[grep("Ombre",cla)]<-NA

s$cover<-factor(cla[match(s$COVTYPE,clan)])

lcc<-levels(s$cover)

r<-raster(nrow=2000,ncol=3000,ext=extent(s))

rx<-fasterize(s,r,field="cover",fun="max")
rn<-fasterize(s,r,field="cover",fun="min")

rlcc<-ratify(rx)

levels(s$cover)
hn<-levels(rlcc)[[1]][,1]
cols<-c("lightgoldenrod","skyblue","grey30","forestgreen","grey60","brown")
#cols<-sample(colors(),length(hn))
par(mar=c(0,0,0,0))
plot(rlcc,legend=FALSE,col=cols)
legend("bottomright",legend=levels(s$cover),fill=cols,bty="n",border=NA,cex=2)
points(ds,col="red",pch=16,cex=0.8)
dbuffer<-1
buffer<-st_transform(st_buffer(st_as_sf(ds),dist=dbuffer),4326)
plot(st_geometry(buffer),add=TRUE)

vrlcc<-velox(rlcc)
l<-vrlcc$extract(buffer)
l<-lapply(l,function(i){
  x<-i[,1]
  x<-x[!is.na(x)]
  sum(x==4L)/length(x)
})

ds$natural<-unlist(l)


grid<-st_make_grid(rlcc,cellsize=c(0.0200,0.0200),what="centers")
pr<-raster(SpatialGrid((points2grid(as_Spatial(grid)))))
g<-st_buffer(st_as_sf(xyFromCell(pr,1:ncell(pr),spatial=TRUE)),dist=0.01)
l<-v$extract(g)

l<-lapply(l,function(i){
  x<-i[,1]
  x<-x[!is.na(x)]
  sum(x==4L)/length(x)
})
pr<-setValues(pr,unlist(l))
plot(pr,col=rev(viridis(100)))

#tmap_mode("view")
#tm_shape(ds) +
#  tm_dots() +
#  tm_layout(basemaps = c("Esri.WorldImagery", "Esri.WorldShadedRelief", "Esri.NatGeoWorldMap"))

can<-raster::getData("GADM", country = "CAN", level = 2)
que<-can[can$NAME_1=="Québec",]

Q<-as(ms_simplify(st_as_sf(spTransform(que,CRS(proj4string(ds)))),keep=0.03),"Spatial")

plot(ds)
plot(Q,add=TRUE)
l<-locator()
h<-gConvexHull(SpatialPoints(cbind(l$x,l$y),proj4string=CRS(proj4string(ds))))
plot(h,add=TRUE)
o<-over(ds,h)
d<-d[!is.na(o),]
ds<-ds[!is.na(o),]

#plot(d$jul,log(d$sp+1))

#m<-gam(sp~s(jul),data=d,family=nb())


par(mar=c(0,0,0,0))
plot(st_geometry(st_as_sfc(st_bbox(buffer))))
plot(rlcc,legend=FALSE,col=cols,add=TRUE)
#points(ds,col="red",pch=16,cex=0.8)
plot(st_geometry(buffer),add=TRUE)
legend("bottomright",legend=levels(s$cover),fill=cols,bty="n",cex=2,bg=gray(1,0.5))



rm(rx,rn,s,can);gc();gc()

# "http://forobs.jrc.ec.europa.eu/products/glc2000/legend/GLC2000_Lccs_110604_export.xls"
#lcc<-as.data.frame(read_xls("C:/Users/rouf1703/Downloads/GLC2000_Lccs_110604_export.xls"))
#cov<-raster("C:/Users/rouf1703/Downloads/CAN_cov/CAN_cov.grd")
#cov<-crop(cov,extent(-76,-70,45,49))
#cov<-crop(cov,xs2)
#que<-spTransform(que,proj4string(ds))
#cov<-projectRaster(cov,crs=CRS(proj4string(ds)))
#mapview(cov)+mapview(ds[1:1000,])
#plot(ds)
#plot(que,add=TRUE)
#plot(cov,add=TRUE)

######################
### weeks
library(reshape)
library(tidyr)
par(mar=c(0,0,0,3))
d$wy<-paste(d$year,d$Week,sep="_")
d$k<-1
x<-reshape(unique(d[,c("Site","wy","k")]), idvar="Site", timevar="wy", direction="wide")
x<-x[,c(names(x)[1],sort(names(x)[-1]))]
x<-raster(t(as.matrix(x[,-1])))
plot(x,asp=NULL)

################################################
### space-time simple from spde tutorial

colSums(xs@data[,5:33])
#xs1<-ds[ds$year%in%c("2015"),]
#xs2<-ds[ds$year%in%c("2016"),]
#xs1$sp<-log(xs1$A9+1)
#xs2$sp<-log(xs2$A9+1)
xs<-ds[ds$year%in%c("2015"),]
#xs$sp<-log(xs$A29+1)
xs$sp<-xs$A3
#xs<-rbind(xs1,xs2)
#xs<-xs1
#xs$Week<-sample(xs$Week)
#xs$sp<-xs$A29
xs<-xs[order(xs$Annee,xs$Week),]
xs$Week<-paste(xs$Annee,xs$Week,sep="_")
xs$week<-as.integer(substr(xs$Week,6,7))
xs$week<-xs$week
xs$week2<-xs$week^2

domain <- inla.nonconvex.hull(coordinates(xs),convex=-0.075, resolution = c(100, 100))
mesh<-inla.mesh.2d(loc=coordinates(xs),max.edge=c(5,10),offset=c(5,5),cutoff=5,boundary=domain,crs=CRS(proj4string(xs)))
plot(mesh,asp=1)
plot(xs,add=TRUE,pch=1,col="red")

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  prior.range=c(5, 0.01), ### P(practic.range<0.05)=0.01
  prior.sigma=c(4, 0.5)) ### P(sigma>1)=0.01

## ----rfindex-------------------------------------------------------------
k<-length(unique(xs$Week))
iset <- inla.spde.make.index('i', n.spde=spde$n.spde, n.group=k)

## ----apred---------------------------------------------------------------
A <- inla.spde.make.A(mesh=mesh, 
                      loc=coordinates(xs), 
                      group=as.integer(factor(xs$Week))) 

## ----stack---------------------------------------------------------------
sdat<-inla.stack(tag='stdata',data=list(y=xs$sp),A=list(A,1),effects=list(c(iset,list(intercept=1)),data.frame(w=xs$Typ,week=xs$week,week2=xs$week2,natural=xs$natural))) 

## ----hbeta---------------------------------------------------------------

co<-seq(-0.99,0.99,by=0.01)
u<-0
alpha<-0.00005
par(mfrow=c(1,3))
hist(inla.pc.rcor0(10000,u=mu,alpha=alpha))
dens<-inla.pc.dcor0(co,u=mu,alpha=alpha)
plot(co,dens,type="l",ylim=c(0,max(dens,na.rm=TRUE)))
sig<-inla.pc.rprec(10000,u=1,alpha=0.05)
hist(1/sig)
table(sig>1)


h.spec <- list(theta=list(prior="pc.cor0", param=c(0.1, NA)))
#h.spec <- list(theta=list(prior="pc.prec", param=c(0.5,0.5)), rho=list(prior="pc.cor1", param=c(0.9,0.9)))
#h.spec <- list(theta = list(prior="pc.prec", param=c(1, NA)),
#               rho = list(prior="pc.cor0", param=c(0.1, NA)))

h.spec <- list(#theta=list(prior='pc.prec', param=c(0.00000000001, 0.00000000005)))#,
  rho = list(prior="pc.cor0", param=c(0.5,0.1)))

#h.spec <- list(theta = list(prior = "betacorrelation",param=c(1,3),initial=-1.098))
#hist(rbeta(10000,1,3))

## ----remote,echo=FALSE---------------------------------------------------
##inla.setOption(inla.call='remote')

## ----ft------------------------------------------------------------------
formulae <- y ~ -1 + intercept + week + week2 + natural + f(i, model=spde, group=i.group,control.group=list(model='ar1', hyper=h.spec)) 
#formulae <- y ~ 0 + w + f(i, model=spde) + f(week,model="rw1")
#formulae <- y ~ 0 + w + f(i, model=spde, group=i.group,control.group=list(model='exchangeable')) 
prec.prior <- list(prior='pc.prec', param=c(1, 0.05))
m <- inla(formulae,  data=inla.stack.data(sdat), 
          control.predictor=list(compute=TRUE, A=inla.stack.A(sdat),link=1), 
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=list(expand.factor.strategy='inla'),
          control.inla = list(int.strategy = "eb"),
          num.threads=3,
          verbose=FALSE,
          family="nbinomial")#"zeroinflatednbinomial1"

## ----sbeta---------------------------------------------------------------
#round(cbind(observed=tapply(dat$y, dat$w, mean), m$summary.fixed), 4) 

## ----echo=FALSE,fig.width=5.5,fig.height=5.5-----------------------------
par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=2:0)
for (j in 1:4) {
  plot(m$marginals.hyper[[j]], type='l', 
       xlab=names(m$marginals.hyper)[j], ylab='Density')
  abline(v=c(1/sd.y^2, sqrt(8)/params[1], 
             params[2]^0.5, rho)[j], col=2)
}

## ----rfidx---------------------------------------------------------------
#str(idat <- inla.stack.index(sdat, 'stdata')$data) 

## ----meanrf--------------------------------------------------------------
#cor(xs$sp, m$summary.linear.predictor$mean[idat])

## ----projgrid------------------------------------------------------------
stepsize <- 5*1/1
coords<-st_coordinates(st_cast(st_buffer(st_as_sf(xs),10),"MULTIPOINT"))
nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
projgrid <- inla.mesh.projector(mesh, xlim=range(coords[,1]),ylim=range(coords[,2]), dims=nxy,crs=CRS(proj4string(xs)))

## ----projpmean-----------------------------------------------------------
xmean <- list()
for (j in 1:k){
  xmean[[j]] <- inla.mesh.project(projgrid,m$summary.random$i$mean[iset$i.group==j])
}

## ----inout---------------------------------------------------------------
b<-gBuffer(gConvexHull(SpatialPoints(domain$loc,p=CRS(proj4string(ds)))),width=0.1,byid=FALSE)
a<-concaveman(coordinates(xs),2)
a<-gBuffer(spPolygons(a,crs=CRS(proj4string(xs))),width=5)
o <- over(SpatialPoints(projgrid$lattice$loc,p=CRS(proj4string(ds))),b)

## ----setNAs---------------------------------------------------------------
for (j in 1:k)   xmean[[j]][is.na(o)] <- NA
r<-stack(lapply(xmean,function(i){
  raster(nrows=nxy[2], ncols=nxy[1], xmn=min(projgrid$x), xmx=max(projgrid$x), ymn=min(projgrid$y), ymx=max(projgrid$y),crs=CRS(proj4string(xs)),vals=as.vector(i[,ncol(i):1])) ## some crazy ordering in INLA output be careful
  #raster(i)
}))
names(r)<-unique(xs$Week)
cols<-colo.scale(200,c("steelblue3","orange","red3","darkred","grey20"))

# voir argument panel.number ou packets de layer
xsbuff<-st_coordinates(st_cast(st_buffer(st_as_sf(xs),7),"MULTIPOINT"))[,1:2]
buf<-concaveman(xsbuff,10)
buf<-spPolygons(buf,crs=CRS(proj4string(xs)))
buf<-gBuffer(buf,width=1)
r<-mask(r,buf)

xxs<-split(xs,xs$Week)
p.strip<-list(cex=0.65,lines=1,col="black")
levelplot(r,col.regions=cols,cuts=199,par.strip.text=p.strip,par.settings = list(axis.line = list(col = "grey90"),strip.background = list(col = 'transparent'),strip.border = list(col = 'grey90')),scales = list(col = "black")) +
  layer(sp.points(xxs[[panel.number()]],col=gray(0,0.5),pch=1,cex=scales:::rescale(c(max(xs$sp),identity(xxs[[panel.number()]]$sp)),to=c(0.2,3))[-1]))+
  layer(sp.points(xxs[[panel.number()]],col=gray(0,0.5),pch=3,cex=0.25))+
  layer(sp.polygons(Q,col=gray(0,0.3)))
par(mfrow=c(1,1))

#plot(m$summary.fitted.values$mean[1:nrow(xs)],xs$sp)

plot(m$marginals.hyperpar$`GroupRho for i`,type="l",xlim=c(0,1))
vv<-seq(0.01,0.99,by=0.01)
lines(vv,dbeta(vv,1,5))

#levelplot(exp(r)-1,zscaleLog=10)

x<-seq(min(xs$week),max(xs$week),by=0.1)
y<-exp(m$summary.fixed[1,1]+m$summary.fixed[2,1]*x+m$summary.fixed[3,1]*x^2)
#xx<-(20:40-mean(xs$week))/sd(xs$week)
xx<-20:40
plot(x,y,type="l",xaxt="n",xlim=range(xx))
axis(1,at=xx,label=20:40)


###################################################################
### SHOW PREDICTONS
###################################################################

###################################################################
### build prediction grid
###################################################################

plot(st_geometry(st_transform(st_as_sf(inla.mesh2sp(mesh)$triangles),4326)))

pr<-raster(ext=extent(a),res=c(2,2),crs=CRS(proj4string(a)))
g<-st_buffer(st_as_sf(xyFromCell(pr,1:ncell(pr),spatial=TRUE)),dist=1)

plot(mesh,asp=1)
plot(g,add=TRUE)

plot(rlcc)
plot(st_transform(g,4326),add=TRUE)

l<-vrlcc$extract(st_transform(g,4326))
l<-lapply(l,function(i){
  x<-i[,1]
  x<-x[!is.na(x)]
  x<-factor(lcc[x],levels=lcc)
  table(x)/sum(table(x))
})
pgrid<-as.data.frame(do.call("rbind",l))
ll<-lapply(pgrid,function(i){
  setValues(pr,i)
})
rgrid<-stack(ll)




###################################################################
### build prediction matrices for the map and the prediction graphs
weeks<-unique(as.integer(factor(xs$Week)))
Ap<-inla.spde.make.A(mesh=mesh,loc=coordinates(rgrid)[rep(1:ncell(rgrid),length(weeks)),][1,,drop=FALSE],group=rep(weeks,each=ncell(rgrid))[1])
n<-5 # number of divisions in generated values for the focus variable
Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(600,5050),ncol=2)[rep(1,n),,drop=FALSE],group=rep(12,n))

# set location for predictions with the values of variables at this location?

################################################
### build newdata with variable values to submit
v<-setdiff(all.vars(formulae),c("y","i","intercept","spde","i.group","h.spec"))
v2<-v[grep("2",v)]
v1<-setdiff(v,v2)
lp<-newdata(x=xs@data[,v1,drop=FALSE],v=v1,n=n,fun=median,list=FALSE)
if(length(v2)){
  lp<-lapply(lp,function(i){
    a<-as.data.frame(lapply(i[,gsub("2","",v2),drop=FALSE],"^",2))
    names(a)<-v2
    cbind(i,a)
  })
}
model<-formula(paste("y~-1+",paste(v,collapse="+")))
#lp<-lapply(lp,function(i){mmatrix(model,i)})
#lpmed<-mmatrix(model,newdata(x=xs[,v,drop=FALSE],v=v,n=1,fun=median,list=FALSE)[[1]][rep(1,length(g)),])[[1]]
vmap<-data.frame(week=sort(unique(xs$week)))
vmap$week2<-vmap$week^2
vmap<-vmap[rep(1:length(weeks),each=ncell(rgrid)),][1,]

vpred<-data.frame(week=seq(20,40,length.out=n))
vpred$week2<-vpred$week^2

########################################################
### bind the data stack for the estimate and for the map
stack.est<-inla.stack(data=list(y=xs$sp),A=list(A,1),effects=list(c(iset,list(intercept=1)),data.frame(xs@data[,v])),tag="est")
stack.map<-inla.stack(data=list(y=NA),A=list(Ap,1),effects=list(c(lapply(iset,"[",1),list(intercept=1)),vmap),tag="map")
stack.pred<-inla.stack(data=list(y=NA),A=list(Apn,1),effects=list(c(lapply(iset,"[",iset$i.group==12),list(intercept=1)),vpred),tag="pred")
full.stack<-inla.stack(stack.est,stack.map,stack.pred)
#full.stack<-inla.stack(stack.est,stack.pred)

#######################################
### add a stack for each focus variable
for(i in seq_along(v1)){
  le<-nrow(lp[[v1[i]]])
  if(le!=n){
    #AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,le),,drop=FALSE]) # for categorical variables
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(600,5050),ncol=2)[rep(1,n),,drop=FALSE],group=rep(12,n))
  }else{
    #AA<-Apn # for numerical variables
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(600,5050),ncol=2)[rep(1,n),,drop=FALSE],group=rep(12,n))
  }
  stack<-inla.stack(data=list(y=NA),A=list(AA,1),effects=list(c(lapply(iset,"[",iset$i.group==12),list(intercept=1)),lp[[v1[i]]]),tag=v1[i])     
  full.stack<-inla.stack(full.stack,stack)
}

##########################################
### extract index for each stack
index.est<-inla.stack.index(full.stack,tag="est")$data
index.map<-inla.stack.index(full.stack,tag="map")$data
#index.pred<-inla.stack.index(full.stack,tag="pred")$data
#index<-list(est=index.est,map=index.map,pred=index.pred)
index<-list(est=index.est,map=index.map)
#index<-list(est=index.est,pred=index.pred)
for(i in seq_along(v1)){
  index<-c(index,list(inla.stack.index(full.stack,tag=v1[i])$data))
}  
names(index)[3:length(index)]<-v1

##################################################
### rerun best model with each variable to predict

#m<-inla(modellmm[[b]],data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.inla=list(strategy='simplified.laplace',int.strategy="eb"),family="gp",control.family=list(list(control.link=list(quantile=q),hyper=hyper.gp)),control.fixed=control.fixed,num.threads=7)

m <- inla(formulae,  
          data=inla.stack.data(full.stack), 
          control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),
          control.compute=list(dic=FALSE,waic=FALSE,cpo=FALSE,config=TRUE),
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=list(expand.factor.strategy='inla'),
          control.inla=list(int.strategy="eb"),#list(strategy='simplified.laplace',int.strategy="eb"),
          num.threads=3,
          verbose=FALSE,
          family="nbinomial")#"zeroinflatednbinomial1"

### from haakon bakka, BTopic112
nsims<-500
samples<-inla.posterior.sample(nsims,m,num.threads="6.6")
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor" # not sure if should use APredictor or Predictor
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-do.call("cbind",samples.effect)
xi.eff<-sapply(samples, function(x) x$hyperpar[grep("Rho",names(x$hyperpar))])
#unique(sapply(strsplit(rownames(m$summary.fitted.values),"\\."),function(i){paste(i[1:min(2,length(i))],collapse=" ")}))

#####################################
### visualize spatial fields

xlim<-range(coordinates(r)[,1])
ylim<-range(coordinates(r)[,2])

proj<-inla.mesh.projector(mesh,xlim=xlim,ylim=ylim,dims=c(300,300))

mfield<-inla.mesh.project(projector=proj,field=m$summary.random[['i']][['mean']])
sdfield<-inla.mesh.project(projector=proj,field=m$summary.random[['i']][['sd']])

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


##################################
### build a relative frequency map

p<-m$summary.fitted.values[index[["map"]],"mean"] # the lambda is to back-transform on the original scale)
rgp<-setValues(r,p)
rr<-mask(rgp,buf)
brks <- seq(min(p),max(p),by=0.1)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(rgp,col=cols,axes=FALSE,legend.shrink=1, legend.width=4,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),legend.args=list(text='Abundance', side=4, font=2, line=2.3))
plot(swediv,add=TRUE,border=gray(0,0.25),lwd=0.01)

### sd
p<-m$summary.fitted.values[index[["map"]],"sd"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swediv,axes=TRUE)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='sd of predicted fire size', side=4, font=2, line=2.3))
plot(swediv,add=TRUE,border=gray(0,0.25),lwd=0.01)


####################################################
### graphical predictions with spatial uncertainty

### this section is not that useful because it is a prediction for a given location, hence it includes uncertainty in the spatial field

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
for(i in seq_along(v1)){
  p<-m$summary.fitted.values[index[[v1[i]]],c("0.025quant","0.5quant","0.975quant")]
  #p[]<-lapply(p,transI)
  dat<-data.frame(lp[[v1[i]]])
  if(nrow(p)==n){
    plot(dat[[v1[i]]],p[,2],type="l",ylim=c(0,max(c(max(p[,3]),max(xs$sp)))),xlab=v1[i],font=2,ylab="",lty=1,yaxt="n")
    lines(dat[[v1[i]]],p[,1],lty=3,lwd=1)
    lines(dat[[v1[i]]],p[,3],lty=3,lwd=1)
    points(xs@data[,v1[i]],xs$sp,pch=16,col=gray(0,0.07))
  }else{
    plot(unique(sort(size[,v1[i]])),p[,2],type="l",ylim=c(0,100),xlab=v1[i],font=2,ylab="",lty=1,yaxt="n")
    segments(x0=as.integer(unique(sort(size[,v[i]]))),x1=as.integer(unique(sort(size[,v[i]]))),y0=p[,1],y1=p[,3],lty=3,lwd=2)
    points(jitter(as.integer(size[,v[i]]),fac=2),transI(size$tTotal),pch=16,col=gray(0,0.07))
  }
  axis(2,las=2)
}
mtext("Fire size in ha",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)


##########################################################
### generate predictions without spatial uncertainty

# page 263 in Zuur

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  grep(paste0(i,":"),row.names(samples[[1]]$latent))  
}) 
nweights<-grep("i",row.names(samples[[1]]$latent))

### this is to compare with a quantile model
#par(mfrow=c(round(sqrt(2*length(v)),0),ceiling(sqrt(2*length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
#mq <- rq (tTotal ~ Road1k + Pp_1000 + urbwtr1k + Ag_1000 + h__1000 + FWI, data = na.omit(size),tau=q)
#na<-all.vars(mq$formula[[3]])
#grobs<-lapply(na,function(i){
#  if(is.factor(size[,i])){
#    plot(ggpredict(mq,terms=i),limits=c(0,100),raw=TRUE)  
#  }else{
#    plot(ggpredict(mq,terms=paste(i,"[n=50]")),limits=c(0,100),raw=TRUE) 
#  }
#})
#grid.arrange(grobs=grobs,ncol=3)

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,3,2,2),oma=c(0,10,0,0))
for(k in seq_along(v1)){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    fixed<-cbind(intercept=1,as.matrix(lp[[v1[k]]])) %*% betas
    #fixed<-cbind(intercept=1,as.matrix(dat)) %*% betas
    ### this if we want a spatial part
    #wk<-samples[[i]]$latent[nweights]
    #if(is.factor(size[,v[k]])){
    #  spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    #  spatial<-as.matrix(Apn) %*% wk
    #}
    #p<-exp(fixed+spatial)
    p<-exp(fixed)
    p
  })
  p<-do.call("cbind",p)
  p<-t(apply(p,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  if(nrow(lp[[v1[k]]])==n){
    vals<-lp[[v1[k]]][,v1[k]]
    plot(vals,p[,2],type="l",ylim=c(0,100),xlab=v1[k],font=2,ylab="",lty=1,yaxt="n")
    points(xs@data[,v1[k]],xs$sp,pch=1,col=gray(0,0.15))
    lines(vals,p[,2],lwd=2)
    lines(vals,p[,1],lty=3)
    lines(vals,p[,3],lty=3)
  }else{
    plot(unique(sort(size[,v[k]])),p[,2],type="l",ylim=c(0,100),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    points(jitter(as.integer(size[,v[k]])),size$tTotal,pch=1,col=gray(0,0.15))
    segments(x0=as.integer(unique(sort(size[,v[k]]))),x1=as.integer(unique(sort(size[,v[k]]))),y0=p[,1],y1=p[,3],lty=3)
  }
  axis(2,las=2)
}
mtext(paste("Fire size at the",q,"quantile (ha)"),outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)

###############################################
### range and sigma

### this is to show the posteriors of the spatial field

res <- inla.spde.result(m, "spatial", spde)
par(mfrow=c(2,1))
plot(res$marginals.range.nominal[[1]],
     type="l", main="Posterior density for range")
plot(inla.tmarginal(sqrt, res$marginals.variance.nominal[[1]]),
     type="l", main="Posterior density for std.dev.")
par(mfrow=c(1,1))





