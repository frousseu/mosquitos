
## MODELS #####################################################

load("data.RData")

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
  SMG=list(
    ~ agriculture50 + forest50,
    ~ agriculture1000+ forest1000,
    ~ agriculture50 + forest50 + agriculture1000+ forest1000
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
spcode<-c("VEX_","CPR_","CQP_","SMG_")[3] # change index to change species
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

#vs<-unique(unlist(lapply(unlist(models),all.vars)))
vs<-names(xs)[grep("jul|50|1000|anom|prcp|tmean",names(xs))]
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
  
  set.seed(1234)
  xs3pts<-as(st_sample(st_as_sf(mappingzone),10000),"Spatial")
  xs2pts<-as(st_cast(predmap,"MULTIPOINT"),"Spatial")
  xs2<-as(xs2,"Spatial")
  xs2map<-as(xs2map,"Spatial")
  
  plot(st_geometry(predmap),add=TRUE,border="red")
  
  plot(xs2pts,add=TRUE)
  plot(xs3pts,add=TRUE)
  
}


edge<-0.5
#domain <- inla.nonconvex.hull(coordinates(ds),convex=-0.015, resolution = c(100, 100))
#mesh<-inla.mesh.2d(loc.domain=coordinates(ds),max.edge=c(edge,3*edge),offset=c(edge,1*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
#domain <- inla.nonconvex.hull(coordinates(xs2pts),convex = -0.15, concave = 0.5, resolution = c(340,340))
domain <- inla.nonconvex.hull(rbind(coordinates(xs2pts),coordinates(xs3pts)),convex = -0.015,resolution=75)
#domain <- inla.nonconvex.hull(rbind(coordinates(xs2pts)),convex = -0.12)
#ims<-inla.mesh.segment(loc=coordinates(xs2pts))
mesh<-inla.mesh.2d(loc.domain=NULL,max.edge=c(edge,3*edge),offset=c(edge,3*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
#mesh<-inla.mesh.2d(loc=coordinates(xs),max.edge=c(2,8),cutoff=2)
#mesh<-inla.mesh.2d(boundary=domain,max.edge=c(edge,2*edge),offset=NULL,cutoff=0.5*edge,crs=CRS(proj4string(xs)))
plot(mesh,asp=1)
plot(Q,col="grey90",border="white",add=TRUE,lwd=2)
plot(mesh,asp=1,add=TRUE)
plot(xs[!xs$db%in%"map",],add=TRUE,pch=16,col="forestgreen")
plot(mappingzone,add=TRUE)
mtext(side=3,line=-2,text=paste("edge =",edge,"km"),font=2,cex=1.5)


#### Restrict predictions ######################################
xsmap<-xs[xs$db=="map",]
xs<-xs[xs$db!="map",]

#### SPDE #################################################
spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  constr = FALSE, # not exactly sure what this does
  prior.range=c(5, 0.1), ### P(practic.range<0.05)=0.01
  prior.sigma=c(0.5,0.1)) ### P(sigma>1)=0.01

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
knots<-seq(min(xs$jul)+0.25,max(xs$jul)-0.25,length.out=9)

v<-setdiff(unique(unlist(lapply(unlist(models),all.vars))),c("y","spatial","intercept","spde","year","knots","lognights")) 

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


#### Index #################################################
index.est<-inla.stack.index(stackfull,tag="est")$data
index.map<-inla.stack.index(stackfull,tag="map")$data
index<-list(est=index.est,map=index.map)
for(i in seq_along(v1)){
  index<-c(index,list(inla.stack.index(stackfull,tag=v1[i])$data))
}  
names(index)[3:length(index)]<-v1

save.image(paste0(spcode,"model_parameters.RData"))

