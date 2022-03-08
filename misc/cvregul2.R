

### Packages #########################################
#library(sp)
#library(rgdal)
#library(scales)
library(INLA)
#library(rgeos)
#library(FRutils)
#library(readxl)
#library(raster) 
#library(rasterVis)
#library(data.table)
#library(alphahull)
#library(concaveman)
#library(mapview)
#library(fasterize)
#library(sf)
#library(velox) 
#library(viridis)
#library(foreach)
#library(doParallel)
#library(rmapshaper)
#library(ncdf4)
#library(daymetr)
#library(gdalUtils)
#library(rgdal)
#library(exactextractr)
library(foreach)
library(doParallel)
#library(DHARMa)
#library(corrplot)
#library(abind)
#library(mgcv)
library(surveillance)
#library(doFuture)

registerDoParallel(detectCores()-0) 
getDoParWorkers()
#registerDoMC(cores=4)
#registerDoFuture()
#cl <- makeCluster(4)
#plan(cluster, workers = cl)

load("mosquitos.RData")

### Czado et al. 2009 
#n<-1000
#x<-rep(5,n)
#mu<-seq(0,100,length.out=n)
#size=0.5
rpsmc<-function(x,mu,size=NULL,n=10000){
  sapply(seq_along(mu),function(i){
    X<-rnbinom(n,mu=mu[i],size=size)
    Xp<-rnbinom(n,mu=mu[i],size=size)
    mean(abs(X-x[i]))-((1/2)*(mean(abs(X-Xp))))
  })
}
#rpsmc(x=x,mu=mu,size=size,n=10000)



inla.setOption(inla.mode="experimental")
year<-c(2013:2016)#c(2003:2006,2013:2016);
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
cat(paste(sort(unique(xs$week)),collapse=", "))
xs<-xs[xs$db!="map",]

vals<-list(intercept=1/30^2,jul=1/30^2,julsquare=1/30^2,default=1/30^2) #5-30
control.fixed<-list(prec=vals,mean=list(intercept=0,default=0),expand.factor.strategy = "inla")

vars<-names(ds)[grep("jul|tmax|prcp|anom|forest|agriculture|water|urban|pond|swamp|pasture|crop|wet|barren",names(ds))]
vars<-vars[-grep("CQ|tmax90|tmax30",vars)]
model<-formula(paste0("y ~ -1 + intercept +",paste(vars,collapse="+"),"+ f(spatial,model=spde,group=spatial.group,control.group=list(model=\"ar1\", hyper=h.spec))"))
#model <- sp ~ jul + julsquare + forest50 + urban50 + urban1000 + agriculture1000  + tmax7 + tmax2 + prcp30

#### Mesh #####################################################
edge<-10
domain <- inla.nonconvex.hull(coordinates(ds),convex=-0.015, resolution = c(100, 100))
mesh<-inla.mesh.2d(loc.domain=coordinates(ds),max.edge=c(edge,3*edge),offset=c(edge,2*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))


#sigma<-sort(c(0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.5,1,5,10,20))
sigma<-sort(c(0.005,0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.5,1,2,5))
sigma<-rep(0.05,length(sigma))
#sigma<-sort(c(0.000001,0.001))
#sigma<-sigma[seq(1,length(sigma),by=2)]
ssigma<-c(0.00001,0.00005,0.0001,0.00025,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1)
#ssigma<-seq(0.01,0.10,length.out=length(sigma))
#ssigma<-seq(0.05,0.2,length.out=length(sigma))
#ssigmap<-c(NA,rep(0.01,length(sigma)-1))
#sigma<-c(0.005,0.02,0.05,0.1,0.2,0.5)
#ssigma<-c(0.005,0.01,0.025,0.05,0.1,0.5)
#sigma<-sigma[seq(1,length(sigma),by=2)]
#ssigma<-ssigma[seq(1,length(ssigma),by=2)]

#sigma<-sigma[seq(1,length(sigma),by=2)]
#cv<-c(2003:2005,2013:2015)
cv<-c("Nord","Sud","Laval","Montréal")
#cv<-c("Montréal")
ex<-expand.grid(1:length(sigma),1:length(cv))
iter<-data.frame(sigma=sigma[ex[,1]],ssigma=ssigma[ex[,1]],cv=cv[ex[,2]])
iter<-iter[order(iter$sigma,iter$cv),]
cf<-control.fixed
#cvsco<-vector(mode="list",length=length(sigma))
#sco<-vector(mode="list",length=nrow(iter))
#co<-vector(mode="list",length=nrow(iter))

#mcoptions <- list(preschedule = FALSE)
#opts <- list(preschedule=FALSE)
#opts <- list()
cvsco<-foreach(i=1:nrow(iter),.packages=c("INLA","surveillance"),.options.multicore=list(preschedule=FALSE)) %dopar% {

  cf$prec$default<-1/(iter$sigma[i]^2)
  
  #### SPDE #################################################
  spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  constr = FALSE, # not exactly sure what this does
  #prior.range=c(5000, NA), ### P(practic.range<0.05)=0.01
  prior.range=c(5, 0.01), ### P(practic.range<0.05)=0.01
  #prior.sigma=c(0.000001, NA)) ### P(sigma>1)=0.01
  prior.sigma=c(iter$ssigma[i], 0.01)) ### P(sigma>1)=0.01
  
  
  print(paste(iter$sigma[i],"-",iter$cv[i],"Start"))
  xxs<-xs
  xxs@data[vars]<-lapply(xxs@data[vars],scale)
  xxs$julsquare<-xxs$jul^2
  obs<-xxs$sp[xxs$region%in%iter$cv[i]]
  xxs$sp[xxs$region%in%iter$cv[i]]<-NA
    
  #plot(xxs,pch=16,col="white")
  #plot(xxs[is.na(xxs$sp),],pch=16,col="dodgerblue",add=TRUE)
  #plot(xxs[!is.na(xxs$sp),],pch=16,col="forestgreen",add=TRUE)
  #plot(Q,add=TRUE)

    
  #### Priors on hyperpar ##################################
  h.spec <- list(#theta=list(prior='pc.prec', param=c(0.5, 0.5)))#,
    #rho = list(prior="pc.cor0", param=c(0.7,0.3)))
    rho = list(prior="pc.cor1", param=c(0.9,0.25))) #0.9 0.25
  #h.spec <- list(theta = list(prior = "betacorrelation",param=c(1,3),initial=-1.098))
  #hist(rbeta(10000,1,3))
    
  #### Make index ##########################################
  k<-length(unique(xxs$week))
  iset<-inla.spde.make.index("spatial",n.spde=spde$n.spde,n.group=k)
    
  #### A matrix ##############################################
  Aest<-inla.spde.make.A(mesh=mesh,loc=coordinates(xxs),group=as.integer(factor(xxs$week))) 
    
  #### Stacks ################################################
  stackest<-inla.stack(tag='est',data=list(y=xxs$sp),A=list(Aest,1),effects=list(c(iset,list(intercept=1)),xxs@data)) 
    
  #### Index #################################################
  index.est<-inla.stack.index(stackest,tag="est")$data

  #### Model ##################################################
  m <- inla(model,data=inla.stack.data(stackest), 
              control.predictor=list(compute=TRUE, A=inla.stack.A(stackest),link=1), 
              control.fixed=cf,
              control.inla = list(strategy='gaussian',int.strategy = "eb"),
              num.threads="1:1",
              verbose=TRUE,
              control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
              family="nbinomial")#"zeroinflatednbinomial1"
    
  pred<-m$summary.fitted.values[index.est,][xxs$region%in%iter$cv[i],1]
  size<-m$summary.hyperpar[grep("size",row.names(m$summary.hyperpar)),1]
  
  sco<-scores(obs,pred,size=size,which=c("logs","dss","ses","rps"),rpsmc(obs,pred,size=size,n=10000))  
  colnames(sco)[5]<-"rpsmc"
  sco<-colSums(sco)
  co<-cbind(t(m$summary.fixed[,1,drop=FALSE]),t(m$summary.hyper[,1,drop=FALSE]))
    
  #print(m$summary.fixed[,1:5])
  #print(m$summary.hyperpar[,1:5])
  print(range(pred))
  print(paste(iter$sigma[i],"-",iter$cv[i],"Done"))
    
  
  #cvsco[[i]]<-colSums(do.call("rbind",sco))
  #list(colSums(do.call("rbind",sco)),colMeans(do.call("rbind",co)),sigma,ssigma)
  list(sco,co)
  
}
#sigma<-sigma[c(1:13)]

cvsco<-list(scores=cbind(iter,do.call("rbind",lapply(cvsco,"[[",1))),coefs=cbind(iter,do.call("rbind",lapply(cvsco,"[[",2))))

saveRDS(cvsco,"cvsco.rds")




if(FALSE){
  
cvsco<-readRDS("C:/Users/God/Documents/mosquitos/cvsco.rds")
res<-aggregate(.~sigma+ssigma,data=cvsco$scores[,-3],FUN=sum)
res[3:ncol(res)]<-lapply(res[3:ncol(res)],scale)
cos<-aggregate(.~sigma+ssigma,data=cvsco$coef[,-3],FUN=mean)
cos$"Range for spatial"<-cos$"Range for spatial"/20
cos<-cos[,!names(cos)%in%c("intercept","Range for spatial")]
#sigma<-lapply(cvsco,"[[",3)
#ssigma<-lapply(cvsco,"[[",4)
#cos<-cos[-(1:2),]

#matplot(res[,-1],type="b",lwd=2,cex=2,xaxt="n")
#axis(1,at=seq_along(res$sigma),labels=res$sigma,las=2)
#legend("top",legend=names(res)[-1],lwd=2,col=1:ncol(res[,-1]))

par(mfrow=c(1,2),mar=c(4,3,0.5,0.5))
param<-"ssigma"
f<-function(i){identity(i)}
plot(0.1,0.1,xlim=range(f(res[,param])),ylim=range(res[,-(1:2)]),type="n",xaxt="n",log="x")
lapply(3:ncol(res),function(i){
  lines(f(res[,param]),res[,i],type="b",cex=2,lwd=2,col=i-2)
})
axis(1,at=f(res[,param]),labels=res[,param],las=2)
legend("top",legend=names(res)[-(1:2)],lwd=2,col=1:ncol(res[,-(1:2)]),cex=2,bty="n")
f<-function(i){identity(i)}
plot(0.1,0.1,xlim=range(f(cos[,param]))*c(1,1.2),ylim=range(cos[,-(1:2)]),type="n",xaxt="n",log="x")
lapply(3:ncol(cos),function(i){
  lines(f(cos[,param]),cos[,i],type="l",cex=2,lwd=2,col=i-1)
  text(tail(f(cos[,param]),1),tail(cos[,i],1),label=names(cos)[i],adj=c(-0.1,0.5),cex=0.5,lwd=1,col=i-1)
})
axis(1,at=f(cos[,param]),labels=cos[,param],las=2)
legend("top",legend=substr(names(cos)[-(1:2)],1,10),lwd=2,col=1:ncol(cos[,-1]),cex=1,ncol=4,bty="n")


n<-50
v<-all.vars(model[[3]])
v2<-v[grep("square",v)]
v1<-setdiff(v,v2)
lp<-newdata(x=x[,v1,drop=FALSE],v=v1,n=n,fun=mean,list=FALSE)
if(length(v2)){
  lp<-lapply(lp,function(i){
    a<-as.data.frame(lapply(i[,gsub("square","",v2),drop=FALSE],"^",2))
    names(a)<-v2
    res<-cbind(i,a)
    res[,order(names(res))]
  })
}

par(mfrow=n2mfrow(length(v1),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
for(k in seq_along(v1)){
  
  betas<-m$summary.fixed[,1,drop=FALSE]
  newd<-cbind(intercept=1,as.matrix(lp[[v1[k]]][,v]))
  dimnames(newd)[[2]][1]<-"(Intercept)"
  if(!all(dimnames(newd)[[2]]==row.names(betas))){stop("betas do not correspond to newd")}
  fixed<-newd %*% betas[,1]
  
  y<-exp(fixed)
  plot(lp[[v1[k]]][,v1[k]],y,type="l",ylim=c(0,500))
  mtext(side=3,line=-1.5,adj=0.05,text=v1[k],cex=0.7)
}



library(surveillance)
library(scoringRules)


size<-0.5
#x<-rnbinom(1,mu=50,size=size)
#mu<-rnbinom(1,mu=50,size=size)

#x<-c(seq(0,1,by=0.1),1:50)*2
x<-rep(5,50)
#mu<-rep(0.01,length(x))
mu<-seq(0,20,length.out=50) #c(seq(0,1,by=0.1),1:50)*1
#sigma<-mu*(1+(mu/size))
#trunc<-ceiling(mu+40*sigma)
#by<-1
#sum((pnbinom(seq(0,trunc,by=by),mu=mu,size=size)-as.integer(mu<=seq(0,trunc,by=by)))^2)
s1<-identity(scores(x,mu,size=size,which="rps"));s1  

### Equation in Czado et al. 2009
s2<-identity(sapply(seq_along(mu),function(i){
  X<-rnbinom(100000,mu=mu[i],size=size)
  Xp<-rnbinom(100000,mu=mu[i],size=size)
  mean(abs(X-x[i]))-((1/2)*(mean(abs(X-Xp))))
}));s2

par(mfrow=c(1,2))
plot(mu,s1,type="l",ylim=range(c(s1,s2)),yaxs="i",xaxs="i")
plot(mu,s2,type="l",ylim=range(c(s1,s2)),yaxs="i",xaxs="i")


rpsmc<-function(x,mu,size=NULL,n=10000){
  sapply(seq_along(mu),function(i){
    X<-rnbinom(n,mu=mu[i],size=size)
    Xp<-rnbinom(n,mu=mu[i],size=size)
    mean(abs(X-x[i]))-((1/2)*(mean(abs(X-Xp))))
  })
}
x<-rnbinom(n,mu=10,size=size);hist(x)
mu<-rgamma(n,0.1,0.01);hist(mu)

iter<-10
l<-lapply(1:iter,function(i){
  size<-0.5
  n<-50
  x<-rnbinom(n,mu=5,size=size)
  mu<-rgamma(n,0.1,0.01)
  c(sum(scores(x,mu,size=size,which="rps")),sum(rpsmc(x,mu,size=size,n=10000)))
})
m<-do.call("rbind",l)
matplot(m,pch=c(1,3),xlab="iter",ylab="score")





### scoringRules
crps_nbinom(y=x,size=size,mu=mu)


library(vegan)

#p<-rda(xs@data[xs$db!="map",vars],scale=TRUE)
p<-rda(ds@data[,vars],scale=TRUE)
s<-scores(p,display="sites")

plot(s[,1],s[,2],pch=1,col="white")
points(s[ds$db=="map",1],s[ds$db=="map",2],pch=1,col="red")
points(s[ds$db!="map" & ds$region=="Montréal",1],s[ds$db!="map" & ds$region=="Montréal",2],pch=1,col="blue")


}






