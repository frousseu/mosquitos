

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
library(surveillance)

registerDoParallel(detectCores()-2) 
getDoParWorkers()

load("mosquitos.RData")

inla.setOption(inla.mode="experimental")
year<-c(2003:2016)#c(2003:2006,2013:2016);
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

vals<-list(intercept=1/5^2,jul=1/30^2,julsquare=1/30^2,default=1/30^2) #5-30
control.fixed<-list(prec=vals,mean=list(intercept=-20,default=0),expand.factor.strategy = "inla")

vars<-names(ds)[grep("jul|tmax|prcp|anom|forest|agriculture|water|urban|pond|swamp|pasture|crop|wet|barren",names(ds))]
vars<-vars[-grep("CQ|tmax90|tmax30",vars)]

model<-formula(paste0("sp~",paste(vars,collapse="+")))

#model <- sp ~ jul + julsquare + forest50 + urban50 + urban1000 + agriculture1000  + tmax7 + tmax2 + prcp30

#sigma<-c(0.1,0.2,0.5,0.7,1,2,3,5,10,30,100)
sigma<-sort(c(0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.5,1,5,10,20))
#sigma<-sort(c(0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.5,1,5,10,20))
#sigma<-sigma[seq(1,length(sigma),by=2)]
#cv<-c(2003:2005,2013:2015)
cv<-c("Nord","Sud","Montréal","Laval")
cf<-control.fixed
cvsco<-vector(mode="list",length=length(sigma))

cvsco<-foreach(i=seq_along(sigma),.packages=c("INLA","surveillance")) %dopar% {
  sco<-vector(mode="list",length=length(cv))
  co<-vector(mode="list",length=length(cv))
  cf$prec$default<-1/(sigma[i]^2)
  for(j in seq_along(cv)){
    
    x<-xs@data
    x[vars]<-lapply(x[vars],scale)
    x$julsquare<-x$jul^2
    obs<-x$sp[x$region%in%cv[j]]
    x$sp[x$region%in%cv[j]]<-NA

    m <- inla(model,data=x, 
          control.predictor=list(compute=TRUE,link=1), 
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=cf,
          control.inla = list(strategy='gaussian',int.strategy = "eb"),
          num.threads="1:1",
          verbose=FALSE,
          control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
          family="nbinomial")#"zeroinflatednbinomial1"
    
    pred<-m$summary.fitted.values[x$region%in%cv[j],1]
    size<-m$summary.hyperpar[1,1]

    sco[[j]]<-colSums(scores(obs,pred,size=size,which=c("logs","dss","ses")))
    co[[j]]<-t(m$summary.fixed[,1,drop=FALSE])
    
    #print(m$summary.fixed[,1:5])
    #print(m$summary.hyperpar[,1:5])
    print(range(pred))
    print(paste(sigma[i],"-",cv[j]))
    
  }
  #cvsco[[i]]<-colSums(do.call("rbind",sco))
  list(colSums(do.call("rbind",sco)),colMeans(do.call("rbind",co)))
  
}

#cvsco<-cvsco[c(1:13)]
#sigma<-sigma[c(1:13)]
res<-scale(do.call("rbind",lapply(cvsco,"[[",1)))
res<-cbind(sigma,res)
res<-as.data.frame(res)
cos<-do.call("rbind",lapply(cvsco,"[[",2))
cos<-cbind(sigma,cos)
cos<-as.data.frame(cos)
saveRDS(res,"res.rds")
saveRDS(cos,"cos.rds")



res<-readRDS("C:/Users/God/Documents/mosquitos/res.rds")
cos<-readRDS("C:/Users/God/Documents/mosquitos/cos.rds")
cos<-cos[,!names(cos)%in%c("(Intercept)")]
#cos<-cos[-(1:2),]

#matplot(res[,-1],type="b",lwd=2,cex=2,xaxt="n")
#axis(1,at=seq_along(res$sigma),labels=res$sigma,las=2)
#legend("top",legend=names(res)[-1],lwd=2,col=1:ncol(res[,-1]))

par(mfrow=c(1,2),mar=c(4,3,0.5,0.5))
f<-function(i){identity(i)}
plot(0.1,0.1,xlim=range(f(res$sigma)),ylim=range(res[,-1]),type="n",xaxt="n",log="x")
lapply(2:ncol(res),function(i){
  lines(f(res$sigma),res[,i],type="b",cex=2,lwd=2,col=i-1)
})
axis(1,at=f(res$sigma),labels=res$sigma,las=2)
legend("top",legend=names(res)[-1],lwd=2,col=1:ncol(res[,-1]),cex=2,bty="n")
f<-function(i){identity(i)}
plot(0.1,0.1,xlim=range(f(cos$sigma)),ylim=range(cos[,-1]),type="n",xaxt="n",log="x")
lapply(2:ncol(cos),function(i){
  lines(f(cos$sigma),cos[,i],type="l",cex=2,lwd=2,col=i-1)
  text(tail(f(cos$sigma),1),tail(cos[,i],1),label=names(cos)[i],adj=c(-0.1,0.5),cex=0.5,lwd=1,col=i-1)
})
axis(1,at=f(cos$sigma),labels=cos$sigma,las=2)
legend("top",legend=names(cos)[-1],lwd=2,col=1:ncol(cos[,-1]),cex=1,ncol=4,bty="n")


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
x<-rnbinom(100,mu=50,size=size)
mu<-rnbinom(100,mu=50,size=size)
sigma<-mu*(1+(mu/size))
trunc<-ceiling(mu+40*sigma)
by<-1
#sum((pnbinom(seq(0,trunc,by=by),mu=mu,size=size)-as.integer(mu<=seq(0,trunc,by=by)))^2)
sum(scores(x,mu,size=size,which="rps"))  

### Equation in Czado et al. 2009
sum(sapply(seq_along(mu),function(i){
  X<-rnbinom(100000,mu=mu[i],size=size)
  Xp<-rnbinom(100000,mu=mu[i],size=size)
  mean(abs(X-x[i]))-((1/2)*(mean(abs(X-Xp))))
}))

### scoringRules
crps_nbinom(y=x,size=size,mu=mu)











