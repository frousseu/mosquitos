
nb<-c(2,3,4,6,8,10,15,20,25,30,37,50)
times<-NULL
for(j in nb){

colSums(xs@data[,5:33])
#xs1<-ds[ds$year%in%c("2015"),]
#xs2<-ds[ds$year%in%c("2016"),]
#xs1$sp<-log(xs1$A9+1)
#xs2$sp<-log(xs2$A9+1)
xs<-ds[ds$year%in%c("2014","2015","2016"),]
#xs$sp<-log(xs$A29+1)
xs$sp<-xs$A29
#xs<-rbind(xs1,xs2)
#xs<-xs1
#xs$Week<-sample(xs$Week)
#xs$sp<-xs$A29
xs<-xs[order(xs$Annee,xs$Week),]
xs$Week<-paste(xs$Annee,xs$Week,sep="_")
xs$week<-as.integer(substr(xs$Week,6,7))
xs$week<-xs$week
xs$week2<-xs$week^2
#xs<-xs[xs$week%in%c(20:(20+j-1)),]
xs<-xs[xs$Week%in%levels(factor(xs$Week))[1:j],]

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
x1<-Sys.time()
m <- inla(formulae,  data=inla.stack.data(sdat), 
          control.predictor=list(compute=TRUE, A=inla.stack.A(sdat),link=1), 
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=list(expand.factor.strategy='inla'),
          control.inla = list(int.strategy = "eb"),
          num.threads=3,
          verbose=FALSE,
          family="nbinomial")#"zeroinflatednbinomial1"
x2<-Sys.time()
times<-c(times,as.numeric(difftime(x2,x1,units="sec")))

}

plot(nb,times,type="l")




