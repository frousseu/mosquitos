

# scp C:/Users/God/Documents/mosquitos/models.R root@do:/root

# nohup Rscript models.R --no-save > verbose.out 2>&1 &

# top -o %MEM

library(INLA)

#load("mosquitos.RData")

load("mosquito_models_priors.RData")

#### Subset data #############################################
inla.setOption(inla.mode="experimental")
year<-c(2015)#c(2003:2006,2013:2016);
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
xs<-xs[xs$nbtraps>=10,]
cat(paste(sort(unique(xs$week)),collapse=", "))

#xs$jul<-(xs$jul-mean(xs$jul[xs$db!="map"]))/sd(xs$jul[xs$db!="map"])
#xs$julsquare<-xs$julsquare^2

#### Mesh #####################################################
edge<-6
domain <- inla.nonconvex.hull(coordinates(ds),convex=-0.015, resolution = c(100, 100))
mesh<-inla.mesh.2d(loc.domain=coordinates(ds),max.edge=c(edge,3*edge),offset=c(edge,2*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
#par(mfrow=c(1,2))
plot(mesh,asp=1)
plot(Q,col="grey90",border="white",add=TRUE,lwd=2)
plot(mesh,asp=1,add=TRUE)
plot(xs[!xs$db%in%"map",],add=TRUE,pch=16,col="forestgreen")
#plot(xs,add=TRUE,pch=16,col=alpha("red",0.5))
plot(mappingzone,add=TRUE)
mtext(side=3,line=-2,text=paste("edge =",edge,"km"),font=2,cex=1.5)


#### Restrict predictions ######################################
xsmap<-xs[xs$db=="map",]
xs<-xs[xs$db!="map",]
keep<-expand.grid(year=year[length(year)],week=32)
keep<-sort(sapply(1:nrow(keep),function(i){paste(keep$year[i],keep$week[i],sep="_W")}))
xsmap<-xsmap[xsmap$week%in%keep,]

stdevprior<-c(12)#c(0.00001,0.0001,0.001,0.01,0.03,0.05,0.07,0.1,0.2,1,2,10)
#stdev<-vector(mode="list",length=length(stdevprior))
#names(stdev)<-stdevprior  

for(ii in seq_along(stdevprior)){

print(stdevprior[[ii]])
  
#### SPDE #################################################
spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  constr = FALSE, # not exactly sure what this does
  prior.range=c(5, 0.01), ### P(practic.range<0.05)=0.01
  prior.sigma=c(stdevprior[ii], 0.01)) ### P(sigma>1)=0.01

#### Priors on hyperpar ##################################
#h.spec <- list(theta=list(prior="pc.prec", param=c(0.5,0.5)), rho=list(prior="pc.cor1", param=c(0.9,0.9)))
#h.spec <- list(theta = list(prior="pc.prec", param=c(1, NA)),
#               rho = list(prior="pc.cor0", param=c(0.1, NA)))

h.spec <- list(#theta=list(prior='pc.prec', param=c(0.5, 0.5)))#,
  #rho = list(prior="pc.cor0", param=c(0.7,0.3)))
  rho = list(prior="pc.cor1", param=c(0.9,0.25))) #0.9 0.25
prec.prior <- list(prior='pc.prec', param=c(1, 0.05))
#h.spec <- list(theta = list(prior = "betacorrelation",param=c(1,3),initial=-1.098))
#hist(rbeta(10000,1,3))


#### Model formula ########################################
#model <- y ~ -1 + intercept + jul + julsquare + forest50 + urban50 + urban1000 + agriculture1000  + tmax7 + tmax2 + prcp30 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + jul + julsquare + forest50 + forest1000 + prcp30 + tmax30 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
model <- y ~ -1 + intercept + jul + julsquare + forest1000 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + jul + julsquare + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#model <- y ~ -1 + intercept + jul + julsquare + forest100 + urban100 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#model <- y ~ -1 + intercept + jul + julsquare + forest + urban + tmax15 + tmax1
#model <- y ~ -1 + intercept + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#formulae <- y ~ 0 + w + f(spatial, model=spde) + f(week,model="rw1")
#formulae <- y ~ 0 + w + f(spatial, model=spde, group=spatial.group,control.group=list(model='exchangeable')) 
v<-setdiff(all.vars(model),c("y","spatial","intercept","spde","spatial.group","h.spec")) 

#### Priors on fixed effects ##############################
vals<-list(intercept=1/35^2,default=1/30^2) #5-30
control.fixed<-list(prec=vals,mean=list(intercept=0,default=0),expand.factor.strategy = "inla")

#### Newdata ########################################
n<-50
v2<-v[grep("square",v)]
v1<-setdiff(v,v2)
lp<-newdata(x=xs@data[,v1,drop=FALSE],v=v1,n=n,fun=mean,list=FALSE)
if(length(v2)){
  lp<-lapply(lp,function(i){
    a<-as.data.frame(lapply(i[,gsub("square","",v2),drop=FALSE],"^",2))
    names(a)<-v2
    res<-cbind(i,a)
    res[,order(names(res))]
  })
}
#model<-formula(paste("y~-1+",paste(v,collapse="+")))
#lp<-lapply(lp,function(i){mmatrix(model,i)})
#lpmed<-mmatrix(model,newdata(x=xs[,v,drop=FALSE],v=v,n=1,fun=median,list=FALSE)[[1]][rep(1,length(g)),])[[1]]


#### Make index ##########################################
k<-length(unique(xs$week))
iset<-inla.spde.make.index("spatial",n.spde=spde$n.spde,n.group=k)


#### A matrix ##############################################
gs<-sort(unique(xs$week))
gs<-match(xsmap$week,gs)
Aest<-inla.spde.make.A(mesh=mesh,loc=coordinates(xs),group=as.integer(factor(xs$week))) 
Amap<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap),group=gs) 
#Apre<-inla.spde.make.A(mesh=mesh,loc=matrix(c(600,5050),ncol=2)[rep(1,n),,drop=FALSE],group=rep(12,n))
#Apre<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap),group=as.integer(factor(xsmap$week)))


#### Stacks ################################################
isetmap<-lapply(iset,"[",iset$spatial.group%in%gs)
stackest<-inla.stack(tag='est',data=list(y=xs$sp),A=list(Aest,1),effects=list(c(iset,list(intercept=1)),xs@data)) 
stackmap<-inla.stack(tag='map',data=list(y=xsmap$sp),A=list(Amap,1),effects=list(c(isetmap,list(intercept=1)),xsmap@data)) 
#stackpre<-inla.stack(tag='pre',data=list(y=xs$sp),A=list(Aest,1),effects=list(c(iset,list(intercept=1)),xs@data)) 
stackfull<-inla.stack(stackest,stackmap)


#### Stack vars ############################################
fixgroup<-ceiling(length(sort(unique(xs$week)))/2) # which group to use for var pred
for(i in seq_along(v1)){
  le<-nrow(lp[[v1[i]]])
  if(le!=n){
    #AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,le),,drop=FALSE]) # for categorical variables
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(600,5050),ncol=2)[rep(1,n),,drop=FALSE],group=rep(fixgroup,n))
  }else{
    #AA<-Apn # for numerical variables
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(600,5050),ncol=2)[rep(1,n),,drop=FALSE],group=rep(fixgroup,n))
  }
  stack<-inla.stack(tag=v1[i],data=list(y=NA),A=list(AA,1),effects=list(c(lapply(iset,"[",iset$spatial.group==fixgroup),list(intercept=1)),lp[[v1[i]]]))     
  stackfull<-inla.stack(stackfull,stack)
}

#### Full stack ############################################
stackfull<-inla.stack(stackest)


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
          num.threads="2:1",
          verbose=FALSE,
          control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
          #control.mode = list(result = m, restart = TRUE)), # to rerun the model with NA predictions according to https://06373067248184934733.googlegroups.com/attach/2662ebf61b581/sub.R?part=0.1&view=1&vt=ANaJVrHTFUnDqSbj6WTkDo-b_TftcP-dVVwK9SxPo9jmPvDiK58BmG7DpDdb0Ek6xypsqmCSTLDV1rczoY6Acg_Zb0VRPn1w2vRj3vzHYaHT8JMCEihVLbY
          family="zeroinflatednbinomial1")#"zeroinflatednbinomial1"


stdev[[as.character(stdevprior[ii])]]<-m$summary.hyperpar
save.image("mosquito_models_priors.RData")

}

cat(paste("edge =",edge),"\n\n")
cat(paste(sort(unique(xs$week)),collapse=", "),"\n\n")
cat("inla.mode =",inla.getOption("inla.mode"),"\n\n")

#save.image("mosquito_models_priors.RData")


