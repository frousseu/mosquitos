
set.seed(123)

ns<-c(1,5,10,100,500,1000,1200,1300,1350,1400,1450,1500,1550)
li<-vector(mode="list",length=length(ns)+1)
names(li)<-c(0,ns)
li[[as.character(0)]]<-m0

for(j in ns){

#### Subset data #############################################
year<-c(2015);
spcode<-"VEX_"
weeks<-29:33
xs<-ds[ds$year%in%year,]
sp<-names(xs)[grep(spcode,names(xs))]
xs$sp<-xs@data[,sp]
xs<-xs[order(xs$week),]
##xs<-xs[xs$week%in%paste0(year,"_W",weeks),]
xs<-xs[substr(xs$week,7,8)%in%weeks,]

#### Mesh #####################################################
edge<-10
domain <- inla.nonconvex.hull(coordinates(ds),convex=-0.015, resolution = c(100, 100))
mesh<-inla.mesh.2d(loc.domain=coordinates(ds),max.edge=c(edge,3*edge),offset=c(edge,2*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
#par(mfrow=c(1,2))
plot(mesh,asp=1)
plot(Q,col="grey90",border="white",add=TRUE,lwd=2)
plot(mesh,asp=1,add=TRUE)
plot(xs[xs$db%in%"map",],add=TRUE,pch=16,col="firebrick")
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

#plot(xsmap)
rem<-concaveman(st_as_sf(xsmap),1)
#plot(st_geometry(rem),add=T,border="red")
rem<-as(as(st_buffer(rem,-20),"Spatial"),"SpatialPolygons")
#plot(rem,add=T,border="red")
ooo<-over(xsmap,rem)
xsmap<-xsmap[!is.na(ooo),]
xsmap<-xsmap[sample(1:nrow(xsmap),j),]
#xsmap<-xsmap[coordinates(xsmap)[,1]>=660,]

plot(mesh,asp=1)
plot(Q,col="grey90",border="white",add=TRUE,lwd=2)
plot(mesh,asp=1,add=TRUE)
plot(xsmap,add=TRUE,pch=16,col="firebrick")
plot(xs,add=TRUE,pch=16,col="forestgreen")
#plot(xs,add=TRUE,pch=16,col=alpha("red",0.5))
plot(mappingzone,add=TRUE)
mtext(side=3,line=-2,text=paste("edge =",edge,"km"),font=2,cex=1.5)


#### SPDE #################################################
spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  constr = FALSE, # not exactly sure what this does
  prior.range=c(5, 0.1), ### P(practic.range<0.05)=0.01
  prior.sigma=c(1, 0.05)) ### P(sigma>1)=0.01

#### Priors on hyperpar ##################################
#h.spec <- list(theta=list(prior="pc.prec", param=c(0.5,0.5)), rho=list(prior="pc.cor1", param=c(0.9,0.9)))
#h.spec <- list(theta = list(prior="pc.prec", param=c(1, NA)),
#               rho = list(prior="pc.cor0", param=c(0.1, NA)))

h.spec <- list(#theta=list(prior='pc.prec', param=c(0.5, 0.5)))#,
  #rho = list(prior="pc.cor0", param=c(0.7,0.3)))
  rho = list(prior="pc.cor1", param=c(0.9,0.25))) #0.9 0.25
prec.prior <- list(prior='pc.prec', param=c(1, 0.5))
#h.spec <- list(theta = list(prior = "betacorrelation",param=c(1,3),initial=-1.098))
#hist(rbeta(10000,1,3))


#### Model formula ########################################
model <- y ~ -1 + intercept + jul + julsquare + forest50 + urban50 + urban1000 + agriculture1000  + tmax7 + tmax2 + prcp30 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + jul + julsquare + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec))
#model <- y ~ -1 + intercept + jul + julsquare + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#model <- y ~ -1 + intercept + jul + julsquare + forest100 + urban100 + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#model <- y ~ -1 + intercept + jul + julsquare + forest + urban + tmax15 + tmax1
#model <- y ~ -1 + intercept + f(spatial, model=spde, group=spatial.group,control.group=list(model='ar1', hyper=h.spec)) 
#formulae <- y ~ 0 + w + f(spatial, model=spde) + f(week,model="rw1")
#formulae <- y ~ 0 + w + f(spatial, model=spde, group=spatial.group,control.group=list(model='exchangeable')) 
v<-setdiff(all.vars(model),c("y","spatial","intercept","spde","spatial.group","h.spec")) 

#### Priors on fixed effects ##############################
vals<-list(intercept=1/5^2,default=1/30^2) #5-30
control.fixed<-list(prec=vals,mean=list(intercept=-20,default=0),expand.factor.strategy = "inla")

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
#stackfull<-inla.stack(stackest)
stackfull<-inla.stack(stackest,stackmap)


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
         family="zeroinflatednbinomial1")#"zeroinflatednbinomial1"
print(j)

li[[as.character(j)]]<-m
}



pars<-do.call("rbind",lapply(seq_along(li),function(i){
  print(i)
  cbind(param=row.names(li[[i]]$summary.hyper),n=as.integer(names(li))[i],data.frame(li[[i]]$summary.hyper))
}))
row.names(pars)<-NULL
pars<-pars[,c("param","n","mean","X0.025quant","X0.975quant")]
#pars<-NaRV.omit(pars)

pars<-split(pars,~param)

f<-function(x){x^(1/2)}
par(mfrow=n2mfrow(length(pars)),mar=c(2,2,1,1),oma=c(4,4,4,1))
lapply(pars,function(k){
  pe<-c("mean","X0.025quant","X0.975quant")
  plot(0,0,type="n",xlim=f(range(ns)),ylim=c(-1,1)*0.5*diff(range(unlist(NaRV.omit(k[k$n<=100,])[-(1:2)])))+range(unlist(NaRV.omit(k[k$n<=100,])[-(1:2)])),yaxt="n",xaxt="n",log="")
  axis(1,at=f(as.integer(names(li))),label=as.integer(names(li)),las=2,cex.axis=0.75,mgp=c(3,0.45,0),tcl=-0.3)
  axis(2,las=2,mgp=c(3,0.45,0),tcl=-0.3)
  lapply(pe,function(j){
    i<-NaRV.omit(k[,c("param","n",j)])
    lines(f(i$n),i[,j],lty=ifelse(j=="mean",1,2))
    points(f(i$n),i[,j],cex=1.5)
  })
  abline(v=f(as.integer(names(li))),lty=3,col=gray(0,0.2))
  mtext(side=3,adj=0.03,line=-2,text=k$param[1])
})
mtext(side=3,line=1,text=paste0("inla.mode = \"",inla.getOption("inla.mode"),"\""),cex=2,font=2,outer=TRUE)
mtext(side=1,line=1,text="Number of NA values to predict",outer=TRUE,cex=1.5)
mtext(side=2,line=1,text="Hyper parameters",outer=TRUE,cex=1.5)

