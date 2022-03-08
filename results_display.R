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