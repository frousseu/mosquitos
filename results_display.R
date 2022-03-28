
library(RcppArmadillo)
library(Rcpp)
library(matrixStats)
library(knitr)
library(kableExtra)
library(webshot)
library(png)
library(grid)
library(magick)
library(pracma)

options(device = "X11")
grDevices::windows.options(record=TRUE)
Sys.setlocale("LC_ALL","English")
# first set working directory
# all files should be in this folder

## RESULTS ##################################################

load("VEX_model_outputs.RData")

ls()[sapply(ls(),function(i){
  obj<-paste0("\\b",i,"\\b")
  as.logical(length(grep(obj,readLines("C:/Users/God/Documents/mosquitos/results_display.R"),value=T)))
})]

# for master matrix * vector multiplications when combining posterior samples
# code stolen from https://stackoverflow.com/questions/51054227/why-is-this-naive-matrix-multiplication-faster-than-base-rs
arma_code <- 
  "arma::vec arma_mm(const arma::mat& m, const arma::vec& v) {
       return m * v;
   };"
arma_mm = cppFunction(code = arma_code, depends = "RcppArmadillo")


niceround<-function(x){
  res<-lapply(x,function(e){
    #e<-exp(a)
    if(e<=0.00001){return(round(e,6))}
    if(e<=0.0001){return(round(e,5))}
    if(e<=0.001){return(round(e,4))}
    if(e<=0.01){return(round(e,3))}
    if(e<=0.1){return(round(e,2))}
    if(e<2){return(formatC(e,format="f",digits=1))}
    if(e>=2){return(round(e,0))}
  })
  sapply(res,as.character)
}

getlabels<-function(x){
  v<-c(jul="Date", agriculture50="% of agriculture within 50 m", forest50="% of forests within 50 m", anom2="Mean temperature anomalies in previous 2 days (\u00B0C)", prcp2="Mean daily precipitations in previous 2 days (mm)", agriculture1000="% of agriculture within 1 km", 
    forest1000="% of forests within 1 km", anom7="Mean temperature anomalies in previous 7 days (\u00B0C)", prcp7="Mean daily precipitations in previous 7 days (mm)", anom30="Mean temperature anomalies in previous 30 days (\u00B0C)", prcp30="Mean daily precipitations in previous 30 days (mm)", anom90="Mean temperature anomalies in previous 90 days (\u00B0C)", 
    prcp90="Mean daily precipitations in previous 90 days (mm)", urban50="% of urban within 50 m", urban1000="% of urban within 1 km")
  #print(v)
  paste(v[match(x,names(v))],x,sep=" - ")
}

# this determines the order of variables presented
ov<-c("jul","anom2","anom7","anom30","anom90","prcp2","prcp7","prcp30","prcp90","agriculture50","agriculture1000","forest50","forest1000","urban50","urban1000","spatial")

#### DIC ####################################################

msel<-c(dics,mfixed$dic$dic)
delta_dics<-round(msel-min(msel),2)
mods<-sapply(c(spmodels,fixed),function(i){
  av<-all.vars(i)
  av<-av[order(match(av,ov))]
  paste(av[!av%in%c("y","knots","lognights","spde")],collapse=" + ")
})
resdics<-data.frame(mods,delta_dics)
resdics<-resdics[order(resdics$delta_dics),]
resdics$Model<-row.names(resdics)
resdics$Model[nrow(resdics)]<-resdics$Model[1]
resdics$Model<-paste0(substr(resdics$Model,1,3),formatC(as.integer(sapply(strsplit(resdics$Model,"_"),"[",2)),width=2,flag="0",format="d"))
resdics$Model[nrow(resdics)]<-paste0(resdics$Model[nrow(resdics)],"nonspatial")

names(resdics)<-c("Variables","\u0394DIC","Model")
resdics$Variables<-gsub("jul","ns(jul)",resdics$Variables)
resdics<-resdics[,c(3,1,2)]

resdics %>%
  kable(row.names = FALSE, align=c("l","l","r")) %>%
  kable_classic(full_width = FALSE, html_font = "Helvetica") %>% 
  row_spec(0, bold = T, background = "#EEEEEE",align="c") %>% 
  #kable_styling(full_width = FALSE, font_size = 12) %>% 
  save_kable(file = file.path("C:/Users/God/Downloads",paste0(spcode,"dics.png")),density=500,zoom=2)
file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"dics.png")))

#resdics %>%
#  kbl(row.names = FALSE) %>%
#  kable_classic(full_width = FALSE, html_font = "Helvetica") %>%
#  kable_styling(latex_options = c("striped", "scale_down")) %>%
#  row_spec(0, bold = T, background = "#EEEEEE",align="c") %>% 
#  #kable_styling(full_width = FALSE, font_size = 12) %>% 
#  #as_image()
#  save_kable(file = file.path("C:/Users/God/Downloads",paste0(spcode,"dics.png")),density=500)
#file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"dics.png")))


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


#### Combine posteriors ########################################

posteriors<-c(m$marginals.fixed,m$marginals.hyper)
posteriors<-posteriors[!names(posteriors)%in%1:50 & !grepl("size|Stdev|Range",names(posteriors))]
posteriors<-posteriors[rev(order(abs(sapply(posteriors,function(i){i[which.max(i[,2]),1]}))))]
posteriors<-c(m$marginals.hyper,posteriors)
names(posteriors)[1:3]<-c("Size for nbinom","Range","SD")
posteriors<-posteriors[4:length(posteriors)]

## remove tail values that are not seeabje and that are extending xlim values for nothing
par(mfrow=n2mfrow(length(posteriors),asp=1.49),mar=c(3,2.5,1,1))
for (j in 1:length(posteriors)) {
  k<-posteriors[[j]][,2]>=1e-2*max(posteriors[[j]][,2])
  posteriors[[j]]<-posteriors[[j]][k,]
}

topoly<-function(x,y,...){
  polygon(c(x,rev(x)),c(y,rep(0,length(y))),...)  
}

posneg<-function(x){
  neg<-trapz(x[x[,1]<0,1],x[x[,1]<0,2])  
  pos<-trapz(x[x[,1]>0,1],x[x[,1]>0,2]) 
  if(neg>0.975){return(adjustcolor("steelblue",0.5))}
  if(pos>0.975){return(adjustcolor("firebrick4",0.5))}
  adjustcolor("ivory3",0.8)
}

png(file.path("C:/Users/God/Downloads",paste0(spcode,"posteriors.png")),width=4,height=6,units="in",res=300,pointsize=11)
xlim<-range(do.call("rbind",posteriors)[,1])
#xlimh<-NULL
xlim<-c(-1.0,0.8)
ylim<-c(0,range(do.call("rbind",posteriors)[,2])[2])
par(mfrow=c(length(posteriors),1),mar=c(0,0,0,0),oma=c(3,1,2,1))
for (j in 1:length(posteriors)) {
  plot(posteriors[[j]][,1],posteriors[[j]][,2],type='n',xlab="",ylab='Density',xlim=xlim,ylim=ylim,yaxt="n",xaxt="n",bty="n",xpd=TRUE)#
  abline(h=seq(0,max(ylim),length.out=5),lty=3,col="grey80")
  #grid(col="grey90")
  #box(col="grey99",lwd=5)
  lines(c(0,0),ylim,lty=3,lwd=1,col="grey45")
  at<-pretty(xlim,8)
  lapply(at,function(x){lines(c(x,x),ylim,lty=3,col="grey80")})
  lines(c(0,0),ylim,lty=1,lwd=1,col="grey70")
  if(j==length(posteriors)){
    axis(1,at=at,labels=at,tcl=-0.2,mgp=c(1.5,0.25,0),cex.axis=1,font=2,col="grey70")
  }
  text(par("usr")[1]+abs(diff(par("usr")[c(1,3)]))*0.0,par("usr")[2]+diff(par("usr")[c(2,4)])*0.90,names(posteriors)[j],cex=1.5,font=2,adj=c(0,1),col=gray(0,0.9))
  topoly(posteriors[[j]][,1],posteriors[[j]][,2],border=NA,col=posneg(posteriors[[j]]))#
  mtext(outer=TRUE,side=1,line=2,text="Coefficient",font=2)
  mtext(outer=TRUE,side=3,line=0.25,text=substr(spcode,1,3),font=2,cex=1.5,adj=0.01)
}
dev.off()
file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"posteriors.png")))


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
r<-field[["mean"]]
r<-mask(r,mappingzone)


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

v1m<-c("jul",v1[v1%in%row.names(m$summary.fixed)])
v1m<-v1m[order(match(v1m,ov))]

png(file.path("C:/Users/God/Downloads",paste0(spcode,"marginal_effects.png")),width=10,height=8,units="in",res=300,pointsize=11)

par(mfrow=n2mfrow(length(v1m),asp=3.45/2),mar=c(3,2.5,1,1),oma=c(4,4,0,0))
for(k in seq_along(v1m)){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
    fixed<-as.matrix(lp[[v1m[k]]][,names(betas)]) %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
    #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    spatial<-as.matrix(AA) %*% wk # stack was Apn in fire
    #}
    p1<-fixed
    p2<-fixed+spatial
    #p<-fixed # ignores spatial part
    list(p1,p2)
  })
  p1<-do.call("cbind",lapply(p,"[[",1))
  p2<-do.call("cbind",lapply(p,"[[",2))
  p1<-t(apply(p1,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  p2<-t(apply(p2,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  p1<-exp(p1)
  p2<-exp(p2)
  if(nrow(lp[[v1m[k]]])==n){
    vals<-bscale(lp[[v1m[k]]][,v1m[k]],v=v1m[k])
    if(TRUE){ # log y scale or not
      if(v1m[k]=="jul"){xlim<-c(135-13,288+8)}else{xlim<-NULL}
      plot(bscale(xs@data[,v1m[k]],v=v1m[k]),xs$sp+1,xlab=getlabels(v1m[k]),font=2,ylab="",yaxt="n",xaxt="n",pch=16,col=gray(0.5,0.2),log="y",xlim=xlim,mgp=c(1.5,0.45,0))
      polygon(c(vals,rev(vals),vals[1]),c(p2[,1],rev(p2[,3]),p2[,1][1])+1,col=alpha("green4",0.25),border=NA)
      polygon(c(vals,rev(vals),vals[1]),c(p1[,1],rev(p1[,3]),p1[,1][1])+1,col=alpha("green4",0.35),border=NA)
      #lines(vals,p2[,2]+1,lwd=1,col="red")
      lines(vals,p1[,2]+1,lwd=1,col="darkgreen")
      at<-c(1,6,11,51,101,501,1001,5001,10001,max(xs$sp)+1)
      axis(2,at=at,label=at-1,mgp=c(2,0.45,0),tcl=-0.2,las=2,font=2,cex.axis=1,gap.axis=0)
      if(v1m[k]=="jul"){
        labd<-paste0("2021-",c("05-01","06-01","07-01","08-01","09-01","10-01","11-01"))
        at<-as.integer(format(as.Date(labd),"%j"))
        lab<-format(as.Date(labd),"%b-%d")
        axis(1,at=at,labels=lab,mgp=c(1.5,0.45,0),font=2,tcl=-0.2)
      }else{
        axis(1,mgp=c(1.5,0.45,0),font=2,tcl=-0.2)
      }
    }else{
      plot(vals,p1[,2],type="l",ylim=c(0,min(c(max(p1[,3]),max(xs$sp))))*1.3,xlab=v1m[k],font=2,ylab="",lty=1,yaxt="n",mgp=c(1.5,0.45,0),tcl=-0.3)
      points(bscale(xs@data[,v1m[k]],v=v1m[k]),xs$sp,pch=1,col=gray(0,0.1))
      lines(vals,p1[,2],lwd=1.5)
      polygon(c(vals,rev(vals),vals[1]),c(p1[,1],rev(p1[,3]),p1[,1][1]),col=alpha("black",0.1),border=NA)
      axis(2,las=2,font=2)
    }
  }else{
    #plot(unique(sort(size[,v[k]])),p[,2],type="l",ylim=c(0,100),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    #points(jitter(as.integer(size[,v[k]])),size$tTotal,pch=16,col=gray(0,0.1))
    #segments(x0=as.integer(unique(sort(size[,v[k]]))),x1=as.integer(unique(sort(size[,v[k]]))),y0=p[,1],y1=p[,3],lty=3)
  }
  
}
mtext(paste("No. of mosquitos per trap"),outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=1)
mtext(paste("Explanatory variables"),outer=TRUE,cex=1.2,side=1,xpd=TRUE,line=1)

dev.off()
file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"marginal_effects.png")))

#### Marginal effects with spatial uncertainty ##########################

# this section is not that useful because it is a prediction for a given location, hence it includes uncertainty in the spatial field

png(file.path("C:/Users/God/Downloads",paste0(spcode,"marginal_effects_spatial.png")),width=12,height=8,units="in",res=300,pointsize=11)

par(mfrow=n2mfrow(length(v1m),asp=1.5),mar=c(3,2.5,1,1),oma=c(0,10,0,0))
for(k in seq_along(v1m)){
  p<-m$summary.fitted.values[index[[v1m[k]]],c("0.025quant","0.5quant","0.975quant")]
  #p[]<-lapply(p,transI)
  dat<-data.frame(lp[[v1m[k]]])
  if(nrow(p)==n){
    vals<-bscale(dat[[v1m[k]]],v=v1m[k])
    if(TRUE){ # log y scale or not
      plot(bscale(xs@data[,v1m[k]],v=v1m[k]),xs$sp+1,xlab=v1m[k],font=2,ylab="",yaxt="n",pch=16,col=gray(0,0.05),mgp=c(1.5,0.45,0),log="y")
      lines(vals,p[,2]+1,lwd=1.5)
      polygon(c(vals,rev(vals),vals[1]),c(p[,1],rev(p[,3]),p[,1][1])+1,col=alpha("black",0.1),border=NA)
      at<-c(1,6,11,51,101,501,1001,5001,10001,max(xs$sp)+1)
      axis(2,at=at,label=at-1,mgp=c(2,0.45,0),tcl=-0.3,las=2,font=2)
    }else{
      plot(vals,p[,2],type="l",ylim=c(0,min(c(max(p[,3]),max(xs$sp)))),xlab=v1m[k],font=2,ylab="",lty=1,yaxt="n",mgp=c(1.5,0.45,0),tcl=-0.3,lwd=1.5)
      #plot(dat[[v1m[k]]],p[,2],type="l",ylim=c(0,300),xlab=v1m[k],font=2,ylab="",lty=1,yaxt="n",mgp=c(2,0.45,0),tcl=-0.3)
      #lines(vals,p[,1],lty=3,lwd=1)
      #lines(vals,p[,3],lty=3,lwd=1)
      points(bscale(xs@data[,v1m[k]],v=v1m[k]),xs$sp,pch=16,col=gray(0,0.07))
      polygon(c(vals,rev(vals),vals[1]),c(p[,1],rev(p[,3]),p[,1][1]),col=alpha("black",0.1),border=NA)
      axis(2,las=2,font=2)
    }
  }else{
    #plot(unique(sort(size[,v1m[k]])),p[,2],type="l",ylim=c(0,100),xlab=v1m[k],font=2,ylab="",lty=1,yaxt="n")
    #segments(x0=as.integer(unique(sort(size[,v[k]]))),x1=as.integer(unique(sort(size[,v[k]]))),y0=p[,1],y1=p[,3],lty=3,lwd=2)
    #points(jitter(as.integer(size[,v[k]]),fac=2),transI(size$tTotal),pch=16,col=gray(0,0.07))
  }
  
}
mtext("No. of mosquitos per trap",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)

dev.off()
file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"marginal_effects_spatial.png")))

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
#xsbuff<-st_coordinates(st_cast(st_buffer(st_as_sf(xs),7),"MULTIPOINT"))[,1:2]
#buf<-concaveman(xsbuff,10)
#buf<-spPolygons(buf,crs=CRS(proj4string(xs)))
#buf<-gBuffer(buf,width=1)
#pred<-mask(pred,buf)
pred<-mask(pred,mappingzone)

cols<-alpha(colo.scale(200,c("darkblue","steelblue3","lightgoldenrod","orange","red3","darkred","grey10")),0.80)
#colssd<-colo.scale((seq(0.01,5,length.out=200))^3,rev(cividis(200)))
colsint<-c("darkblue","steelblue","ivory3","firebrick3","firebrick4")
#colsint<-c("darkgreen","palegreen","ivory3","purple4","grey10")
colsfield<-colo.scale(seq(range(values(pred[["mean.spatial.field"]]),na.rm=TRUE)[1],range(values(pred[["mean.spatial.field"]]),na.rm=TRUE)[2],length.out=200),colsint,center=TRUE)
colssd<-colo.scale((seq(0.01,5,length.out=200))^3,colsfield)

titles<-c("Mean","CI 2.5%","CI 97.5%","SD","Mean Spatial Field","SD Spatial Field")
legtitles<-c("No. of Mosquitos / trap","No. of Mosquitos / trap","No. of Mosquitos / trap","No. of Mosquitos / trap (on the link scale)","No. of Mosquitos / trap (on the link scale)","No. of Mosquitos / trap (on the link scale)")
names(titles)<-names(pred)
names(legtitles)<-names(pred)  



png(file.path("C:/Users/God/Downloads",paste0(spcode,"maps_all.png")),width=16,height=8,units="in",res=300,pointsize=11)

par(mfrow=n2mfrow(nlayers(pred),asp=1.5),mar=c(1,0.5,1,7),oma=c(0,0,0,0),bty="n")
lapply(names(pred),function(i){
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
    zlim<-NULL # limits specific to graph
    if(i %in% quantities[1:3] && TRUE){ # limits determined by CI if TRUE
      zlim<-range(values(pred2[[quantities[1:3]]]),na.rm=TRUE) 
    }
    axis.args=list(at=frange(values(pred2[[i]])),labels=niceround(f(frange(values(pred2[[i]])))),cex.axis=1.2,lwd=0,tck=-0.2,mgp=c(3,0.3,0),lwd.ticks=1)
    legend.args=list(text=legtitles[i], side=4, font=2, line=-2.5, cex=1)
  }else{
    if(i%in%names(meansd)){
      zlim<-NULL
      axis.args=list(at=sort(c(0,frange(values(pred2[[i]])))),labels=sort(round(c(0,f(frange(values(pred2[[i]])))),1)),cex.axis=1.2,lwd=0,tck=-0.2,mgp=c(3,0.3,0),lwd.ticks=1)
      legend.args=list(text=legtitles[i], side=4, font=2, line=-2.5, cex=1)
    }else{
      #zlim<-range(values(pred2[[quantities[1:3]]]),na.rm=TRUE)
      #axis.args=list(at=c(frange(zlim),range(values(pred2[[i]]),na.rm=TRUE)),labels=round(f(c(frange(zlim),range(values(pred2[[i]]),na.rm=TRUE))),0),cex.axis=0.8,lwd=0,tck=-0.2,mgp=c(3,0.3,0),lwd.ticks=1)
      #legend.args=list(text='Nb of mosquitos / trap', side=4, font=2, line=-2.5, cex=0.8)
    }
  }
  plot(pred2[[i]],col=col,zlim=zlim,legend.width=2.5, legend.shrink=1,axis.args=axis.args,legend.args=legend.args,axes=FALSE,box=FALSE,tcl=0.2,mgp=c(1.5,0.0,0),cex.axis=0.7)
  plot(st_geometry(water),border=NA,col="white",add=TRUE)
  show<-unique(xsmap$week) # show only obs from a given year in the week
  show<-paste0(2003:2016,substr(unique(xsmap$week),5,8)) # show obs from all years in a given week
  obs<-xs$sp[xs$week%in%show]
  cexminmax<-c(1,10)  
  ocex<-scales::rescale(c(ifelse(obs==0,1,identity(obs))),to=cexminmax)
  opch<-ifelse(obs==0,4,1)
  plot(xs[xs$week%in%show,],add=TRUE,cex=ocex,pch=opch,lwd=0.2,col=gray(0.2,1))
  #plot(Q,add=TRUE,border=gray(0,0.25))
  #plot(mappingzone,add=TRUE)
  #plot(mesh,add=TRUE)
  mtext(side=3,line=-3,text=titles[i],adj=0.05,font=2,cex=2)
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
  text(posx,posy,label=vleg,adj=c(0.5,4),cex=1.2,xpd=TRUE)
})

dev.off()
file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"maps_all.png")))


#### Map abundance across season #############################

# not done yet and not general enough

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  #grep(paste0(i,":"),row.names(samples[[1]]$latent))  
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
yearpred<-"2003"
#nweights<-grep("spatial",row.names(samples[[1]]$latent))
Amapp<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap)) 
Amapmatrix<-as.matrix(Amapp)

days<-seq(min(xs$jul[xs$year==yearpred]),max(xs$jul[xs$year==yearpred]),length.out=20)
lpr<-foreach(j=seq_along(days),.packages=c("raster")) %do% {
  juls<-lp[["jul"]][which.min(abs(lp[["jul"]]$jul-days[j])),,drop=FALSE]
  standardv<-names(nparams)[!names(nparams)%in%c("intercept","jul","julsquare",1:50)]
  dat<-as.matrix(xsmap@data[,standardv])
  dat<-cbind(dat,data.frame(jul=days[j],julsquare=days[j]^2)[rep(1,nrow(dat)),])
  dat<-cbind(dat,juls[,names(juls)%in%paste0("X",1:50)][rep(1,nrow(dat)),])
  dat<-cbind(intercept=1,dat)

  #betas<-samples[[1]]$latent[nparams]
  betas<-m$summary.fixed[,1]
  names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
  fixed<-as.matrix(dat[,names(betas)]) %*% betas # make sure betas and vars are in the same order
  #wk<-samples[[1]]$latent[nweights]
  wk<-m$summary.random$spatial[,"mean"]
  #spatial<-Amapmatrix %*% wk
  spatial<-arma_mm(Amapmatrix,wk) # ~ 2 times faster than %*%
  #}
  p<-fixed+spatial
  #p<-fixed # ignores spatial part
  #p<-spatial
  #cat("\r",paste(i,"/",max(nsims)," - ",j,"/",length(days)))
  #p<-do.call("cbind",p)
  #p<-cbind(rowQuantiles(p,probs=c(0.0275,0.975),na.rm=TRUE),rowMeans(p,na.rm=TRUE))[,c(1,3,2)]
  #p<-exp(p)
  xsmap$preds<-p[,1]
  pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
  pr<-mask(pr,mappingzone)
  pr<-exp(pr)
  #pr<-disaggregate(pr,fact=1,method="bilinear")
  #print(j)
  cat("\r",j,"/",length(days))
  pr
}
lpr<-lapply(lpr,rast)


jul<-round(days*vscale[["jul"]]["sd"]+vscale[["jul"]]["mean"],0)
datelim<-range(as.Date(format(as.Date(jul,origin=paste0(yearpred,"-01-01")),"%Y-%m-%d")))+c(-5,5)

zlim1<-range(sapply(lpr,function(i){range(values(i),na.rm=TRUE)}))
zlim2<-range(c(sapply(lpr,function(i){range(values(i),na.rm=TRUE)}),xs$sp[xs$year==yearpred]))
zlim<-c(zlim1[1],zlim2[2])

xs2<-xs[xs$year==yearpred,]
climate<-aggregate(.~date,data=xs2@data[,c("date",names(xs2)[grep("anom|tmean|prcp",names(xs2))])],mean)
climate$date<-as.Date(climate$date)
climate[names(climate)[-1]]<-lapply(names(climate)[-1],function(i){bscale(climate[,i],i)})

img <- image_graph(1500, 1200, res = 150)
lapply(seq_along(lpr),function(i){
  j<-round(days[i]*vscale[["jul"]]["sd"]+vscale[["jul"]]["mean"],0)
  xdate<-as.Date(format(as.Date(j,origin=paste0(yearpred,"-01-01")),"%Y-%m-%d"))
  gw<-layout(matrix(c(rep(1,18),2,3,4,5,5,5),ncol=1))
  par(mar=c(0,0,0,0),oma=c(0,0,0,12))
  at<-seq(min(values(log(lpr[[i]])),na.rm=TRUE),max(values(log(lpr[[i]])),na.rm=TRUE),length.out=5)
  #lab<-ifelse(round(exp(at),0)==0,round(exp(at),2),round(exp(at),0))
  lab<-niceround(exp(at))
  labels<-paste(lab,c("min pred.",rep("",length(at)-2),"max pred."))
  plot(log(lpr[[i]]),range=log(zlim),col=cols,asp=1,axes=FALSE,bty="n",plg=list(at=at,labels=labels,cex=1.5,xpd=TRUE),mar=c(1,0,0.5,0))
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
  mtext(side=3,line=-1.6,text=paste(gsub("_","",spcode),yearpred,"  observations:",paste(format(as.Date(range(rd)),"%b-%d"),collapse=" to "),sep="  "),adj=0.025,cex=1.5)
  mtext(side=4,line=-1,text="Number of mosquitos per trap (observed and predicted)",adj=0.5)
  xp<-xmin(lpr[[1]])+((xmax(lpr[[1]])-xmin(lpr[[1]]))*c(0.05,0.13))
  yp<-rep(ymin(lpr[[1]])+((ymax(lpr[[1]])-ymin(lpr[[1]]))*0.90),2)  
  points(xp,yp,pch=21,cex=2,bg=colobs[c(which.min(xxs$sp),which.max(xxs$sp))],col="grey10",lwd=0.4)
  text(xp,yp,label=xxs$sp[c(which.min(xxs$sp),which.max(xxs$sp))],cex=0.7,col="grey10",adj=c(0.5,-1))
  text(xp,yp,label=c("min obs.","max obs."),cex=1,col="grey10",adj=c(1.2,0.5))
  rec<-function(ybottom=-1000,xpd=TRUE){
    rect(xleft=range(as.Date(rd))[1],ybottom=ybottom,xright=range(as.Date(rd))[2],ytop=3000000,border=NA,col=adjustcolor("darkgreen",0.2),xpd=xpd)
  }
  recb<-function(){
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col=adjustcolor("darkgreen",0.1),border=NA)
  }
  par(mar=c(0.5,3,0,0))
  plot(climate$date,climate$tmean2,xlim=datelim,type="n",xaxt="n",yaxt="n",bty="n")
  recb()
  lines(climate$date,climate$tmean2,lwd=2,col="tomato")
  abline(h=mean(climate$tmean2),lty=3,lwd=1,col="tomato")
  axis(2,mgp=c(2,0.45,0),tcl=-0.2,las=2,font=2,cex.axis=0.75,gap.axis=0)
  rec()
  
  par(mar=c(0.5,3,0,0))
  plot(climate$date,climate$anom2,xlim=datelim,type="n",xaxt="n",yaxt="n",bty="n")
  recb()
  lines(climate$date,climate$anom2,lwd=2,col="tomato4")
  abline(h=0,lty=3,lwd=1,col="tomato4")
  axis(2,mgp=c(2,0.45,0),tcl=-0.2,las=2,font=2,cex.axis=0.75,gap.axis=0)
  rec()
  
  par(mar=c(0.5,3,0,0))
  plot(climate$date,climate$prcp2,xlim=datelim,type="n",xaxt="n",yaxt="n",bty="n")
  recb()
  lines(climate$date,climate$prcp2,xlim=datelim,lwd=2,col="blue")
  axis(2,mgp=c(2,0.45,0),tcl=-0.2,las=2,font=2,cex.axis=0.75,gap.axis=0)
  abline(h=0,lty=3,lwd=1,col="blue")
  rec()
  
  par(mar=c(1.5,3,0.25,0))
  plot(xdate,0.75,pch=25,xlim=datelim,yaxt="n",xaxt="n",cex=2,bty="n",col=1,bg=1,ylim=c(0,max(xs2$sp))+1,log="y",type="n")
  recb
  at<-c(0.00000001,1,6,11,51,101,501,1001,5001,10001,max(xs2$sp)+1)
  axis(2,at=at,label=at-1,mgp=c(2,0.45,0),tcl=-0.2,las=2,font=2,cex.axis=0.75,gap.axis=0)
  labd<-paste0(yearpred,c("-05-01","-06-01","-07-01","-08-01","-09-01","-10-01","-11-01"))
  labd<-sort(c(labd,gsub("01","15",labd)))
  at<-as.Date(labd)
  lab<-format(at,"%b-%d")
  axis(1,at=at,labels=lab,mgp=c(1.5,0.45,0),font=2,tcl=-0.2)
  points(as.Date(xs2$date),xs2$sp+1,pch=16,cex=1,col=gray(0,0.1),xpd=TRUE)
  rec(ybottom=0.8,xpd=TRUE)

})
dev.off()
animation <- image_animate(img, fps = 2, optimize = TRUE)
image_write(animation,file.path("C:/Users/God/Downloads",paste0(paste0(spcode,yearpred),"predicted_abundance.gif")))
file.show(file.path("C:/Users/God/Downloads",paste0(paste0(spcode,yearpred),"predicted_abundance.gif")))


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

#par(mfrow=n2mfrow(length(v1m),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
#for(k in seq_along(v1m)){
p<-lapply((1:nsims)[1:20],function(i){
  dat<-xsmap@data
  juls<-lp[["jul"]][which.min(abs(lp[["jul"]]$jul-dat$jul[1])),,drop=FALSE] # finds the closest jul value to get the corresponding spline basis
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
  #spatial<-Amapmatrix %*% wk
  spatial<-arma_mm(Amapmatrix,wk) # ~ 3 times faster than %*%
  #}
  p<-fixed+spatial
  #p<-fixed # ignores spatial part
  #p<-spatial
  print(i)
  p
})
p<-do.call("cbind",p)
#p<-t(apply(p,1,function(i){c(quantile(i,0.0275,na.rm=TRUE),mean(i),quantile(i,0.975,na.rm=TRUE))}))
p<-cbind(rowQuantiles(p,probs=c(0.0275,0.975),na.rm=TRUE),rowMeans(p,na.rm=TRUE))[,c(1,3,2)]
#p<-exp(p)

xsmap$preds<-p[,2]

pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
pr<-mask(pr,mappingzone)
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



#### Map predictions across season with uncertainty #############################

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

#par(mfrow=n2mfrow(length(v1m),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
#for(k in seq_along(v1m)){

days<-seq(min(xs$jul[xs$year==yearpred]),max(xs$jul[xs$year==yearpred]),length.out=10)
lpr<-foreach(j=seq_along(days),.packages=c("raster")) %do% {
  juls<-lp[["jul"]][which.min(abs(lp[["jul"]]$jul-days[j])),,drop=FALSE]
  standardv<-names(nparams)[!names(nparams)%in%c("intercept","jul","julsquare",1:50)]
  dat<-as.matrix(xsmap@data[,standardv])
  dat<-cbind(dat,data.frame(jul=days[j],julsquare=days[j]^2)[rep(1,nrow(dat)),])
  dat<-cbind(dat,juls[,names(juls)%in%paste0("X",1:50)][rep(1,nrow(dat)),])
  dat<-cbind(intercept=1,dat)
  p<-lapply((1:nsims)[1:20],function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
    fixed<-as.matrix(dat[,names(betas)]) %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
    #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    #spatial<-Amapmatrix %*% wk
    spatial<-arma_mm(Amapmatrix,wk) # ~ 2 times faster than %*%
    #}
    p<-fixed+spatial
    #p<-fixed # ignores spatial part
    #p<-spatial
    cat("\r",paste(i,"/",max(nsims)," - ",j,"/",length(days)))
    p[,1]
  })
  p<-do.call("cbind",p)
  #p<-t(apply(p,1,function(i){c(quantile(i,0.0275,na.rm=TRUE),mean(i,na.rm=TRUE),quantile(i,0.975,na.rm=TRUE))}))
  p<-cbind(rowQuantiles(p,probs=c(0.0275,0.975),na.rm=TRUE),rowMeans(p,na.rm=TRUE))[,c(1,3,2)]
  #p<-exp(p)
  xsmap$preds<-p[,2]
  pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
  pr<-mask(pr,mappingzone)
  pr<-exp(pr)
  #pr<-disaggregate(pr,fact=1,method="bilinear")
  #print(j)
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
  lab<-niceround(exp(at))
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


#### Custom map predictions with posterior samples ############################

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

#days<-c("06-05","06-07","06-09","06-11")
days<-c("07-01")
days<-as.integer(format(as.Date(paste(yearpred,days,sep="-")),"%j"))
days<-(days-vscale$jul[["mean"]])/vscale$jul[["sd"]]

lprr<-foreach(j=seq_along(days),.packages=c("raster")) %do% {
  juls<-lp[["jul"]][which.min(abs(lp[["jul"]]$jul-days[j])),,drop=FALSE]
  standardv<-names(nparams)[!names(nparams)%in%c("intercept","jul","julsquare",1:50)]
  dat<-as.matrix(xsmap@data[,standardv])
  dat<-cbind(dat,data.frame(jul=days[j],julsquare=days[j]^2)[rep(1,nrow(dat)),])
  dat<-cbind(dat,juls[,names(juls)%in%paste0("X",1:50)][rep(1,nrow(dat)),])
  dat<-cbind(intercept=1,dat)
  p<-lapply((1:nsims)[1:200],function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-ifelse(names(nparams)%in%1:50,paste0("X",names(nparams)),names(nparams))
    fixed<-as.matrix(dat[,names(betas)]) %*% betas # make sure betas and vars are in the same order
    wk<-samples[[i]]$latent[nweights]
    #spatial<-Amapmatrix %*% wk
    spatial<-arma_mm(Amapmatrix,wk) # ~ 2 times faster than %*%
    #}
    p<-fixed+spatial
    #p<-fixed # ignores spatial part
    #p<-spatial
    cat("\r",paste(i,"/",max(nsims)," - ",j,"/",length(days)))
    p[,1]
  })
  p<-do.call("cbind",p)
  p<-cbind(rowQuantiles(p,probs=c(0.0275,0.975),na.rm=TRUE),rowMeans(p,na.rm=TRUE))[,c(1,3,2)]
  xsmap$preds<-p[,2]
  pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
  pr<-mask(pr,mappingzone)
  pr<-exp(pr)
  #pr<-disaggregate(pr,fact=1,method="bilinear")
  pr
}

lpr<-lapply(lprr,rast)


png(file.path("C:/Users/God/Downloads",paste0(spcode,"color_obs.png")),width=13,height=8,units="in",res=300)
par(mfrow=n2mfrow(length(days),asp=3.5/2),mar=c(0,0,3,0),oma=c(0,0,0,0))

jul<-round(days*vscale[["jul"]]["sd"]+vscale[["jul"]]["mean"],0)
datelim<-range(as.Date(format(as.Date(jul,origin=paste0(yearpred,"-01-01")),"%Y-%m-%d")))+c(-5,5)

zlim1<-range(sapply(lpr,function(i){range(values(i),na.rm=TRUE)}))
zlim2<-range(c(sapply(lpr,function(i){range(values(i),na.rm=TRUE)}),xs$sp[xs$year==yearpred]))
zlim<-c(zlim1[1],zlim2[2])

lapply(seq_along(lpr),function(i){
  j<-round(days[i]*vscale[["jul"]]["sd"]+vscale[["jul"]]["mean"],0)
  xdate<-as.Date(format(as.Date(j,origin=paste0(yearpred,"-01-01")),"%Y-%m-%d"))
  at<-seq(min(values(log(lpr[[i]])),na.rm=TRUE),max(values(log(lpr[[i]])),na.rm=TRUE),length.out=5)
  #lab<-ifelse(round(exp(at),0)==0,round(exp(at),2),round(exp(at),0))
  lab<-niceround(exp(at))
  labels<-paste(lab,c("min pred. > 0",rep("",length(at)-2),"max pred."))
  plot(log(lpr[[i]]),range=log(zlim),col=cols,asp=1,axes=FALSE,bty="n",plg=list(at=at,labels=labels,cex=0.75))
  plot(st_geometry(water),border=NA,col="white",add=TRUE)
  rd<-as.character(seq.Date(xdate-3,xdate+3,by=1))
  xxs<-st_transform(st_as_sf(xs[xs$date%in%rd,]),crs=crs(lpr[[1]]))
  colobs<-c(log(zlim),log(xxs$sp))
  colobs<-ifelse(is.infinite(colobs),NA,colobs)
  colobs<-colo.scale(colobs,cols)[-(1:2)]
  colobs<-ifelse(is.na(colobs),"#FFFFFF",colobs)
  #plot(st_geometry(xxs),pch=1,cex=scales::rescale(identity(c(0,max(xs$sp[xs$year==yearpred],na.rm=TRUE),xxs$sp)+0.01),to=c(0.25,10))[-(1:2)],add=TRUE)
  plot(st_geometry(xxs),pch=21,cex=1,add=TRUE,bg=colobs,col="grey10",lwd=0.4)
  #text(st_coordinates(st_geometry(xxs)),label=xxs$sp,cex=0.6,col="grey10",adj=c(0.5,-0.8))
  mtext(side=3,line=1,text=paste(gsub("_","",spcode),yearpred,sep="  "),adj=0.0)
  mtext(side=3,line=0,text=paste("predictions:",paste(format(xdate,"%b-%d"),collapse=" to "),sep="  "),adj=0.0)
  mtext(side=3,line=-1,text=paste("observations:",paste(format(as.Date(range(rd)),"%b-%d"),collapse=" to "),sep="  "),adj=0.0)

  mtext(side=4,line=-1,text="Number of mosquitos per trap (observed and predicted)",adj=0.5,cex=0.7)
  xp<-xmin(lpr[[1]])+((xmax(lpr[[1]])-xmin(lpr[[1]]))*c(0.85,0.85))
  yp<-ymin(lpr[[1]])+((ymax(lpr[[1]])-ymin(lpr[[1]]))*c(1.05,1.025)) 
  points(xp,yp,pch=21,cex=1,bg=colobs[c(which.min(xxs$sp),which.max(xxs$sp))],col="grey10",lwd=0.4,xpd=TRUE)
  #text(xp,yp,label=,cex=0.7,col="grey10",adj=c(0.5,-1))
  text(xp,yp,label=paste(c("min obs.","max obs."),xxs$sp[c(which.min(xxs$sp),which.max(xxs$sp))]),cex=0.7,col="grey10",adj=c(1.2,0.5),xpd=TRUE)

  
})
dev.off()
file.show(file.path("C:/Users/God/Downloads",paste0(spcode,"color_obs.png")))

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

ff<-function(x){
  if(FALSE){
    log(x+1)
  }else{
    identity(x)
  }
}

quant<-c("mean","0.5quant")[1]

par(mfrow=c(1,2))
plot(ff(m$summary.fitted.values[index[["est"]],quant]),ff(xs$sp),asp=1)
plot(ff(mfixed$summary.fitted.values[,quant]),ff(xs$sp),asp=1)
par(mfrow=c(1,1))
cor(ff(m$summary.fitted.values[index[["est"]],quant]),ff(xs$sp))^2
cor(ff(mfixed$summary.fitted.values[,quant]),ff(xs$sp))^2




#### SPDE posteriors ##########################################

res <- inla.spde.result(m, "spatial", spde)
par(mfrow=c(2,1))
plot(res$marginals.range.nominal[[1]],
     type="l", main="Posterior density for range")
plot(inla.tmarginal(sqrt, res$marginals.variance.nominal[[1]]),
     type="l", main="Posterior density for std.dev.")
par(mfrow=c(1,1))


### Combine maps ##########################################

images<-list.files("C:/Users/God/Downloads",pattern="*maps_all.png",full.names=TRUE)
ims<-do.call("c",lapply(images,function(x){
  im<-image_read(x)
  i<-image_info(im)
  print(i)
  im1<-image_crop(im,"x1200+0+0")
  im2<-image_crop(im,"x1200+1600+1200")
  res<-image_append(c(im1,im2),stack=FALSE)
  res<-image_draw(res)
  w<-image_info(res)$width
  h<-image_info(res)$height
  pos<-seq(0,w,by=w/5)
  te<-c("Mean","CI 2.5%","CI 97.5%","Mean Spatial Field","SD Spatial Field")
  for(i in seq_along(pos)){
    rect(0+pos[i],0,900+pos[i],150,col="white",border=NA,lty="dashed",lwd=5)
    if(x==images[1]){text(pos[i]+(w/5*0.4),50, te[i],family="Helvetica",cex=10,font=2)}
  }
  dev.off()
  border<-150
  res<-image_border(res,"#FFFFFF",border)
  res<-image_crop(res,paste0(image_info(res)$width-border,"x",image_info(res)$height,"+0+0"))
  print(image_info(res))
  code<-substr(sapply(strsplit(x,"/"),tail,1),1,3)
  res<-image_annotate(res,code,size=150,gravity="west",location="+70+200",degrees=270,weight=900)
  res
}))
im<-image_append(image_scale(ims, "6000"), stack = TRUE)
image_write(im,"C:/Users/God/Downloads/mosquito_maps.png")
file.show("C:/Users/God/Downloads/mosquito_maps.png")
im<-image_read("C:/Users/God/Downloads/mosquito_maps.png")
im<-image_scale(im[1],"x700")
image_write(im,"C:/Users/God/Downloads/reduced_mosquito_maps.png")
#file.show("C:/Users/God/Downloads/reduced_mosquito_maps.png")


### Combine marginal effects ######################################

images<-list.files("C:/Users/God/Downloads",pattern="*marginal_effects.png",full.names=TRUE)
ims<-do.call("c",lapply(images,function(x){
  im<-image_read(x)
  im<-image_border(im,"#FFFFFF","x100")
  code<-substr(sapply(strsplit(x,"/"),tail,1),1,3)
  im<-image_draw(im)
    text(350,50,code,family="Helvetica",cex=10,font=2)
  dev.off()
  im
}))
res1<-image_append(c(ims[1],ims[2]),stack=FALSE)
res2<-image_append(c(ims[3],ims[4]),stack=FALSE)  
res<-image_append(c(res1,res2),stack=TRUE)
image_write(res,"C:/Users/God/Downloads/mosquito_effects.png")
file.show("C:/Users/God/Downloads/mosquito_effects.png")
im<-image_read("C:/Users/God/Downloads/mosquito_effects.png")
im<-image_scale(im[1],"x700")
image_write(im,"C:/Users/God/Downloads/reduced_mosquito_effects.png")
#file.show("C:/Users/God/Downloads/reduced_mosquito_effects.png")


### Combine DIC tables ######################################

images<-list.files("C:/Users/God/Downloads",pattern="*dics.png",full.names=TRUE)
ims<-do.call("c",lapply(images,function(x){
  im<-image_read(x)
  im<-image_border(im,"#FFFFFF","25x90")
  code<-substr(sapply(strsplit(x,"/"),tail,1),1,3)
  im<-image_draw(im)
  text(110,50,code,family="Helvetica",cex=6,font=2)
  dev.off()
  im
}))
res1<-image_append(c(ims[1],ims[2]),stack=FALSE)
res2<-image_append(c(ims[3],ims[4]),stack=FALSE)  
res<-image_append(c(res1,res2),stack=TRUE)
image_write(res,"C:/Users/God/Downloads/mosquito_dic_tables.png")
file.show("C:/Users/God/Downloads/mosquito_dic_tables.png")
im<-image_read("C:/Users/God/Downloads/mosquito_dic_tables.png")
im<-image_scale(im[1],"x700")
image_write(im,"C:/Users/God/Downloads/reduced_mosquito_dic_tables.png")
#file.show("C:/Users/God/Downloads/reduced_mosquito_dic_tables.png")

### Combine posteriors ############################################

images<-list.files("C:/Users/God/Downloads",pattern="*posteriors.png",full.names=TRUE)
ims<-do.call("c",lapply(images,function(x){
  im<-image_read(x)
  im
}))
#res1<-image_append(c(ims[1],ims[2]),stack=FALSE)
#res2<-image_append(c(ims[3],ims[4]),stack=FALSE)  
res<-image_append(ims,stack=FALSE)
image_write(res,"C:/Users/God/Downloads/mosquito_posteriors_all.png")
file.show("C:/Users/God/Downloads/mosquito_posteriors_all.png")
im<-image_read("C:/Users/God/Downloads/mosquito_posteriors_all.png")
im<-image_scale(im[1],"x700")
image_write(im,"C:/Users/God/Downloads/reduced_mosquito_posteriors_all.png")
#file.show("C:/Users/God/Downloads/reduced_mosquito_dic_tables.png")


### Combine animations ######################################

images<-list.files("C:/Users/God/Downloads",pattern="*predicted_abundance.gif",full.names=TRUE)
ni<-nrow(image_info(image_read(images[1])))
#images<-images[c(1,1,1,1)]
ims<-lapply(images,function(x){
  im<-image_read(x)
  im<-image_border(im,"#FFFFFF","0x20")
  im
})
ims2<-lapply(1:ni,function(i){
  im<-do.call("c",lapply(seq_along(ims),function(j){
    ims[[j]][i]  
  }))  
  res1<-image_append(c(im[1],im[2]),stack=FALSE)
  res2<-image_append(c(im[3],im[4]),stack=FALSE)  
  res<-image_append(c(res1,res2),stack=TRUE)
  res
})
#image_montage(ims,tile="2x2")
animation<-image_animate(do.call("c",ims2),fps=2,optimize=FALSE)
image_write(animation,"C:/Users/God/Downloads/mosquito_animations.gif")
file.show("C:/Users/God/Downloads/mosquito_animations.gif")
rm(images,im,ni,ims,ims2,animation)
im<-image_read("C:/Users/God/Downloads/mosquito_animations.gif")
im<-image_scale(im[1],"x700")
image_write(im,"C:/Users/God/Downloads/reduced_mosquito_animations.gif")
#file.show("C:/Users/God/Downloads/reduced_mosquito_animations.gif")


### Study area map with lcc ################################################

# add a scale on map
lscale<-function(w=10000,lab="20 km"){
  wi<-diff(par("usr")[c(1,2)])
  he<-diff(par("usr")[c(3,4)])
  xleft<-par("usr")[1]+wi*0.5
  xright<-xleft+w
  h<-1500   
  ybottom<-par("usr")[3]+he*0.05
  ytop<-ybottom+h
  print(c(xleft,ybottom,xright,ytop))
  rect(xleft,ybottom,xright,ytop,col="white",border="black")
  rect(xright,ybottom,xright+w,ytop,col="black",border="black")
  text(xright+w,ybottom+(ytop-ybottom)/2,labels=lab,adj=c(-0.2,0.5),cex=0.55,font=2)
}

l<-list.files(file.path(path,"LULC/LULC/"),pattern=".tif",full.names=TRUE)
l<-l[substr(l,nchar(l)-3,nchar(l))==".tif"]
lulc<-stack(l)
water<-lulc[["LULC2011"]]
water<-rast(water)
water<-crop(water,st_transform(st_as_sf(mappingzone),crs=crs(water)))
water[water!=20]<-NA
water<-ms_simplify(st_union(st_transform(st_as_sf(as.polygons(water)),crs=crs(ds))),0.05)
z<-lulc[["LULC2011"]]
z<-mask(z,st_transform(st_as_sf(mappingzone),crs(z)))


can<-st_as_sf(raster::getData("GADM", country = "CAN", level = 1))
usa<-st_as_sf(raster::getData("GADM", country = "USA", level = 1))
ne<-rbind(can,usa)
plot(st_geometry(ne),axes=TRUE)
ne<-st_crop(ne,c(xmin=-100,ymin=30,xmax=-45,ymax=70))
ne<-st_make_valid(ne)
ne<-st_union(ne)
ne<-ms_simplify(ne,keep=0.03)


png("C:/Users/God/Downloads/location_map.png",width=8,height=8,units="in",res=200,pointsize=11)
lim<-c(-85,40,-55,55)
plot(st_geometry(ne),border=NA,col="grey90",axes=TRUE,bg="lightblue",xlim=lim[c(1,3)],ylim=lim[c(2,4)])
plot(lakes110,col="lightskyblue",border=NA,add=TRUE)
plot(st_geometry(ne),border=gray(0,0.1),col=NA,add=TRUE,lwd=0.1)
plot(st_geometry(st_transform(st_as_sf(mappingzone),st_crs(ne))),col=gray(0,0.15),border="red",lwd=0.5,add=TRUE)
mtext(side=3,line=0.25,text="Location of study area",font=2,cex=2,adj=0.025)
box(col="grey60")
dev.off()
file.show("C:/Users/God/Downloads/location_map.png")

lccnames$used<-c("other","water","other","urban","wetland","wetland","wetland","wetland","agriculture","agriculture","shrub","forest")
cols<-c(other="grey10",water="lightskyblue",urban="grey50",wetland="tan4",agriculture="khaki",shrub="yellowgreen",forest="forestgreen")
lccnames$cols<-adjustcolor(cols[match(lccnames$used,names(cols))],0.8)

png("C:/Users/God/Downloads/lcc_map.png",width=7,height=4,units="in",res=300,pointsize=11)
par(mar=c(0,0,0,8))
plot(z,col=lccnames$cols,breaks=c(0,lccnames$classn),axis.args=arg,xlim=xlim,ylim=ylim,zlim=c(0,220),legend=FALSE,legend.width=1,legend.shrink=1.25,axes=TRUE,bty="n")
locs<-ds[ds$id!="map" & !duplicated(ds@data[,c("longitude","latitude")]),]
plot(st_geometry(st_transform(st_as_sf(locs),proj4string(lulc))),cex=1,pch=16,col=adjustcolor("firebrick",0.6),add=TRUE)
name<-unique(lccnames$used)
legend("topright",legend=paste0(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name))),fill=unique(lccnames$cols),bty="n",border=NA,cex=1.2,xpd=TRUE,y.intersp=1.15,inset=c(-0.27,0))
legend("topright",legend=" Trap location",pch=16,col=adjustcolor("firebrick",0.6),bty="n",border=NA,cex=1.2,xpd=TRUE,inset=c(-0.305,0.44))
lscale()
box(col="white")
mtext(side=3,line=-1.75,text="Study area and land cover classes",adj=0.05,font=2,cex=1.25)
dev.off()
file.show("C:/Users/God/Downloads/lcc_map.png")

im1<-image_read("C:/Users/God/Downloads/location_map.png")
im2<-image_read("C:/Users/God/Downloads/lcc_map.png")
im<-image_composite(im2, image_scale(im1,"x565"),gravity="southeast",offset="+150+0")
#im3<-image_read("C:/Users/God/Downloads/trap_year.png")
#im<-image_append(c(im,image_scale(im3,"x1350")),stack=TRUE)
image_write(im,"C:/Users/God/Downloads/location.png")
file.show("C:/Users/God/Downloads/location.png")



### Show traps in mapping zone #######################

l<-split(ds[!ds$db%in%"map",],ds$year[!ds$db%in%"map"])

png("C:/Users/God/Downloads/trap_year.png",width=6,height=4,units="in",res=200,pointsize=11)
par(mfrow=n2mfrow(length(l)),mar=c(0,0,0,0),oma=c(0,0,0,0))
lapply(l,function(i){
  x<-i[!duplicated(i$longitude),]
  plot(mappingzone,col="grey90",border=NA)
  plot(water,col="lightblue",border=NA,add=TRUE)
  bb<-st_as_sfc(st_bbox(mappingzone))
  plot(st_difference(bb,st_as_sf(mappingzone)),col="white",border=NA,add=TRUE)
  plot(x,pch=16,cex=0.6,col=adjustcolor("firebrick",0.8),add=TRUE)
  mtext(i$year[1],side=3,adj=c(0.1),line=-1.5,cex=0.5,font=2)
})
dev.off()
file.show("C:/Users/God/Downloads/trap_year.png")


### Show no. of traps by week #####################

l<-split(ds[!ds$db%in%"map",],ds$week[!ds$db%in%"map"])
names(l)<-sapply(l,function(i){i$week[1]})
ee<-expand.grid(year=sort(unique(substr(names(l),1,4))),week=sort(unique(substr(names(l),6,8))))
m<-match(apply(ee,1,function(i){paste(i[1],i[2],sep="_")}),names(l))
ee$nbtraps<-sapply(m,function(i){if(is.na(i)){0}else{nrow(l[[i]])}})
ee<-ee[order(ee$year,ee$week),]
nbtraps<-split(ee,ee$year)

png("C:/Users/God/Downloads/trap_weeks.png",width=5,height=7,units="in",res=200,pointsize=11)
par(mfrow=c(length(nbtraps),1),mar=c(0.5,3,0,0),oma=c(3,0,0,0))
lapply(nbtraps,function(i){
  b<-barplot(i$nbtraps,names.arg=if(identical(i,nbtraps[[length(nbtraps)]])){i$week}else{NULL},ylim=c(0,max(ee$nbtraps)),border=NA,col="forestgreen",xpd=FALSE,yaxt="n")
  text(b[,1],rep(-15,nrow(i)),i$nbtraps,cex=0.75,xpd=TRUE)
  mtext(side=3,line=-1.5,text=i$year[1],adj=0.005)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = gray(0,0.05),border=NA)
  axis(2,las=TRUE)
  #grid()
})
dev.off()
file.show("C:/Users/God/Downloads/trap_weeks.png")


### check lp combinations for lcc vars

lpx<-lapply(names(vscale),function(i){
  bscale(lp[["agriculture50"]][,i],i)  
})

vs<-c("agriculture1000","forest1000","urban1000","water1000","shrub1000","barren1000","wetland1000")

hist(rowSums(do.call("cbind",lapply(vs,function(i){
  bscale(xs@data[,i],i)
}))),xlab="")



vs<-c("agriculture1000","forest1000","urban1000","water1000","shrub1000","barren1000","wetland1000")
#vs<-paste0(vs,".y")

hist(rowSums(ds@data[,vs]),xlab="",xlim=0:1,breaks=50)


### Compare with glmmTMB

library(patchwork)
library(ggplot2)
library(ggeffects)

load("CQP_model_outputs.RData")

dat<-xs@data[xs$db!="map",]
dat$y<-as.integer(dat$sp>0)

knots1<-knots
boundary<-NULL
knots<-seq(min(xs$jul)+0.50,max(xs$jul)-0.50,length.out=5)
boundary<-c(min(xs$jul)+0.10,max(xs$jul)-0.10)

mm<-glmmTMB(sp ~ 0 + ns(jul,knots=knots,Boundary=boundary) + agriculture50 + forest50 + agriculture1000 + forest1000 + anom2 + prcp2 + anom30 + prcp30 + offset(lognights), ziformula=~0, data=dat,family=nbinom2())

#mm<-glmmTMB(y ~ ns(jul,knots=knots) + agriculture50 + forest50 + agriculture1000 + forest1000 + anom2 + prcp2 + anom30 + prcp30 + offset(lognights), ziformula=~0, data=dat,family=binomial())

plot(dat$jul,dat$sp+1,log="y",pch=16,col=gray(0,0.05))
lines(lp$jul$jul,predict(mm,newdata=cbind(lp$jul,lognights=0),type="response")+1)
abline(v=knots,lty=3)
abline(v=boundary,lty=3,col="red")


#mm<-gam(sp ~ poly(jul,3) + agriculture50 + forest50 + agriculture1000 + forest1000 + anom2 + prcp2 + anom30 + prcp30 + offset(lognights), data=dat,family="nb")

vs<-all.vars(formula(mm))[-1]
gl<-lapply(vs,function(i){
  g<-ggpredict(mm,terms=paste(i,"[n=100]"))
  plot(g,add=TRUE,jitter=FALSE)+coord_cartesian(ylim=c(0,20)) 
})
wrap_plots(gl,nrow=2)



plot(pred[[1]])

plot(water,border=NA,col="lightblue")
plot(xs,cex=scales::rescale(xs$sp,to=c(0.5,10)),pch=ifelse(xs$sp==0,4,1),add=TRUE)




