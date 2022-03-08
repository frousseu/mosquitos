

Amap<-inla.spde.make.A(mesh=mesh,loc=coordinates(xsmap),group=gs)


params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  grep(paste0(i,":"),row.names(samples[[1]]$latent))  
}) 
#table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
nweights<-grep("spatial",row.names(samples[[1]]$latent))[iset$spatial.group==unique(gs)]
Amapmatrix<-as.matrix(Amap)

#par(mfrow=n2mfrow(length(v1),asp=3.5/2),mar=c(3,2,1,1),oma=c(0,10,0,0))
#for(k in seq_along(v1)){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    fixed<-cbind(intercept=1,as.matrix(xsmap@data[,names(nparams[-1])])) %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
    #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    spatial<-Amapmatrix %*% rep(wk,unique(gs)) # stack was Apn in fire
    #}
    #p<-exp(fixed+spatial)
    p<-exp(fixed) # ignores spatial part
    p<-spatial
    print(i)
    p
  })
  p<-do.call("cbind",p)
  p<-t(apply(p,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  p

xsmap$preds<-p[,2]

pr<-rasterize(xsmap,pgrid,field="preds",fun=mean)
pr<-mask(pr,buf)
pr<-disaggregate(pr,fact=5,method="bilinear")

par(mfrow=c(1,2))
plot(pred[[1]])
plot(pr)




