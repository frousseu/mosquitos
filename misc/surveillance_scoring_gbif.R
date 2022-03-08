
library(surveillance)

### Equation in Czado et al. 2009 https://doi.org/10.1111/j.1541-0420.2009.01191.x
rpsmc<-function(x,mu,size=NULL,n=10000){
  sapply(seq_along(mu),function(i){
    X<-rnbinom(n,mu=mu[i],size=size)
    Xp<-rnbinom(n,mu=mu[i],size=size)
    mean(abs(X-x[i]))-((1/2)*(mean(abs(X-Xp))))
  })
}

### Scoring small mu
n<-50
size<-0.5
x<-rep(5,n)
mu<-seq(0,20,length.out=n)
scores1<-scores(x,mu,size=size,which="rps")  
scores2<-rpsmc(x,mu,size=size,n=100000)
plot(0,0,type="n",xlim=range(mu),ylim=range(c(scores1,scores2)),yaxs="i",xaxs="i",xlab="mu",ylab="score")
lines(mu,scores1,col="red",lwd=2)
lines(mu,scores2,col="blue",lwd=2)
abline(v=x[1],lty=3)


### Comparing scores and the Czado equation for random values
size<-0.5
n<-50
iter<-25
l<-lapply(1:iter,function(i){
  x<-rnbinom(n,mu=5,size=size)
  mu<-rgamma(n,0.1,0.01)
  c(sum(scores(x,mu,size=size,which="rps")),sum(rpsmc(x,mu,size=size,n=10000)))
})
m<-do.call("rbind",l)
matplot(m,pch=c(1,3),xlab="iter",ylab="score")


library(rgbif)

name_lookup(query="Echinochloa colona")$data[,1:2] |> as.data.frame()
name_backbone(name="Echinochloa colona") |> as.data.frame()
name_suggest(q="Echinochloa colona")
occ_search(scientificName="Echinochloa colona",limit=200)

name_lookup(query="Cyperus involucratus")$data[,1:2] |> as.data.frame()
name_backbone(name="Cyperus involucratus") |> as.data.frame()
name_suggest(q="Cyperus involucratus")
occ_search(scientificName="Cyperus involucratus",limit=200)
occ_search(search="Cyperus involucratus",limit=200)

name_lookup(query="Eulalia aurea")$data[,1:2] |> as.data.frame()
name_backbone(name="Eulalia aurea") |> as.data.frame()
name_suggest(q="Eulalia aurea")
occ_search(scientificName="Eulalia aurea",limit=200)
occ_search(search="Eulalia aurea",limit=200)

name_lookup(query="Paspalum conjugatum")$data[,1:2] |> as.data.frame()
name_backbone(name="Paspalum conjugatuma") |> as.data.frame()
name_suggest(q="Paspalum conjugatum")
occ_search(scientificName="Paspalum conjugatum",limit=200)
occ_search(search="Paspalum conjugatum",limit=200)


a<-c(0,90,180,270,360)
(-a+90)%%360



data(SPDEtoy)
coords <- as.matrix(SPDEtoy[, 1:2])
p5 <- coords[1:5, ]

# Boundaries
bound1 <- inla.nonconvex.hull(p5)
bound2 <- inla.nonconvex.hull(p5, convex = 0.5, concave = -0.15)
bound3 <- inla.nonconvex.hull(p5, concave = 0.5)
bound4 <- inla.nonconvex.hull(p5, concave = 0.5,
                              resolution = c(20, 20))

# Meshes
m10 <- inla.mesh.2d(boundary = bound1, cutoff = 0.05, 
                    max.edge = c(0.1, 0.2))
m11 <- inla.mesh.2d(boundary = bound2, cutoff = 0.05, 
                    max.edge = c(0.1, 0.2))
m12 <- inla.mesh.2d(boundary = bound3, cutoff = 0.05, 
                    max.edge = c(0.1, 0.2))
m13 <- inla.mesh.2d(boundary = bound4, cutoff = 0.05, 
                    max.edge = c(0.1, 0.2))

par(mfrow=c(2,2))
plot(m10,asp=1)
plot(m11,asp=1)
plot(m12,asp=1)
plot(m13,asp=1)

alt<-raster("C:/Users/God/Downloads/srtm_48_17/srtm_48_17.tif")
run<-st_read("C:/Users/God/Downloads","Reunion_2015_region")
run<-st_transform(run,st_crs(alt))

alt<-crop(alt,run)

par(bg="black")
plot(st_geometry(run))
plot(alt,col=gray(rev(seq(0.3,0.9999,by=0.01))),add=TRUE)


slope <- terrain(alt, opt='slope')
aspect <- terrain(alt, opt='aspect')
hill <- hillShade(slope, aspect, 40,270)
plot(hill, col=grey(0:100/100), legend=FALSE, main='Switzerland')


matrix(1:24,ncol=6) %*% 1:3







