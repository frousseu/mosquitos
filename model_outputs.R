

library(INLA)
library(splines)


load("VEX_model_selection.RData")

# select the best model from the list
model<-spmodels[[which.min(dics)]]

# the formula for the best model without the spatial effect
fixed<-formula(paste("y ~",as.character(model[[3]])[2]))

inla.setOption(inla.mode="experimental")

### The best model without the spatial effect ################################
#mfixed <- inla(fixed,data=inla.stack.data(stackest), 
#          control.predictor=list(compute=TRUE, A=inla.stack.A(stackest),link=1), 
#          #control.family=list(hyper=list(theta=prec.prior)), 
#          control.fixed=control.fixed,
#          control.inla = list(strategy='gaussian',int.strategy = "eb"),
#          num.threads="2:2",
#          verbose=TRUE,
#          control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
#          family="nbinomial")#"zeroinflatednbinomial1"


mfixed <- inla(fixed,data=cbind(y=xs$sp,xs@data), 
               control.predictor=list(compute=TRUE,link=1), 
               #control.family=list(hyper=list(theta=prec.prior)), 
               control.fixed=control.fixed,
               control.inla = list(strategy='auto',int.strategy = "eb"),
               num.threads="2:2",
               verbose=TRUE,
               control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
               family="nbinomial")#"zeroinflatednbinomial1"
  
  
# get all variables names
#xs@data<-cbind(xs@data,setNames(do.call("data.frame",lapply(1:35,function(i){rnorm(nrow(xs),0,1)})),paste0("yyy",1:35,"1000")))

lccs<-paste(unique(gsub("1000","",names(xs)[grep("1000",names(xs))])),collapse="|")
vss<-names(xs)[grep(paste0(lccs,"|anom|tmean|prcp"),names(xs))]
vss<-vss[!grepl("CQ",vss)]
vss<-vss[grep("anom|prcp|forest|agriculture|urban",vss)] # keep what is potentially in the model only and not other classes or tmean

fixedfull<-formula(paste0("y ~ -1 + ns(jul, knots = knots) + ",paste(vss,collapse="+")))

mfixedfull <- inla(fixedfull,data=cbind(y=xs$sp,xs@data), 
               control.predictor=list(compute=TRUE,link=1), 
               #control.family=list(hyper=list(theta=prec.prior)), 
               control.fixed=control.fixed,
               control.inla = list(strategy='auto',int.strategy = "eb"),
               num.threads="2:2",
               verbose=TRUE,
               control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
               family="nbinomial")#"zeroinflatednbinomial1"


### The best model with the spatial model ####################################
m <- inla(model,data=inla.stack.data(stackfull), 
          control.predictor=list(compute=TRUE, A=inla.stack.A(stackfull),link=1), 
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=control.fixed,
          control.inla = list(strategy='auto',int.strategy = "eb"),
          num.threads="2:2",
          verbose=TRUE,
          control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=TRUE),
          #control.mode = list(result = m, restart = TRUE)), # to rerun the model with NA predictions according to https://06373067248184934733.googlegroups.com/attach/2662ebf61b581/sub.R?part=0.1&view=1&vt=ANaJVrHTFUnDqSbj6WTkDo-b_TftcP-dVVwK9SxPo9jmPvDiK58BmG7DpDdb0Ek6xypsqmCSTLDV1rczoY6Acg_Zb0VRPn1w2vRj3vzHYaHT8JMCEihVLbY
          family="nbinomial")#"zeroinflatednbinomial1"


 #### Posterior samples ####################################

# from haakon bakk a, BTopic112
nsims<-500
samples<-inla.posterior.sample(nsims,m,num.threads="2:2")
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor" # not sure if should use APredictor or Predictor
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-do.call("cbind",samples.effect)
nbsize<-sapply(samples, function(x) x$hyperpar[grep("size",names(x$hyperpar))])
zeroprob<-sapply(samples, function(x) x$hyperpar[grep("probability",names(x$hyperpar))])

#### Hyperpar samples #####################################
# I don't think this works, hyperpars cannot be sampled independently of fixed effects and everything should be sampled jointly, this will overestimate variability
sampleshyper<-inla.hyperpar.sample(nsims,m)
nbsize<-sampleshyper[,grep("size",dimnames(sampleshyper)[[2]])]
zeroprob<-sampleshyper[,grep("probability",dimnames(sampleshyper)[[2]])]


#save.image("VEX_model_outputs.RData")
save.image(paste0(spcode,"model_outputs.RData"))

print("Done !")


