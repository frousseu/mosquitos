

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
          family="nbinomial")#"zeroinflatednbinomial1"


#### Posterior samples ####################################

# from haakon bakk a, BTopic112
nsims<-500
samples<-inla.posterior.sample(nsims,m,num.threads="1:1")
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


#save.image("mosquitos_model.RData")
#load("mosquitos_model.RData")