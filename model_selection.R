

library(INLA)
library(splines)

load("CQP_model_parameters.RData")

# Only the stack for estimation is used here
# This file is mean to be ran on a server so necessary packages need to be called

spmodels<-models[[gsub("_","",spcode)]]
dics<-numeric(length(spmodels))

inla.setOption(inla.mode="experimental")

for(i in seq_along(dics)){

#### Model ##################################################
m <- inla(spmodels[[i]],data=inla.stack.data(stackest), 
          control.predictor=list(compute=FALSE, A=inla.stack.A(stackest),link=1), 
          #control.family=list(hyper=list(theta=prec.prior)), 
          control.fixed=control.fixed,
          control.inla = list(strategy='gaussian',int.strategy = "eb"),
          num.threads="2:2",
          verbose=FALSE,
          control.compute=list(dic=TRUE,waic=FALSE,cpo=FALSE,config=FALSE),
          #control.mode = list(result = m, restart = TRUE)), # to rerun the model with NA predictions according to https://06373067248184934733.googlegroups.com/attach/2662ebf61b581/sub.R?part=0.1&view=1&vt=ANaJVrHTFUnDqSbj6WTkDo-b_TftcP-dVVwK9SxPo9jmPvDiK58BmG7DpDdb0Ek6xypsqmCSTLDV1rczoY6Acg_Zb0VRPn1w2vRj3vzHYaHT8JMCEihVLbY
          family="nbinomial")#"zeroinflatednbinomial1"

dics[i]<-m$dic$dic
print(paste(Sys.time(),"-",round(dics[i],1),"-",i,"/",length(dics)))

}

#save(dics,spmodels,file=paste0(spcode,"model_selection.RData"))
rm(m) # make sure no model is kept to reduce possibility of errors

save.image(paste0(spcode,"model_selection.RData"))

source("model_outputs.R")

print("Done !")
