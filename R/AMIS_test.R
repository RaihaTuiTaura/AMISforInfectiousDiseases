#' Define simple "transmission" model
#'
#' Uses the negative binomial model of worm burden W and overdispersion k.
#' p = 1 - (1+W/k)^-k.
#' See Anderson and May (1992).
#' param seeds a vector of seeds
#' param parameters a matrix of sampled parameter vectors
#' return a vector of prevalences 
transmission_model_simple<-function(seeds,parameters) {
  n_samples<-length(seeds)
  prevalences <- rep(NA,n_samples)
  for (i in 1:n_samples) {
    prevalences[i]<-1 - (1+parameters[i,1]/parameters[i,2])^(-parameters[i,2])
  }
  return(prevalences)
}
# prevalence equals first parameter
transmission_model_identity<-function(seeds,parameters) {
  return(matrix(parameters[,1],length(parameters[,1]),1))
}
# number of locations
L<-3
# number of map samples
M<-1000
prevalence_map<-matrix(NA,L,M)
# produce samples for prevalence map with 3 locations given by B(2,1), B(1,1)=Uniform, B(1,2). 
for (l in 1:L) {
  prevalence_map[l,]<-rbeta(M,max(1,l-1),max(1,3-l))
}
rownames(prevalence_map)<-c("Here","There","Everywhere")
# 2D exponential prior
rprior <- function(n) {
  params<-matrix(NA,n,2)
  colnames(params)<-c("W","k")
  params[,1]<-rexp(n)
  params[,2]<-rexp(n)
  return(params)
}
dprior <- function(x,log=FALSE) {
  if (log) {
    return(sum(dexp(x,log=T)))
  } else {
    return(prod(dexp(x)))
  }
}
prior<-list(rprior=rprior,dprior=dprior)
amis_params<-default_amis_params()
# check log works
amis_params$log<-TRUE
# output using log
logout<-amis(prevalence_map,transmission_model_identity, prior, amis_params,seed=1)
library(weights)
par(mfcol=c(2,3))
x<-0:100/100
# first row of plots should look like samples from the 3 beta distributions. second row should not.
for (i in 1:L) {
  wtd.hist(logout[,2],breaks=100,weight=exp(logout[,4+i]),xlim=c(0,1.5),probability=T)
  lines(x,dbeta(x,max(1,i-1),max(1,3-i)),col="red")
  wtd.hist(logout[,3],breaks=100,weight=exp(logout[,4+i]),xlim=c(0,1.5),probability=T)
}
# compare with unlogged
amis_params<-default_amis_params()
# should achieve identical plots 
output<-amis(prevalence_map,transmission_model_identity, prior, amis_params,seed=1)
par(mfcol=c(2,3))
for (i in 1:L) {
  wtd.hist(output[,2],breaks=100,weight=output[,4+i],xlim=c(0,1.5))
  wtd.hist(output[,3],breaks=100,weight=output[,4+i],xlim=c(0,1.5))
}
# should  be identical to numerical error
summary(exp(logout[,5:7])-output[,5:7])

# Now test in 1D
rprior1D <- function(n) {
  params<-matrix(NA,n,1)
  colnames(params)<-c("W")
  params[,1]<-rexp(n)
  return(params)
}
dprior1D <- function(x,log=FALSE) {
  if (log) {
    return(dexp(x,log=T))
  } else {
    return(dexp(x))
  }
}
prior1D<-list(rprior=rprior1D,dprior=dprior1D)
# should produce similar plots
outone<-amis(prevalence_map,transmission_model_identity,prior1D, amis_params,seed=1)
par(mfcol=c(1,3))
for (i in 1:L) {
  wtd.hist(outone[,2],breaks=100,weight=outone[,3+i],xlim=c(0,1.5))
}
#
# Now try analytical version of same map
#
likelihood<-function(data,prevalence,log) {
  if (is.matrix(data)) {
    lik<-matrix(NA,length(prevalence),length(data[,1]))
    for (i in 1:length(data[,1])) {
      lik[,i]<-dbeta(prevalence,data[i,1],data[i,2],log=log)
    }
    return(lik)
  } else {
    return(dbeta(prevalence,data[1],data[2],log=log))
  }
}
# analytical version of same map.
an_prevalence_map<-list(list(data=matrix(c(1,1,2,2,1,1),3,2),likelihood=likelihood))
# first row of plots should look even closer to beta distributions
outan<-amis(an_prevalence_map,transmission_model_identity,prior,amis_params,seed=1)
par(mfcol=c(2,3))
for (i in 1:L) {
  wtd.hist(outan[,2],breaks=100,weight=outan[,4+i],xlim=c(0,1.5))
  wtd.hist(outan[,3],breaks=100,weight=outan[,4+i],xlim=c(0,1.5))
}
#
# Check case where there is only 1 location
#
an_one_prevalence_map<-list(list(data=matrix(c(1,2),1,2),likelihood=likelihood))
outan_one<-amis(an_one_prevalence_map,transmission_model_identity,prior,amis_params,seed=1)
# With this method, first row of plots should be IDENTICAL to beta distributions, although this is an illusion!
par(mfcol=c(1,2))
wtd.hist(outan_one[,2],breaks=0:100/100,weight=outan_one[,5],xlim=c(0,1))
wtd.hist(outan_one[,3],breaks=100,weight=outan_one[,5],xlim=c(0,1.5))
#
# Now try histogram version
#
brs<-c(0:100/100)
amis_params[["breaks"]]<-brs
outanh<-amis(an_prevalence_map,transmission_model_identity,prior,amis_params,seed=1)
# With this method, first row of plots should be IDENTICAL to beta distributions, although this is an illusion!
par(mfcol=c(2,3))
for (i in 1:L) {
  wtd.hist(outanh[,2],breaks=0:100/100,weight=outanh[,4+i],xlim=c(0,1))
  wtd.hist(outanh[,3],breaks=100,weight=outanh[,4+i],xlim=c(0,1.5))
}
#
# Try final combo - histogram weights with sample-based map
#
outhut<-amis(prevalence_map,transmission_model_identity,prior,amis_params,seed=1)
# first row of plots should look like samples from beta distributions.
par(mfcol=c(2,3))
for (i in 1:L) {
  wtd.hist(outhut[,2],breaks=100,weight=outhut[,4+i],xlim=c(0,1.5))
  wtd.hist(outhut[,3],breaks=100,weight=outhut[,4+i],xlim=c(0,1.5))
}
#
# Now investigate multiple timepoints
#
# this function uses parameter 1 as prevalence at timepoint 1 and so on
transmission_model_identity<-function(seeds,parameters) {
  return(parameters)
}
# don't do histogram weighting
amis_params[["breaks"]]<-NULL
# same analytical map for timepoints 1 and 2.
map_list<-list(an_prevalence_map[[1]],an_prevalence_map[[1]])
multid<-amis(map_list,transmission_model_identity,prior,amis_params,seed=1)
# second row is filtering distribution (should look like 3 betas)
# first row is not smoothing distribution
par(mfcol=c(2,3))
for (i in 1:L) {
  wtd.hist(multid[,2],breaks=100,weight=multid[,5+i],xlim=c(0,1))
  wtd.hist(multid[,3],breaks=100,weight=multid[,5+i],xlim=c(0,1))
}
#
# Check multiple timepoints with histogram weights
#
amis_params[["breaks"]]<-0:100/100
multidh<-amis(map_list,transmission_model_identity,prior,amis_params,seed=1)
par(mfcol=c(2,3))
#plots should be same as previous, but less noisy
for (i in 1:L) {
  wtd.hist(multidh[,2],breaks=100,weight=multidh[,5+i],xlim=c(0,1.5))
  wtd.hist(multidh[,3],breaks=100,weight=multidh[,5+i],xlim=c(0,1.5))
}
#
# Try out Bayesian updating
#
# 2D exponential prior
rprior <- function(n) {
  params<-matrix(NA,n,5)
  colnames(params)<-c("A","B","C","D","E")
  params[,1]<-rexp(n)
  params[,2]<-rexp(n)
  params[,3]<-rexp(n)
  params[,4]<-rexp(n)
  params[,5]<-rexp(n)
  return(params)
}
dprior <- function(x,log=FALSE) {
  if (log) {
    return(sum(dexp(x,log=T)))
  } else {
    return(prod(dexp(x)))
  }
}
prior<-list(rprior=rprior,dprior=dprior)
map_list<-list(an_prevalence_map[[1]],an_prevalence_map[[1]],an_prevalence_map[[1]],an_prevalence_map[[1]],an_prevalence_map[[1]])
map_list[[3]]$data<-matrix(c(1,1,NA,2,1,NA),3,2)
amis_params[["bayesian"]]<-c(TRUE,rep(TRUE,length(map_list)-1))
multibu<-amis(map_list,transmission_model_identity,prior,amis_params,seed=1)
par(mfrow=c(3,5))
#plots should not be same as target, because of influence of prior
for (i in 1:L) {
  wtd.hist(multibu[,2],breaks=0:100/100,weight=multibu[,11+i],xlim=c(0,1))
  wtd.hist(multibu[,3],breaks=0:100/100,weight=multibu[,11+i],xlim=c(0,1))
  wtd.hist(multibu[,4],breaks=100,weight=multibu[,11+i],xlim=c(0,1.5))
  wtd.hist(multibu[,5],breaks=0:100/100,weight=multibu[,11+i],xlim=c(0,1))
  wtd.hist(multibu[,6],breaks=0:100/100,weight=multibu[,11+i],xlim=c(0,1))
}
