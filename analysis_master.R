## apply GBC/sGBC to AD proteomics dataset



## author: Ziyi Li (ziyi.li@emory.edu)
## 2/8/2018



## 1. save all the .R files to the same folder

## 2. load dataset
load("ADdata_300.rda")

## 3. source function files
source("GBC_EM.R")
source("DWL.R")

## 4. provide parameter settings
p = nrow(dat$X)
n = ncol(dat$X)
type = rep(0,p)  # all variables are assumed to follow normal distribution
param = rep(0.25,p)  # hyper parameter for variance of normal priors
L=5  # number of latent biclusters

## 5. nu1 and nu4 can be fixed
nu1=5
nu4=20

## 5.1 apply GBC to the dataset
res = GBC_EM(X=as.matrix(dat$X), L=L, dist.type=type, param=param, use.network=F, init="svd", nu1=nu1, nu4=nu4, cutoff=0)

## 5.2 apply sGBC to the dataset
res = GBC_EM(X=as.matrix(dat$X), L=L, dist.type=type, param=param, use.network=T, edge=dat$edge, init="svd", nu1=nu1, nu4=nu4, cutoff=0)


## 6.1 select nu1 and nu4 by BIC & apply GBC
source("BIC_tuning.R")

nu1_vec = c(7,9,11,13,15,20,25)
nu4_vec = c(20,40,50,60,70,90,110)

BIC_rec = matrix(999999, length(nu4_vec), length(nu1_vec))
for(i1 in 1:length(nu1_vec)){
     print(paste("i1=",i1))
     for(i2 in 1:length(nu4_vec)){
          print(paste("i2=",i2))
          set.seed(123)
          tmp = try(GBC_EM(X=as.matrix(dat$X), L=L, dist.type=dat$type, param=dat$param, use.network=FALSE, init="svd", nu1=nu1_vec[i1], nu4=nu4_vec[i2], cutoff=0),silent=T)
          if(is.list(tmp)){
               BIC_rec[i2,i1] = BIC_tuning(dat,tmp)
          }
     }
}
minind<-which(BIC_rec == min(BIC_rec), arr.ind=TRUE)[1,]
res = GBC_EM(X=as.matrix(dat$X), L=L, dist.type=type, param=param, use.network=F, init="svd", nu1=nu1_vec[minind[1]], nu4=nu4_vec[minind[2]], cutoff=0)
res$S # this is the bicluster result



## 6.2 select nu1 and nu4 by BIC & apply sGBC
source("BIC_tuning.R")

nu1_vec = c(7,9,11,13,15,20,25)
nu4_vec = c(20,40,50,60,70,90,110)

BIC_rec = matrix(999999, length(nu4_vec), length(nu1_vec))
for(i1 in 1:length(nu1_vec)){
     print(paste("i1=",i1))
     for(i2 in 1:length(nu4_vec)){
          print(paste("i2=",i2))
          set.seed(123)
          tmp = try(GBC_EM(X=as.matrix(dat$X), L=L, dist.type=dat$type, param=dat$param, use.network=T, edge=dat$edge, init="svd", nu1=nu1_vec[i1], nu4=nu4_vec[i2], cutoff=0),silent=T)
          if(is.list(tmp)){
               BIC_rec[i2,i1] = BIC_tuning(dat,tmp)
          }
     }
}
minind<-which(BIC_rec == min(BIC_rec), arr.ind=TRUE)[1,]
res = GBC_EM(X=as.matrix(dat$X), L=L, dist.type=type, param=param, use.network=F, init="svd", nu1=nu1_vec[minind[1]], nu4=nu4_vec[minind[2]], cutoff=0)
res$S # this is the bicluster result
