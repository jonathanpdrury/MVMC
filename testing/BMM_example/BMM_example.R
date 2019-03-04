require(ape)
require(phytools)
require(deSolve)

source('createModel_BM_MV.R')
source('PhenotypicModel_MV_BETA.R') 

set.seed(123)
tree = rcoal(3)

#first, for the case where there is no covariance between sigma1 and sigma2

sigma_1 = 0.5
sigma_2 = 0.2
sigma_12 = 0

rate.matrix<-matrix(c(sigma_1,sigma_12,sigma_12,sigma_2),nrow=2)

#the analytical solution to the multivariate variance covariance matrix for BM:

vcv_1<-kronecker(rate.matrix,vcv(tree))

#the "RPANDA" approach

bmm<-createModel_BM_MV(tree)

vcv_2<-getTipDistribution(bmm,params=c(0,0,log(sqrt(sigma_1)),log(sqrt(sigma_2)),sqrt(sigma_12)))$Sigma

##vcv_1 and vcv_2 are the same (ordered differently, but the same otherwise)



##now, for the case where there *is* covariance between sigma1 and sigma2


sigma_1 = 0.5
sigma_2 = 0.2
sigma_12 = 0.1

rate.matrix<-matrix(c(sigma_1,sigma_12,sigma_12,sigma_2),nrow=2)

#the analytical solution to the multivariate variance covariance matrix for BM:

vcv_1<-kronecker(rate.matrix,vcv(tree))

#the "RPANDA" approach

bmm<-createModel_BM_MV(tree)

vcv_2<-getTipDistribution(bmm,params=c(0,0,log(sqrt(sigma_1)),log(sqrt(sigma_2)),sqrt(sigma_12)))$Sigma

##vcv_1 and vcv_2 are no longer the same (the diagonal has changed in the RPANDA version)

