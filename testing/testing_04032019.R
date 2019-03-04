###first fitting according to mvBM example

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Simulate the traits
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(0,0)
data<-mvSIM(tree, param=list(sigma=sigma, ntraits=2, theta=theta,
            names_traits=c("head.size","mouth.size")), model="BM1", nsim=1)

## Fitting the models
# BMM - Analysis with multiple rates
fit_mvBM<-mvBM(tree, data,model="BM1")



bmm<-createModel_BM_MV(tree)

vcv_2<-getTipDistribution(bmm,params=c(0,0,log(sqrt(sigma_1)),log(sqrt(sigma_2)),sqrt(sigma_12)))$Sigma


data.sorted<-c(data[rownames(vcv_2)[1:50],1],data[rownames(vcv_2)[1:50],2])

fit_RPANDA<-fitTipData(bmm,data.sorted)


####appears to be working well (though could update as indicated in Marc email)

load('~/MVMC/testing/DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2.RData')
source('~/MVMC/R/mv_BM_DD_Likelihood_Calculator.R', chdir = TRUE)
source('~/MVMC/R/PhenotypicModel.R') #from RPANDA (Marc Manceau)
source('~/MVMC/R/PhenotypicADiag.R') #from RPANDA (Marc Manceau)
source('~/MVMC/R/DDexp_nogeo_ADiag.R') #from RPANDA 
source('~/MVMC/R/DDlin_nogeo_ADiag.R') #from RPANDA
load('~/MVMC/testing/tree_list.RData')
require(phytools)

sim.results = DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2_complete[[1]]

tree = tree.list[[2]][[1]]
sig2 = sig2.matrices[[2]]
r.term = pars.list[[2]][[2]][[2]]

trait.1 = c()
trait.2 = c()
for (m in 1:length(sim.results)){
  trait.1 = c(trait.1,sim.results[[m]][1])
  trait.2 = c(trait.2,sim.results[[m]][2])
}

sig2.1var = log(var(trait.1)/max(nodeHeights(tree)))
sig2.2var = log(var(trait.2)/max(nodeHeights(tree)))

dd.ob<-createModel_DDexp(tree)

max_likelihood = optim(
 par = c(sig2.1var,sig2.2var,0,0,0,0,0,0),
  fn = log_likelihood_mv_BM_DD,
  tree = tree,
  sim.value = sim.results,
  dd.ob = dd.ob,
  model = "DDexp",
  optim = TRUE,
  control = list(maxit = 2000, fnscale = -1)
)

#return the root:

log_likelihood_mv_BM_DD(tree,par=max_likelihood$par,return.anc=TRUE,model="DDexp",dd.ob=dd.ob,sim.value=sim.value,optim=TRUE)


source('~/MVMC/R/PhenotypicModel_MV_BETA.R')
source('~/MVMC/R/createModel_BM_MV.R')
source('~/MVMC/R/createModel_DDexp_MV.R')

require(deSolve)

#dumb, inefficient hack to get ordering right
bmm<-createModel_BM_MV(tree)
vcv_2<-getTipDistribution(bmm,params=c(0,0,0,0,0))$Sigma

data = matrix(nrow=length(tree$tip.label),ncol=2)
for (j in 1:length(tree$tip.label)){
  data[j,1] = sim.results[[j]][1]
  data[j,2] = sim.results[[j]][2]
}
colnames(data) = c("trait.1","trait.2")
rownames(data) = names(sim.results)

data.sorted<-c(data[rownames(vcv_2)[1:50],1],data[rownames(vcv_2)[1:50],2])


ddm.ob<-createModel_DDexp_MV(tree)

fit_RPANDA_DD<-fitTipData(ddm.ob,data.sorted)
