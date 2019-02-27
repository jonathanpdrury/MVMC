require(phytools)
require(mvMORPH)
require(MASS)
require(msos)
require(parallel)
require(snow)

source('PhenotypicModel.R') #from RPANDA (Marc Manceau)
source('PhenotypicADiag.R') #from RPANDA (Marc Manceau)
source('DDexp_nogeo_ADiag.R') #from RPANDA 
source('DDlin_nogeo_ADiag.R') #from RPANDA
source('mv_BM_DD_Likelihood_Calculator.R')

load('DDexp_root_sig2_values.RData')
load('DDexp_r_term_matrices.RData')
load('tree_list.RData')
load('DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2.RData')


tree = tree.list[[2]][[1]]
sig2 = sig2.matrices[[2]]
r.term = pars.list[[2]][[2]][[2]]
sim.results = DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2_complete[[1]]


# trait.1 = c()
# trait.2 = c()
# for (m in 1:length(sim.results)){
#   trait.1 = c(trait.1,sim.results[[m]][1])
#   trait.2 = c(trait.2,sim.results[[m]][2])
# }
# 
# #double check: does the following keep the traits in the right order?
# traits<-data.frame(trait.1,trait.2)
# rownames(traits)<- names(sim.results)

traits = matrix(nrow=length(tree$tip.label),ncol=2)
for (j in 1:length(tree$tip.label)){
  traits[j,1] = sim.results[[j]][1]
  traits[j,2] = sim.results[[j]][2]
}
colnames(traits) = c("trait.1","trait.2")
rownames(traits) = names(sim.results)



##First off, fit data using mvMORPH

fit = mvBM(tree,traits,model="BM1")
# LogLikelihood: 	 -53.09887 

##this shoudl return mv-vcv matrix
kronecker(fit$sigma,vcv(tree))

##If we constrain the rate matrix to be 0 and pass the fit data from mvMORPH, should get the same likelihood

likelihood = log_likelihood_mv_BM_DD(
  tree = tree,
  sig2.matrix = fit$sigma,
  slope.matrix = matrix(rep(0,4),nrow=2),
  sim.value = sim.results,
  root = fit$theta,
  model = "DDexp"
)

# likelihood
# [,1]
# [1,] -140.7369


#Using 'model="BM"' should return the same result as mvMORPH, and as model = "DDexp" with a rate matrix of 0.

likelihood = log_likelihood_mv_BM_DD(
  tree = tree,
  sig2.matrix = fit$sigma,
  slope.matrix = matrix(rep(0,4),nrow=2),
  sim.value = sim.results,
  root = fit$theta,
  model = "BM"
)


# > likelihood
# [,1]
# [1,] -53.09887 <- identical to mvMORPH
