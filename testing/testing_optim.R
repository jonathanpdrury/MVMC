require(phytools)
require(MASS)
require(msos)

# setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/Likelihood Test')
source('PhenotypicModel.R') #from RPANDA (Marc Manceau)
source('PhenotypicADiag.R') #from RPANDA (Marc Manceau)
source('DDexp_nogeo_ADiag.R') #from RPANDA 
source('DDlin_nogeo_ADiag.R') #from RPANDA
source('mv_BM_DD_Likelihood_Calculator.R')

# setwd('../complete_DDexp_neg')
load('DDexp_root_sig2_values.RData')
load('DDexp_r_term_matrices.RData')
load('tree_list.RData')
load('DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2.RData')

tree = tree.list[[2]][[1]]
sig2 = sig2.matrices[[2]]
r.term = pars.list[[2]][[2]][[2]]
sim.results = DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2_complete[[1]]

likelihood = log_likelihood_mv_BM_DD(
  tree = tree,
  sig2.matrix = sig2,
  slope.matrix = r.term,
  root = c(0,0),
  sim.value = sim.results,
  model = "DDexp"
)

# likelihood
# [,1]
# [1,] -276.5095

par = c(log(sig2[1]),log(sig2[4]),sig2[2],0,0,r.term[1],r.term[4],r.term[2])

likelihood_optim = log_likelihood_mv_BM_DD(
  tree = tree,
  par = par,
  sim.value = sim.results,
  model = "DDexp",
  optim = TRUE
)

# likelihood_optim
# [,1]
# [1,] 276.5095

trait.1 = c()
trait.2 = c()
for (m in 1:length(sim.results)){
  trait.1 = c(trait.1,sim.results[[m]][1])
  trait.2 = c(trait.2,sim.results[[m]][2])
}

sig2.1var = log(var(trait.1)/max(nodeHeights(tree)))
sig2.2var = log(var(trait.2)/max(nodeHeights(tree)))

max_likelihood = optim(
  par = c(sig2.1var,sig2.2var,0,0,0,0,0,0),
  fn = log_likelihood_mv_BM_DD,
  tree = tree,
  sim.value = sim.results,
  model = "DDexp",
  optim = TRUE,
  control = list(maxit = 2000, fnscale = -1)
)

# max_likelihood
# $`par`
# [1] -3.251883553 -2.130652207  0.124717241  0.036555756  0.051304366 -0.006800662  0.001726189 -0.015084684
# 
# $value
# [1] 535835.1
# 
# $counts
# function gradient 
# 1985       NA 
# 
# $convergence
# [1] 10
# 
# $message
# NULL

max_likelihood_BFGS = optim(
  par = c(sig2.1var,sig2.2var,0,0,0,0,0,0),
  fn = log_likelihood_mv_BM_DD,
  tree = tree,
  sim.value = sim.results,
  model = "DDexp",
  optim = TRUE,
  control = list(maxit = 2000, fnscale = -1),
  method = "BFGS"
)

# max_likelihood_BFGS
# $`par`
# [1]  -1.7375020   0.2301150   2.9998500   0.3952839   0.1369419  59.0968320 106.7446311   0.0000000
# 
# $value
# [1] 1e-50
# 
# $counts
# function gradient 
# 4        2 
# 
# $convergence
# [1] 0
# 
# $message
# NULL

max_likelihood_CG = optim(
  par = c(sig2.1var,sig2.2var,0,0,0,0,0,0),
  fn = log_likelihood_mv_BM_DD,
  tree = tree,
  sim.value = sim.results,
  model = "DDexp",
  optim = TRUE,
  control = list(maxit = 2000),
  method = "CG"
)

# max_likelihood_CG
# $`par`
# [1] -4.940483e+00 -4.666197e+00 -2.286461e+03 -4.113091e-01 -1.410137e-01 -6.166320e+01 -1.113802e+02  3.712906e+05
# 
# $value
# [1] 1e-50
# 
# $counts
# function gradient 
# 51       17 
# 
# $convergence
# [1] 0
# 
# $message
# NULL

max_likelihood_SANN = optim(
  par = c(sig2.1var,sig2.2var,0,0,0,0,0,0),
  fn = log_likelihood_mv_BM_DD,
  tree = tree,
  sim.value = sim.results,
  model = "DDexp",
  optim = TRUE,
  control = list(maxit = 2000, fnscale = -1),
  method = "SANN"
)

# max_likelihood_SANN
# $`par`
# [1] -4.134137353 -1.409975988  0.001626135  2.241777577 -0.090846656 -0.437061172 -1.849721964 -0.402537725
# 
# $value
# [1] 238565030
# 
# $counts
# function gradient 
# 2000       NA 
# 
# $convergence
# [1] 0
# 
# $message
# NULL
# 




