.libPaths("/ddn/data/ckkr89/R/x86_64-pc-linux-gnu-library/3.4")
require(phytools)
require(MASS)
require(parallel)
require(snow)

source('Final_Simulation_Parameters.R')
source('sim_t_comp_bivariate.R')
source('Simulation_Code.R')

setwd("/ddn/data/ckkr89/Final Simulation Data")

run_sim_BM = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=sig2.matrices,
    root=root,
    Nsim=sim.vector,
    model="BM",
    save.values=TRUE
  )
}

run_sim_OU = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=sig2.matrices,
    root=root,
    pars.format=pars.format,
    half.lives=half.lives,
    OU.theta=OU.theta,
    Nsim=sim.vector,
    model="OU",
    save.values=TRUE
  )
}

run_sim_MC = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=sig2.matrices,
    root=root,
    pars.format=pars.list,
    half.lives=half.lives,
    Nsim=sim.vector,
    model="MC",
    save.values=TRUE
  )
}

run_sim_DDexp_neg = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=sig2.matrices,
    root=root,
    DD.tip.rate=DDneg.tip.values,
    Nsim=sim.vector,
    model="DDexp",
    save.values=TRUE
  )
}

run_sim_DDexp_pos = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=DDpos.sig2.matrices,
    root=root,
    DD.tip.rate=DDpos.tip.values,
    Nsim=sim.vector,
    model="DDexp",
    save.values=TRUE
  )
}

run_sim_DDlin_neg = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=sig2.matrices,
    root=root,
    DD.tip.rate=DDneg.tip.values,
    Nsim=sim.vector,
    model="DDlin",
    save.values=TRUE
  )
}

run_sim_DDlin_pos = function(sim.vector){
  mv_sim_multiple(
    tree.list=tree.list,
    sig2.matrices=DDpos.sig2.matrices,
    root=root,
    DD.tip.rate=DDpos.tip.values,
    Nsim=sim.vector,
    model="DDlin",
    save.values=TRUE
  )
}

no.cores = detectCores()
cl = makeCluster(no.cores, type=getClusterOption("type"))
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cl, ls())

clusterApply(cl,Nsim,run_sim_BM)
clusterApply(cl,Nsim,run_sim_OU)
clusterApply(cl,Nsim,run_sim_MC)
clusterApply(cl,Nsim,run_sim_DDexp_neg)
clusterApply(cl,Nsim,run_sim_DDexp_pos)
clusterApply(cl,Nsim,run_sim_DDlin_neg)
clusterApply(cl,Nsim,run_sim_DDlin_pos)

stopCluster(cl)


##lapply(tree.list,function(x)return(x[[sim.vector]]))














