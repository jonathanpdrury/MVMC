##setwd('D:/Documents/User/MBiol Project/Test/complete_sim_results')
#setwd('J:/My_Documents/Project/Simulation Study/Preliminary Simulation/complete_sim_results')
# source('Likelihood Test/PhenotypicModel.R') #from RPANDA (Marc Manceau)
# source('Likelihood Test/PhenotypicADiag.R') #from RPANDA (Marc Manceau)
# source('Likelihood Test/DDexp_nogeo_ADiag.R') #from RPANDA 
# source('Likelihood Test/DDlin_nogeo_ADiag.R') #from RPANDA
log_likelihood_mv_BM_DD = function(
  tree,
  par = NULL, ##for optimisation, is a vector of all of the parameters [1:4] being the sig2 matrix, [5:6] being the root, [7:10] being the slope
  sig2.matrix = NULL,
  slope.matrix = NULL,
  root = NULL,
  sim.value,
  model,
  optim = FALSE
){
  require(phytools)
  require(MASS)
  require(msos)
  
  
  if (optim){
    sig2.matrix = matrix(par[1:4],ncol=2)
    root = par[5:6]
    if (model=="DDexp"||model=="DDlin"){
      slope.matrix = matrix(par[7:10],ncol=2)
    }
  }
  
  if (!model %in% c("BM","DDexp","DDlin")){
    stop("model must be BM, DDexp or DDlin")
  }

  sim_values = c()
  sorted.sim = sim.value[tree$tip.label]
  for (i in 1:length(root)){
    for (j in 1:length(tree$tip.label)){
      position = j+((i-1)*length(tree$tip.label))
      sim_values[position] = sorted.sim[[j]][i]
    }
  }
  
  if (model=="BM"){
    vcv.tree = vcv(tree)
    vcv.sig2.product = kronecker(sig2.matrix,vcv.tree)
  } else {
    if (model=="DDexp"){	
      dd.ob = createModel_DDexp(tree)
    } else if (model=="DDlin") {
      dd.ob = createModel_DDlin(tree)
    }
    
    sign.sig12 = sign(sig2.matrix[1,2])
    av.sig12 = abs(sig2.matrix[1,2])
    
    block1 = getTipDistribution(dd.ob,params=c(0,log(sqrt(sig2.matrix[1,1])),slope.matrix[1,1]))$Sigma
    block2 = getTipDistribution(dd.ob,params=c(0,log(sqrt(av.sig12)),slope.matrix[1,2]))$Sigma
    block3 = getTipDistribution(dd.ob,params=c(0,log(sqrt(av.sig12)),slope.matrix[2,1]))$Sigma
    block4 = getTipDistribution(dd.ob,params=c(0,log(sqrt(sig2.matrix[2,2])),slope.matrix[2,2]))$Sigma
    
    block2 = sign.sig12*block2
    block3 = sign.sig12*block3
    
    #compile blocks
    vcv.sig2.product = rbind(cbind(block1,block2),cbind(block3,block4))
  }
  
  D=kronecker(diag(2),matrix(1,ncol=1,nrow=length(tree$tip.label)))
  log.det.product = logdet(vcv.sig2.product)
  data.root.difference = (sim_values - (D%*%matrix(root)))
  
  #log.likelihood = -0.5*(t(data.root.difference)%*%ginv(vcv.sig2.product)%*%data.root.difference) - 0.5*log.det.product - length(tree$tip.label)*length(root)*log(2*pi)*0.5
  test = try(chol2inv(chol(vcv.sig2.product)))
  if (class(test) == "try-error"){
    test = try(ginv(vcv.sig2.product))
    if (class(test) == "try-error"){
      return(1E6)
    } else {
      log.likelihood = -0.5*(t(data.root.difference)%*%ginv(vcv.sig2.product)%*%data.root.difference) - 0.5*log.det.product - length(tree$tip.label)*length(root)*log(2*pi)*0.5
    }
  } else {
    log.likelihood = -0.5*(t(data.root.difference)%*%chol2inv(chol(vcv.sig2.product))%*%data.root.difference) - 0.5*log.det.product - length(tree$tip.label)*length(root)*log(2*pi)*0.5
    if (optim){
      return(-log.likelihood)
    } else {
      return(log.likelihood)
    }
  }
}






