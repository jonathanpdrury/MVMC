##setwd('D:/Documents/User/MBiol Project/Test/complete_sim_results')
#setwd('J:/My_Documents/Project/Simulation Study/Preliminary Simulation/complete_sim_results')
# source('Likelihood Test/PhenotypicModel.R') #from RPANDA (Marc Manceau)
# source('Likelihood Test/PhenotypicADiag.R') #from RPANDA (Marc Manceau)
# source('Likelihood Test/DDexp_nogeo_ADiag.R') #from RPANDA
# source('Likelihood Test/DDlin_nogeo_ADiag.R') #from RPANDA
log_likelihood_mv_BM_DD_BETA= function(
  tree,
  par = NULL, ##for optimisation, is a vector of all of the parameters [1:3] being the sig2 values, [4:5] being the slope
  sig2.matrix = NULL,
  slope.matrix = NULL,
  dd.ob = NULL,
  sim.value,
  model,
  optim = FALSE,
  return.anc = FALSE
  
){
  require(phytools)
  require(MASS)
  require(msos)
  
    
  
if (optim){
    sig2.matrix = matrix(c(exp(par[1]),par[3],par[3],exp(par[2])),ncol=2)
    if (model=="DDexp"||model=="DDlin"){
   		slope.matrix = matrix(c(par[4],0,0,par[5]),ncol=2)
   		no.taxa=length(sim.value)
  		test = sig2.matrix*exp(slope.matrix*no.taxa)
  		if(abs(test[1,2])>test[1,1] || abs(test[1,2])>test[2,2]){
 			 return(-1E6)
  		}
  		
    } else {
  		if(abs(sig2.matrix[1,2])>sig2.matrix[1,1] || abs(sig2.matrix[1,2])>sig2.matrix[2,2]){
 			 return(-1E6)
  		}
    
    }
  
}
  
  
  if (!model %in% c("BM","DDexp","DDlin")){
    stop("model must be BM, DDexp or DDlin")
  }

##there is an error here--what should this be? vector of 100? or 
  sim_values = c()
  sorted.sim = sim.value[tree$tip.label]
  for (i in 1:dim(sig2.matrix)[1]){
    for (j in 1:length(tree$tip.label)){
      position = j+((i-1)*length(tree$tip.label))
      sim_values[position] = sorted.sim[[j]][i]
    }
  }
  
  if (model=="BM"){
    vcv.tree = vcv(tree)
    vcv.sig2.product = kronecker(sig2.matrix,vcv.tree)
  } else {
  if(dd.ob@name!=model){stop("dd.ob is built for incorrect model")}    
    sign.sig12 = sign(sig2.matrix[1,2])
    av.sig12 = abs(sig2.matrix[1,2])
    
    block1 = getTipDistribution(dd.ob,params=c(0,log(sqrt(sig2.matrix[1,1])),slope.matrix[1,1]))$Sigma
    block2 = getTipDistribution(dd.ob,params=c(0,log(sqrt(av.sig12)),slope.matrix[1,2]))$Sigma
    block3 = getTipDistribution(dd.ob,params=c(0,log(sqrt(av.sig12)),slope.matrix[2,1]))$Sigma
    block4 = getTipDistribution(dd.ob,params=c(0,log(sqrt(sig2.matrix[2,2])),slope.matrix[2,2]))$Sigma
    
    block2 = sign.sig12*block2
    block3 = sign.sig12*block3
    
    #compile blocks
    vcv.sig2.product = rbind(cbind(block1[tree$tip.label,tree$tip.label],block2[tree$tip.label,tree$tip.label]),cbind(block3[tree$tip.label,tree$tip.label],block4[tree$tip.label,tree$tip.label]))
  }
  
D=kronecker(diag(2),matrix(1,ncol=1,nrow=length(tree$tip.label)))
  
op <- getOption("show.error.messages")
options(show.error.messages=FALSE)
IV = try(chol2inv(chol(vcv.sig2.product)))  
options(show.error.messages=op)

if (class(IV) == "try-error"){
  	
  	op <- getOption("show.error.messages")
 	options(show.error.messages=FALSE)
    IV = try(ginv(vcv.sig2.product))
    options(show.error.messages=op)

    if (class(IV) == "try-error"){
      
      return(-1E6)
      
    } else {
    
	  theta<-solve(t(D)%*%IV%*%D)%*%t(D)%*%IV%*%sim_values
	  if(return.anc){return(theta)} 
      data.root.difference = (sim_values - (D%*%theta))
  	  log.det.product = logdet(vcv.sig2.product)
      log.likelihood = -0.5*(t(data.root.difference)%*%IV%*%data.root.difference) - 0.5*log.det.product - length(tree$tip.label)*length(root)*log(2*pi)*0.5
    }
    
  } else {
  
	theta<-solve(t(D)%*%IV%*%D)%*%t(D)%*%IV%*%sim_values
	if(return.anc){return(theta)} 
    data.root.difference = (sim_values - (D%*%theta))
  	log.det.product = logdet(vcv.sig2.product)
    log.likelihood = -0.5*(t(data.root.difference)%*%IV%*%data.root.difference) - 0.5*log.det.product - length(tree$tip.label)*length(root)*log(2*pi)*0.5
 
  }
  return(log.likelihood)
}


