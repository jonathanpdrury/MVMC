log_likelihood_univariate_BM = function(tree,sig2,root,sim){
  require(phytools)
  require(matlib)
  
  vcv.tree = vcv(tree)
  column.vector = matrix(rep.int(1,length(tree$tip.label)))
  vcv.sig2 = sig2 * vcv.tree
  det.vcv.sig2 = det(vcv.sig2)
  data.root.difference = (sim - root*column.vector)
  
  exponent = exp(-1/(2*(t(data.root.difference)%*%inv(vcv.sig2)%*%data.root.difference)))
  denominator = sqrt(((2*pi)^length(tree$tip.label))*det.vcv.sig2)
  log.likelihood = log10(exponent / denominator)
  return(log.likelihood)
}

log_likelihood_multivariate_BM = function(tree,sig2,root,sim){
  require(phytools)
  require(matlib)
 
##should be sorted by tree$tip.label  
  sim_values = c()
  for (i in 1:length(root)){
    for (j in 1:length(tree$tip.label)){
      position = j+((i-1)*length(tree$tip.label))
      sim_values[position] = sim[[j]][i]
    }
  }
  
  D = matrix(ncol=length(root),nrow=(length(root)*length(tree$tip.label)))
  for (k in 1:nrow(D)){
    for (l in 1:ncol(D)){
      if ((l-1)*length(tree$tip.label) < k && k <= l * length(tree$tip.label)){
        D[k,l] = 1
      } else {
        D[k,l] = 0
      }
    }
  }
  
  vcv.tree = vcv(tree)
  vcv.sig2.product = kronecker(sig2,vcv.tree)
  det.product = det(vcv.sig2.product)
  data.root.difference = (sim_values - (D%*%matrix(root)))
  
  #exponent = exp(-1/(2*(t(data.root.difference)%*%inv(vcv.sig2.product)%*%data.root.difference)))
  exponent = exp(-0.5*(t(data.root.difference)%*%inv(vcv.sig2.product)%*%data.root.difference))

  denominator = sqrt(((2*pi)^(length(tree$tip.label)*length(root)))*det.product)
  
  
  #log.likelihood = log10(exponent / denominator)
  #log.likelihood = log(exponent / denominator)
  log.likelihood =  -0.5*(t(data.root.difference)%*%inv(vcv.sig2.product)%*%data.root.difference)-0.5*log(det.product)-length(tree$tip.label)*length(root)*log(2*pi)*0.5
  return(log.likelihood)
}


