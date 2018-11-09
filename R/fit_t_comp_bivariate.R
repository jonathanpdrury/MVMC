phylo = pbtree(n=20)
sig2.matrix = matrix(c(0.05,0.025,0.025,0.05),ncol=2)
s.term = matrix(c(0.01,0,0,0.01),ncol=2)
pars = list(sig2.matrix,s.term)
pars
root.values = c(0,0)
Nsegments = 1000
model = "MC"

sim_t_comp<-function(phylo,pars,root.values,Nsegments=1000,model="BM,OU,MC"){
  require(phytools)
  require(MASS)
  
  #return error if non-ultrametric tree
  if(phylo$Nnode!=(length(phylo$tip.label)-1)){stop("phylo object must be ultrametric")}
  
  ##check the parameters provided match up with the model given.
  if (class(pars)!="list"){stop("pars must be a list")}
  if (model=="BM"){
    if (length(pars)!=1||class(pars[[1]])!="matrix"){stop("Brownian Simuation parameters must only contain a sig2 matrix")}
  }
  if (model=="OU"){
    if (length(pars)!=3){stop("OU simulation must use 3 parameters for sig2, theta and mu respectively")}
    if (ncol(pars[[1]])!=ncol(pars[[2]])){stop("sig2 and theta must be matrices of equal size")}
    if (length(pars[[3]])!=ncol(pars[[2]])){stop("mu must be a vector with a length equal to the amount of traits")}
  }
  if (model=="MC"){
    if (length(pars)!=2){stop("MC simulation must use 2 parameters for sig2 and S respectively")}
    if (ncol(pars[[1]])!=ncol(pars[[2]])){stop("MC simulation parameters must be two equal size matrices")}
  }
  
  ##store parameters provided by user depending on model provided
  sig2=pars[[1]]
  if (model=="OU"){
    theta.matrix = pars[[2]]
    mu.vector = pars[[3]]
  }
  if (model=="MC"){
    s.term.matrix = pars[[2]]
  }
  
  ##check that root.value is a vector and assign to a vector
  if (length(root.values)!=ncol(sig2)){stop("There must be one root value for each trait")}
  ROOT<-root.values
  
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
  
  ##define a few useful variables
  ##nodeDist, the distance between each node and the root
  ##nodeDiff, the distance between each node and the previous node (times for integration)
  
  nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
  root <- length(phylo$tip.label) + 1
  heights<-nodeHeights(phylo)
  for (i in 1:dim(phylo$edge)[1]){
    nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
  }
  nodeDist<-c(nodeDist,max(heights))
  nodeDiff<-diff(nodeDist)
  ###label the branches for each segment of tree to be integrated and identify the node at which the branch terminates
  
  if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
    node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
    node.order<-node.order+length(phylo$tip.label)
    old.edge<-phylo$edge
    phylo$edge[,1]<-node.order
    for(j in 1:length(phylo$edge[,2])){
      if(phylo$edge[j,2]>length(phylo$tip.label)){
        #match number order in old edge
        #lookup value in new edge
        #replace with value
        phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
      }
    }
    nodeDist<-vector()
    for (i in 1:dim(phylo$edge)[1]){
      nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
    }
    nodeDist<-c(nodeDist,max(heights))
    nodeDiff<-diff(nodeDist)
  }
  
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:phylo$Nnode){
    other<-phylo$edge[phylo$edge[,1]==i+length(phylo$tip.label), 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+length(phylo$tip.label)
      if(b>length(phylo$tip.label)){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-phylo$tip.label[[b]]
        int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
      }
      mat<-rbind(mat,int)
    }
  }		
  
  
  ##now come up with a list of branches that exist for each time segment (+1 at each branching in this version, which means tree can't have any polytomies or speciation events at exactly the same time)
  nat<-list()
  for(i in 1:length(nodeDiff)){
    if(i==1){
      nat[[i]]<-list(mat[mat[,1]==(length(phylo$tip.label)+i),2])} else {
        IN<-vector()
        P<-mat[as.numeric(mat[,1])<=(length(phylo$tip.label)+i),c(2,3)]
        IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(length(phylo$tip.label)+i),1])
        nat[[i]]<-list(IN)
      }
  }	
  for(i in 2:length(nodeDiff)){			##THIS LOOP checks for an error
    if(length(unlist(nat[[i]]))!=(length(unlist(nat[[i-1]]))+1)){
      print(paste("ERROR at node",i+length(phylo$tip.label)))
    }	
  }
  
  ##which model is being simulated?
  if(is.na(match(model,c("BM","OU","MC")))){stop("model not specified correctly, must be 'BM', 'OU' or 'MC'")}
  
  ##define simulator for each time step of the model
  if (model == "BM") {
    simvalueBM = function(sig2.matrix, start.value.vector, seglen) {
      x = start.value.vector + t(t(mvrnorm(n=2 ,mu=0, Sigma=1))%*%chol(sig2.matrix*seglen))
    }
  }
  if (model == "OU") {
    simvalueOU = function(sig2.matrix, theta.matrix, mu.vector, start.value.vector, seglen) {
      x = start.value.vector + theta.matrix%*%(mu.vector-start.value.vector)*seglen + t(t(mvrnorm(n=2 ,mu=0, Sigma=1))%*%chol(sig2.matrix*seglen))
    }
  }
  if (model == "MC") {
    simvalueMC = function(sig2.matrix, s.term.matrix, mu.vector, start.value.vector, seglen) {
      x = start.value.vector + s.term.matrix%*%(mu.vector-start.value.vector)*seglen + t(t(mvrnorm(n=2 ,mu=0, Sigma=1))%*%chol(sig2.matrix*seglen))
    }
  }
  
  ##define the number of smaller segments to divide a tree into
  N<-Nsegments
  seglength<-tail(nodeDist,n=1)/N
  
  
  #these create the objects where results are stored--so far, these store univarate data; would need to be updated to store multivariate data
  masterbranch.1<-list()
  masterbranch.2<-list()
  masterseg<-list()
  
  for(i in 1:phylo$Nnode){ ##for each node interval
    if(i==1){
      mu1<-ROOT[1]
      mu2<-ROOT[2]
      } #initialize mu value at ROOT 
    
    branchespresent<-length(unlist(nat[[i]]))
    masterbranch.1[[i]]<-as.list(rep(ROOT[1],branchespresent))
    names(masterbranch.1[[i]])<-unlist(nat[[i]])
    masterbranch.2[[i]]<-as.list(rep(ROOT[2],branchespresent))
    names(masterbranch.2[[i]])<-unlist(nat[[i]])
    
    if(i>1){
      #update starting values to be the ending values from last time
      for(l in 1:length(masterbranch.1[[i]])){
        if(!is.na(match(names(masterbranch.1[[i]][l]),names(masterbranch.1[[i-1]])))){ #if the name of the element matches branch present in previous iteration
          masterbranch.1[[i]][l]<-sapply(masterbranch.1[[i-1]][match(names(masterbranch.1[[i]][l]),names(masterbranch.1[[i-1]]))],tail,n=1)
          masterbranch.2[[i]][l]<-sapply(masterbranch.2[[i-1]][match(names(masterbranch.2[[i]][l]),names(masterbranch.2[[i-1]]))],tail,n=1)
        } else {
          ##look up branch from which descended
          prev.branch<-mat[mat[,3]==mat[mat[,2]==names(masterbranch.1[[i]][l]),1],2]
          ##look up starting value
          masterbranch.1[[i]][l]<-sapply(masterbranch.1[[i-1]][match(prev.branch,names(masterbranch.1[[i-1]]))],tail,n=1)
          masterbranch.2[[i]][l]<-sapply(masterbranch.2[[i-1]][match(prev.branch,names(masterbranch.2[[i-1]]))],tail,n=1)
        }
      }
    }
    
    masterseg[[i]]<-as.list(rep(0,branchespresent))
    names(masterseg[[i]])<-unlist(nat[[i]])
    for(m in 1:length(masterseg[[i]])){
      if(i==1){
        if(seglength<nodeDiff[i]){
          masterseg[[i]][[m]]<-c(masterseg[[i]][[m]],seq(seglength,nodeDiff[i],seglength))} else{
            masterseg[[i]][[m]]<-c(masterseg[[i]][[m]],nodeDiff[i])	
          }} else{
            masterseg[[i]][[m]]<-seq(nodeDist[i],nodeDist[i+1],seglength)
          }
      if(tail(masterseg[[i]][[m]],n=1)!=nodeDiff[i]){
        masterseg[[i]][[m]]<-c(masterseg[[i]][[m]],nodeDist[i+1])
      }
    }
    
    for(k in 1:length(diff(masterseg[[i]][[1]]))){ ##for each of the segments between intervals
      for(j in 1:branchespresent){## for each of the branches present
        ##extract starting value
        start.value<-c(tail(masterbranch.1[[i]][[j]],n=1),tail(masterbranch.2[[i]][[j]],n=1)) #replace 1 with matching function
        #run simvalue function with mu (where to get value of mu and update?)
        segsize<-diff(masterseg[[i]][[1]])[k]
        mu.value<-c(tail(mu1,n=1),tail(mu2,n=1))
        if(model=="BM"){
          sv<-simvalueBM(sig2,start.value,segsize)
        }
        if(model=="OU"){
          sv<-simvalueOU(sig2,theta.matrix,mu.vector,start.value,segsize)
        }
        if(model=="MC"){
          sv<-simvalueMC(sig2,s.term.matrix,mu.value,start.value,segsize)}
        masterbranch.1[[i]][[j]]<-c(masterbranch.1[[i]][[j]],sv[1])
        masterbranch.2[[i]][[j]]<-c(masterbranch.2[[i]][[j]],sv[2])
      }
      ##THINK about whether the placement of mu was wrong in original simulations
      mu1=c(mu1,mean(sapply(masterbranch.1[[i]],tail,n=1)))# the last element in each element of branchlist/branches present ##this assumes that the first element is updated for each i
      mu2=c(mu2,mean(sapply(masterbranch.2[[i]],tail,n=1)))# the last element in each element of branchlist/branches present ##this assumes that the first element is updated for each i
    }
  }
  final.values.1 = lapply(c(masterbranch.1[[phylo$Nnode]]),tail,n=1)
  final.values.2 = lapply(c(masterbranch.2[[phylo$Nnode]]),tail,n=1)
  final.values = list()
  for (i in 1:length(phylo$tip.label)){
    final.values[[i]] = c(final.values.1[[i]],final.values.2[[i]])
    names(final.values)[i] = names(final.values.1[i])
  }
  
  return(final.values)
}




