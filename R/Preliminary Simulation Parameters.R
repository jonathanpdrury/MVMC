library(phytools)
Nsim = 20
root = c(0,0)
tree.sizes = c(50,100)
tree.list = list()
for (i in 1:length(tree.sizes)){
  tree.list[[i]] = list()
  for (j in 1:Nsim){
    tree = pbtree(n=tree.sizes[i])
    tree$edge.length = 10/max(nodeHeights(tree))*tree$edge.length
    tree.list[[i]][[j]] = tree
  }
}

matrix1 = matrix(c(1,0,0,0.5),ncol=2)
matrix2 = matrix(c(1,0.5,0.5,0.5),ncol=2)
sig2.matrices = list(matrix1,matrix2)

##when DD with a positive secondary term
DDpos.matrix1 = matrix(c(0.01,0,0,0.005),ncol=2)
DDpos.matrix2 = matrix(c(0.01,0.005,0.005,0.005),ncol=2)
DDpos.sig2.matrices = list(matrix1,matrix2)

pars1 = matrix(c("pars.1",0,0,"pars.2"))
pars2 = matrix(c("pars.1","cov.pars","cov.pars","pars.2"))
pars.list = list(pars1,pars2)
half.lives = c(0.5,1)
OU.theta = c(1,2)

##DD with negative secondary terms
matrix11 = matrix(c(0.5,0,0,0.25),ncol=2)
matrix12 = matrix(c(0.1,0,0,0.05),ncol=2)
matrix21 = matrix(c(0.5,0.25,0.25,0.25),ncol=2)
matrix22 = matrix(c(0.1,0.05,0.05,0.05),ncol=2)
DD.tip.values = list(list(matrix11,matrix12),list(matrix21,matrix22))

##DD with positive secondary terms
DDpos.matrix12 = matrix(c(0.5,0,0,0.25),ncol=2)
DDpos.matrix13 = matrix(c(1,0,0,0.5),ncol=2)
DDpos.matrix22 = matrix(c(0.5,0.25,0.25,0.25),ncol=2)
DDpos.matrix23 = matrix(c(1,0.5,0.5,0.5),ncol=2)
DD.root.values = list(list(DDpos.matrix12,DDpos.matrix13),list(DDpos.matrix22,DDpos.matrix23))