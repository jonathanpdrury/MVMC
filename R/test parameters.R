library(phytools)
matrix1 = matrix(c(1,0,0,0.5),ncol=2)
matrix2 = matrix(c(1,0.5,0.5,1),ncol=2)
sig2.matrices = list(matrix1,matrix2)
root = c(0,0)
pars1 = matrix(c("pars.1",0,0,"pars.2"),ncol=2)
pars2 = matrix(c("pars.1","cov.pars","cov.pars","pars.2"),ncol=2)
pars.format = list(pars1,pars2)
half.lives = c(1,5)
OU.theta = c(1,2)
tree.sizes = c(5,10)
Nsim=2
tree.list = list()
for (i in 1:length(tree.sizes)){
  tree.list[[i]] = list()
  for (j in 1:Nsim){
    tree = pbtree(n=tree.sizes[i])
    tree$edge.length = 10/max(nodeHeights(tree))*tree$edge.length
    tree.list[[i]][[j]] = tree
  }
}
names(tree.list) = tree.sizes
DD.matrix11 = matrix(c(0.5,0,0,0.25),ncol=2)
DD.matrix12 = matrix(c(0.25,0,0,0.125),ncol=2)
DD.matrix21 = matrix(c(0.5,0.25,0.25,0.5),ncol=2)
DD.matrix22 = matrix(c(0.25,0.25,0.25,0.5),ncol=2)
DD.tip.rate = list(list(DD.matrix11,DD.matrix12),list(DD.matrix21,DD.matrix22))
