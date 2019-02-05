library(phytools)
Nsim = 1:100
root = c(0,0)
tree.sizes = c(20,50,100,150)
tree.list = list()
for (i in 1:length(tree.sizes)){
  tree.list[[i]] = list()
  for (j in Nsim){
    tree = pbtree(n=tree.sizes[i])
    tree$edge.length = 10/max(nodeHeights(tree))*tree$edge.length
    while(min(diff(sort(branching.times(tree)))) >= 10/1000000){##length of the tree/Nsim in sim_t_comp_bivariate
      tree = NULL
      tree = (pbtree(n=tree.sizes[i]))
    }
    tree.list[[i]][[j]] = tree
  }
}

matrix1 = matrix(c(1,0,0,0.5),ncol=2)
matrix2 = matrix(c(1,0.5,0.5,0.5),ncol=2)
matrix3 = matrix(c(1,-0.5,-0.5,0.5),ncol=2)
sig2.matrices = list(matrix1,matrix2,matrix3)

##when DD with a positive secondary term
DDpos.matrix1 = matrix(c(0.01,0,0,0.005),ncol=2)
DDpos.matrix2 = matrix(c(0.01,0.005,0.005,0.005),ncol=2)
DDpos.matrix3 = matrix(c(0.01,-0.005,-0.005,0.005),ncol=2)
DDpos.sig2.matrices = list(DDPos.matrix1,DDpos.matrix2,DDpos.matrix3)

pars1 = matrix(c("pars.1",0,0,"pars.2"),ncol=2)
pars2 = matrix(c("pars.1","cov.pars","cov.pars","pars.2"),ncol=2)
pars3 = matrix(c("pars.1","-cov.pars","-cov.pars","pars.2"),ncol=2)
pars4 = matrix(c("pars.1","cov.pars",0,"pars.2"),ncol=2)
pars.format = list(pars1,pars2,pars3,pars4)
half.lives = c(0,0.5,1,5)
OU.theta = c(1,2)

##DD with negative secondary terms
matrix11 = matrix(c(0.5,0,0,0.25),ncol=2)
matrix12 = matrix(c(0.1,0,0,0.05),ncol=2)
matrix13 = matrix(c(0.01,0,0,0.005),ncol=2)
matrix21 = matrix(c(0.5,0.25,0.25,0.25),ncol=2)
matrix22 = matrix(c(0.1,0.05,0.05,0.05),ncol=2)
matrix23 = matrix(c(0.01,0.005,0.005,0.005),ncol=2)
matrix31 = matrix(c(0.5,-0.25,-0.25,0.25),ncol=2)
matrix32 = matrix(c(0.1,-0.05,-0.05,0.05),ncol=2)
matrix33 = matrix(c(0.01,-0.005,-0.005,0.005),ncol=2)
DD.tip.values = list(list(matrix11,matrix12,matrix13),list(matrix21,matrix22,matrix23),list(matrix31,matrix32,matrix33))

##DD with positive secondary terms
DDpos.matrix11 = matrix(c(0.1,0,0,0.05),ncol=2)
DDpos.matrix12 = matrix(c(0.5,0,0,0.25),ncol=2)
DDpos.matrix13 = matrix(c(1,0,0,0.5),ncol=2)
DDpos.matrix21 = matrix(c(0.1,0.05,0.05,0.05),ncol=2)
DDpos.matrix22 = matrix(c(0.5,0.25,0.25,0.25),ncol=2)
DDpos.matrix23 = matrix(c(1,0.5,0.5,0.5),ncol=2)
DDpos.matrix31 = matrix(c(0.1,-0.05,-0.05,0.05),ncol=2)
DDpos.matrix32 = matrix(c(0.5,-0.25,-0.25,0.25),ncol=2)
DDpos.matrix33 = matrix(c(1,-0.5,-0.5,0.5),ncol=2)
DD.root.values = list(list(DDpos.matrix11,DDpos.matrix12,DDpos.matrix13),list(DDpos.matrix21,DDpos.matrix22,DDpos.matrix23),list(DDpos.matrix31,DDpos.matrix32,DDpos.matrix33))
