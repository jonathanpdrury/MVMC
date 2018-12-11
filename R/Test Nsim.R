##gives probability with max 1, the number of trees which have a minimum branching time greater than the accepted difference
##accepted difference calculated by tree.length/Nsim (as in sim_t_comp_bivariate)
test_min_branching_point = function(
  tree.size,
  tree.length,
  Nsim,
  Ntrees
){
  require(phytools)
  successful.trees = 0
  for (i in 1:Ntrees){
    tree = pbtree(n=tree.size)
    tree$edge.length = tree.length/max(nodeHeights(tree))*tree$edge.length
    if (min(diff(sort(branching.times(tree)))) >= tree.length/Nsim){
      successful.trees = successful.trees + 1
    }
  }
  chance.success = successful.trees/Ntrees
  return(chance.success)
}
