sim_data_matrix = matrix(nrow=length(tree$tip.label),ncol=2)
for (j in 1:length(tree$tip.label)){
  sim_data_matrix[j,1] = sim[[j]][1]
  sim_data_matrix[j,2] = sim[[j]][2]
}
colnames(sim_data_matrix) = c("trait.1","trait.2")
rownames(sim_data_matrix) = names(sim)

fit_estimate = mvBM(tree,sim_data_matrix[tree$tip.label,],model="BM1")