require(phytools)
require(matlib)
require(msos)

setwd('J:/My_Documents/Project/Simulation Study/Preliminary Simulation/complete_sim_results/Likelihood Test')
source('PhenotypicModel.R') #from RPANDA (Marc Manceau)
source('PhenotypicADiag.R') #from RPANDA (Marc Manceau)
source('DDexp_nogeo_ADiag.R') #from RPANDA 
source('DDlin_nogeo_ADiag.R') #from RPANDA
source('Likelihood_Calculator.R')

max_likelihood_calculator_BM_DD = function(
  model,
  tree.size,
  sig2.matrix
){
  if (model == "BM"){
    setwd('J:/My_Documents/Project/Simulation Study/Preliminary Simulation/complete_sim_results/complete_BM')
    load('BM_sig2_values.RData')
    load('tree_list.RData')
    
    ## load data and save as sim.results to avoid excessive eval parse paste
    x = eval(
      parse(
        text=paste(
          "load('BM_sim_tree_",
          tree.size,
          "_sig2_",
          sig2.matrix,
          ".RData')",
          sep = ""
        )
      )
    )
    sim.results = get(x)
    
    sig2.1.estimate = c()
    sig2.2.estimate = c()
    sig2.3.estimate = c()
    sig2.4.estimate = c()
    root.1.estimate = c()
    root.2.estimate = c()
    
    for (i in length(sim.results)){
      result = optim(
        par = c(sig2.matrices[[sig2.matrix]],c(0,0)),
        fn = log_likelihood_mv_BM_DD,
        tree = tree.list[[toString(tree.size)]][[i]],
        sim.value = sim.results[[i]],
        model = model,
        optim = TRUE
      )
      
      sig2.1.estimate = c(sig2.1.estimate, result[[1]][1])
      sig2.2.estimate = c(sig2.2.estimate, result[[1]][2])
      sig2.3.estimate = c(sig2.3.estimate, result[[1]][3])
      sig2.4.estimate = c(sig2.4.estimate, result[[1]][4])
      root.1.estimate = c(root.1.estimate, result[[1]][5])
      root.2.estimate = c(root.2.estimate, result[[1]][6])
    }
    
    sig2.1.average = mean(sig2.1.estimate)
    sig2.2.average = mean(sig2.2.estimate)
    sig2.3.average = mean(sig2.3.estimate)
    sig2.4.average = mean(sig2.4.estimate)
    root.1.average = mean(root.1.estimate)
    root.2.average = mean(root.2.estimate)
    
    final.estimate = list(
      matrix(
        c(
          sig2.1.average,
          sig2.2.average,
          sig2.3.average,
          sig2.4.average 
        ),
        ncol = 2
      ),
      c(
        root.1.average,
        root.2.average
      )
    )
    names(final.estimate) = c("sig2","root")
    return(final.estimate)
  }
}









































