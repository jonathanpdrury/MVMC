.libPaths("/ddn/data/ckkr89/R/x86_64-pc-linux-gnu-library/3.4")
setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/Likelihood Test')
#setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/Likelihood Test')
source('PhenotypicModel.R') #from RPANDA (Marc Manceau)
source('PhenotypicADiag.R') #from RPANDA (Marc Manceau)
source('DDexp_nogeo_ADiag.R') #from RPANDA 
source('DDlin_nogeo_ADiag.R') #from RPANDA
source('mv_BM_DD_Likelihood_Calculator.R')
require(phytools)
require(mvMORPH)
require(MASS)
require(msos)
require(parallel)
require(snow)


parallel_max_likelihood_calculator_BM_OU_DD = function(
  model,
  Nsim
){
  require(phytools)
  require(mvMORPH)
  require(matlib)
  require(msos)
  require(parallel)
  require(snow)
  
  if (model == "BM"){
    setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/complete_BM')
    #setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/complete_BM')
    load('BM_sig2_values.RData')
    load('tree_list.RData')
    
    if (!dir.exists('likelihood_results')){
      dir.create('likelihood_results')
    }
    
    rand_string_create = function(n=1,length=12){
      randomString = c(1:n) # initialize vector
      for (i in 1:n){
        randomString[i] = paste(
          sample(
            c(0:9,letters,LETTERS),
            length,
            replace=TRUE
          ),
          collapse=""
        )
      }
      return(randomString)
    }
    rand.string = rand_string_create()
    
    eval(
      parse(
        text=paste(
          "dir.create('likelihood_results/BM_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    
    tree.size = c()
    tree.number = c()
    sig2.1 = c()
    sig2.2 = c()
    sig2.3 = c()
    sig2.4 = c()
    root.1 = c()
    root.2 = c()
    est.sig2.1 = c()
    est.sig2.2 = c()
    est.sig2.3 = c()
    est.sig2.4 = c()
    est.root.1 = c()
    est.root.2 = c()
    log.likelihood = c()
    convergence = c()
    
    for (i in 1:length(tree.list)){
      for (j in 1:length(sig2.matrices)){
        
        ## load data and save as sim.results to avoid excessive eval parse paste
        x = eval(
          parse(
            text=paste(
              "load('BM_sim_tree_",
              names(tree.list[i]),
              "_sig2_",
              j,
              ".RData')",
              sep = ""
            )
          )
        )
        sim.results = get(x)
        
        for (k in Nsim){
          sim.results.matrix = matrix(nrow=length(tree.list[[i]][[k]]$tip.label),ncol=2)
          
          for (l in 1:length(tree.list[[i]][[k]]$tip.label)){
            sim.results.matrix[l,1] = sim.results[[k]][[l]][1]
            sim.results.matrix[l,2] = sim.results[[k]][[l]][2]
          }
          
          colnames(sim.results.matrix) = c("trait.1","trait.2")
          rownames(sim.results.matrix) = names(sim.results[[k]])
          
          result = mvBM(
            tree = tree.list[[i]][[k]],
            data = sim.results.matrix[tree.list[[i]][[k]]$tip.label,],
            model = "BM1"
          )
          
          #add results to vectors for the data frame
          tree.size = c(tree.size,names(tree.list[i]))
          tree.number = c(tree.number,k)
          sig2.1 = c(sig2.1,sig2.matrices[[j]][1])
          sig2.2 = c(sig2.2,sig2.matrices[[j]][2])
          sig2.3 = c(sig2.3,sig2.matrices[[j]][3])
          sig2.4 = c(sig2.4,sig2.matrices[[j]][4])
          root.1 = c(root.1,0)
          root.2 = c(root.2,0)
          est.sig2.1 = c(est.sig2.1,result[[5]][1])
          est.sig2.2 = c(est.sig2.2,result[[5]][2])
          est.sig2.3 = c(est.sig2.3,result[[5]][3])
          est.sig2.4 = c(est.sig2.4,result[[5]][4])
          est.root.1 = c(est.root.1,result[[4]][1])
          est.root.2 = c(est.root.2,result[[4]][2])
          log.likelihood = c(log.likelihood,result[[1]])
          convergence = c(convergence,result[[6]])
        }
      }
    }
    
    BM_data = data.frame(
      tree.size,
      tree.number,
      sig2.1,
      sig2.2,
      sig2.3,
      sig2.4,
      root.1,
      root.2,
      est.sig2.1,
      est.sig2.2,
      est.sig2.3,
      est.sig2.4,
      est.root.1,
      est.root.2,
      log.likelihood,
      convergence
    )
    
    eval(
      parse(
        text=paste(
          "setwd('likelihood_results/BM_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    write.csv(BM_data, file = "BM_Data.csv")
    
  } else if (model == "OU"){
    
    setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/complete_OU')
    #setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/complete_OU')
    load('OU_sig2_values.RData')
    load('OU_alpha_matrices.RData')
    load('OU_theta_values.RData')
    load('tree_list.RData')
    
    if (!dir.exists('likelihood_results')){
      dir.create('likelihood_results')
    }
    
    rand_string_create = function(n=1,length=12){
      randomString = c(1:n) # initialize vector
      for (i in 1:n){
        randomString[i] = paste(
          sample(
            c(0:9,letters,LETTERS),
            length,
            replace=TRUE
          ),
          collapse=""
        )
      }
      return(randomString)
    }
    rand.string = rand_string_create()
    
    eval(
      parse(
        text=paste(
          "dir.create('likelihood_results/OU_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    
    tree.size = c()
    tree.number = c()
    sig2.1 = c()
    sig2.2 = c()
    sig2.3 = c()
    sig2.4 = c()
    theta.1 = c()
    theta.2 = c()
    alpha.1 = c()
    alpha.2 = c()
    alpha.3 = c()
    alpha.4 = c()
    est.sig2.1 = c()
    est.sig2.2 = c()
    est.sig2.3 = c()
    est.sig2.4 = c()
    est.theta.1 = c()
    est.theta.2 = c()
    est.alpha.1 = c()
    est.alpha.2 = c()
    est.alpha.3 = c()
    est.alpha.4 = c()
    log.likelihood = c()
    convergence = c()
    
    for (i in 1:length(tree.list)){
      for (j in 1:length(sig2.matrices)){
        for (k in 1:length(pars.list)){
          for (l in 1:length(pars.list[[k]])){
            x = eval(
              parse(
                text=paste(
                  "load('OU_sim_tree_",
                  names(tree.list)[i],
                  "_sig2_",
                  j,
                  "_half_life_",
                  names(pars.list)[k],
                  "_pars_",
                  l,
                  ".RData')",
                  sep = ""
                )
              )
            )
            sim.results = get(x)
            
            for (m in Nsim){
              sim.results.matrix = matrix(nrow=length(tree.list[[i]][[m]]$tip.label),ncol=2)
              
              for (n in 1:length(tree.list[[i]][[m]]$tip.label)){
                sim.results.matrix[n,1] = sim.results[[m]][[n]][1]
                sim.results.matrix[n,2] = sim.results[[m]][[n]][2]
              }
              
              colnames(sim.results.matrix) = c("trait.1","trait.2")
              rownames(sim.results.matrix) = names(sim.results[[m]])
              
              result = mvOU(
                tree = tree.list[[i]][[m]],
                data = sim.results.matrix[tree.list[[i]][[k]]$tip.label,],
                model = "OU1"
              )
              
              tree.size = c(tree.size,names(tree.list[i]))
              tree.number = c(tree.number,k)
              sig2.1 = c(sig2.1,sig2.matrices[[j]][1])
              sig2.2 = c(sig2.2,sig2.matrices[[j]][2])
              sig2.3 = c(sig2.3,sig2.matrices[[j]][3])
              sig2.4 = c(sig2.4,sig2.matrices[[j]][4])
              theta.1 = c(theta.1,OU.theta)
              theta.2 = c(theta.2,OU.theta)
              alpha.1 = c(alpha.1,pars.list[[k]][[l]][1])
              alpha.2 = c(alpha.2,pars.list[[k]][[l]][2])
              alpha.3 = c(alpha.3,pars.list[[k]][[l]][3])
              alpha.4 = c(alpha.4,pars.list[[k]][[l]][4])
              est.sig2.1 = c(est.sig2.1,result[[6]][1])
              est.sig2.2 = c(est.sig2.2,result[[6]][2])
              est.sig2.3 = c(est.sig2.3,result[[6]][3])
              est.sig2.4 = c(est.sig2.4,result[[6]][4])
              est.theta.1 = c(est.theta.1,result[[4]][1])
              est.theta.2 = c(est.theta.2,result[[4]][2])
              est.alpha.1 = c(est.alpha.1,result[[5]][1])
              est.alpha.2 = c(est.alpha.1,result[[5]][2])
              est.alpha.3 = c(est.alpha.1,result[[5]][3])
              est.alpha.4 = c(est.alpha.1,result[[5]][4])
              log.likelihood = c(log.likelihood,result[[1]])
              convergence = c(convergence,result[[7]])
            }
          }
        }
      }
    }
    
    OU_data = data.frame(
      tree.size,
      tree.number,
      sig2.1,
      sig2.2,
      sig2.3,
      sig2.4,
      theta.1,
      theta.2,
      alpha.1,
      alpha.2,
      alpha.3,
      alpha.4,
      est.sig2.1,
      est.sig2.2,
      est.sig2.3,
      est.sig2.4,
      est.theta.1,
      est.theta.2,
      est.alpha.1,
      est.alpha.2,
      est.alpha.3,
      est.alpha.4,
      log.likelihood,
      convergence
    )
    
    eval(
      parse(
        text=paste(
          "setwd('likelihood_results/OU_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    write.csv(OU_data, file = "OU_data.csv")
    
  } else if (model == "MC"){
    ##TBD
  } else if (model == "DDexp_neg"){
    setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/complete_DDexp_neg')
    #setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/complete_DDexp_neg')
    load('DDexp_root_sig2_values.RData')
    load('DDexp_r_term_matrices.RData')
    load('tree_list.RData')
    
    if (!dir.exists('likelihood_results')){
      dir.create('likelihood_results')
    }
    
    rand_string_create = function(n=1,length=12){
      randomString = c(1:n) # initialize vector
      for (i in 1:n){
        randomString[i] = paste(
          sample(
            c(0:9,letters,LETTERS),
            length,
            replace=TRUE
          ),
          collapse=""
        )
      }
      return(randomString)
    }
    rand.string = rand_string_create()
    
    eval(
      parse(
        text=paste(
          "dir.create('likelihood_results/DDexp_neg_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    
    tree.size = c()
    tree.number = c()
    sig2.1 = c()
    sig2.2 = c()
    sig2.3 = c()
    sig2.4 = c()
    root.1 = c()
    root.2 = c()
    r.term.1 = c()
    r.term.2 = c()
    r.term.3 = c()
    r.term.4 = c()
    est.sig2.1 = c()
    est.sig2.2 = c()
    est.sig2.3 = c()
    est.sig2.4 = c()
    est.root.1 = c()
    est.root.2 = c()
    est.r.term.1 = c()
    est.r.term.2 = c()
    est.r.term.3 = c()
    est.r.term.4 = c()
    log.likelihood = c()
    convergence = c()
    
    for (i in 1:length(tree.list)){
      for (j in 1:length(sig2.matrices)){
        for (k in 1:length(pars.list[[i]][[j]])){
          x = eval(
            parse(
              text=paste(
                "load('DDexp_neg_sim_tree_",
                names(tree.list)[i],
                "_root_sig2_",
                j,
                "_tip_sig2_",
                k,
                ".RData')",
                sep = ""
              )
            )
          )
          sim.results = get(x)
          
          for (l in Nsim){
            result = optim(
              par = c(0.5,0,0,1,1,0.5,-0.1,-0.05,-0.05,-0.1),
              fn = log_likelihood_mv_BM_DD,
              tree = tree.list[[i]][[l]],
              sim.value = sim.results[[l]],
              model = "DDexp",
              optim = TRUE
            )
            
            tree.size = c(tree.size,names(tree.list[i]))
            tree.number = c(tree.number,k)
            sig2.1 = c(sig2.1,sig2.matrices[[j]][1])
            sig2.2 = c(sig2.2,sig2.matrices[[j]][2])
            sig2.3 = c(sig2.3,sig2.matrices[[j]][3])
            sig2.4 = c(sig2.4,sig2.matrices[[j]][4])
            root.1 = c(root.1,0)
            root.2 = c(root.2,0)
            r.term.1 = c(r.term.1,pars.list[[i]][[j]][[k]][1])
            r.term.2 = c(r.term.2,pars.list[[i]][[j]][[k]][2])
            r.term.3 = c(r.term.3,pars.list[[i]][[j]][[k]][3])
            r.term.4 = c(r.term.4,pars.list[[i]][[j]][[k]][4])
            est.sig2.1 = c(est.sig2.1,result[[1]][1])
            est.sig2.2 = c(est.sig2.2,result[[1]][2])
            est.sig2.3 = c(est.sig2.3,result[[1]][3])
            est.sig2.4 = c(est.sig2.4,result[[1]][4])
            est.root.1 = c(est.root.1,result[[1]][5])
            est.root.2 = c(est.root.2,result[[1]][6])
            est.r.term.1 = c(est.r.term.1,result[[1]][7])
            est.r.term.2 = c(est.r.term.2,result[[1]][8])
            est.r.term.3 = c(est.r.term.3,result[[1]][9])
            est.r.term.4 = c(est.r.term.4,result[[1]][10])
            log.likelihood = c(log.likelihood,result[[2]])
            convergence = c(convergence,result[[4]])
          }
        }
      }
    }
    
    DDexp_neg_data = data.frame(
      tree.size,
      tree.number,
      sig2.1,
      sig2.2,
      sig2.3,
      sig2.4,
      root.1,
      root.2,
      r.term.1,
      r.term.2,
      r.term.3,
      r.term.4,
      est.sig2.1,
      est.sig2.2,
      est.sig2.3,
      est.sig2.4,
      est.root.1,
      est.root.2,
      est.r.term.1,
      est.r.term.2,
      est.r.term.3,
      est.r.term.4,
      log.likelihood,
      convergence
    )
    
    eval(
      parse(
        text=paste(
          "setwd('likelihood_results/DDexp_neg_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    write.csv(DDexp_neg_data, file = "DDexp_neg_data.csv")
    
  } else if (model == "DDexp_pos"){
    setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/complete_DDexp_pos')
    #setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/complete_DDexp_pos')
    load('DDexp_root_sig2_values.RData')
    load('DDexp_r_term_matrices.RData')
    load('tree_list.RData')
    
    if (!dir.exists('likelihood_results')){
      dir.create('likelihood_results')
    }
    
    rand_string_create = function(n=1,length=12){
      randomString = c(1:n) # initialize vector
      for (i in 1:n){
        randomString[i] = paste(
          sample(
            c(0:9,letters,LETTERS),
            length,
            replace=TRUE
          ),
          collapse=""
        )
      }
      return(randomString)
    }
    rand.string = rand_string_create()
    
    eval(
      parse(
        text=paste(
          "dir.create('likelihood_results/DDexp_pos_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    
    tree.size = c()
    tree.number = c()
    sig2.1 = c()
    sig2.2 = c()
    sig2.3 = c()
    sig2.4 = c()
    root.1 = c()
    root.2 = c()
    r.term.1 = c()
    r.term.2 = c()
    r.term.3 = c()
    r.term.4 = c()
    est.sig2.1 = c()
    est.sig2.2 = c()
    est.sig2.3 = c()
    est.sig2.4 = c()
    est.root.1 = c()
    est.root.2 = c()
    est.r.term.1 = c()
    est.r.term.2 = c()
    est.r.term.3 = c()
    est.r.term.4 = c()
    log.likelihood = c()
    convergence = c()
    
    for (i in 1:length(tree.list)){
      for (j in 1:length(sig2.matrices)){
        for (k in 1:length(pars.list[[i]][[j]])){
          x = eval(
            parse(
              text=paste(
                "load('DDexp_pos_sim_tree_",
                names(tree.list)[i],
                "_root_sig2_",
                j,
                "_tip_sig2_",
                k,
                ".RData')",
                sep = ""
              )
            )
          )
          sim.results = get(x)
          
          for (l in Nsim){
            result = optim(
              par = c(0.5,0,0,1,1,0.5,-0.1,-0.05,-0.05,-0.1),
              fn = log_likelihood_mv_BM_DD,
              tree = tree.list[[i]][[l]],
              sim.value = sim.results[[l]],
              model = "DDexp",
              optim = TRUE
            )
            
            tree.size = c(tree.size,names(tree.list[i]))
            tree.number = c(tree.number,k)
            sig2.1 = c(sig2.1,DDpos.sig2.matrices[[j]][1])
            sig2.2 = c(sig2.2,DDpos.sig2.matrices[[j]][2])
            sig2.3 = c(sig2.3,DDpos.sig2.matrices[[j]][3])
            sig2.4 = c(sig2.4,DDpos.sig2.matrices[[j]][4])
            root.1 = c(root.1,0)
            root.2 = c(root.2,0)
            r.term.1 = c(r.term.1,pars.list[[i]][[j]][[k]][1])
            r.term.2 = c(r.term.2,pars.list[[i]][[j]][[k]][2])
            r.term.3 = c(r.term.3,pars.list[[i]][[j]][[k]][3])
            r.term.4 = c(r.term.4,pars.list[[i]][[j]][[k]][4])
            est.sig2.1 = c(est.sig2.1,result[[1]][1])
            est.sig2.2 = c(est.sig2.2,result[[1]][2])
            est.sig2.3 = c(est.sig2.3,result[[1]][3])
            est.sig2.4 = c(est.sig2.4,result[[1]][4])
            est.root.1 = c(est.root.1,result[[1]][5])
            est.root.2 = c(est.root.2,result[[1]][6])
            est.r.term.1 = c(est.r.term.1,result[[1]][7])
            est.r.term.2 = c(est.r.term.2,result[[1]][8])
            est.r.term.3 = c(est.r.term.3,result[[1]][9])
            est.r.term.4 = c(est.r.term.4,result[[1]][10])
            log.likelihood = c(log.likelihood,result[[2]])
            convergence = c(convergence,result[[4]])
          }
        }
      }
    }
    
    DDexp_pos_data = data.frame(
      tree.size,
      tree.number,
      sig2.1,
      sig2.2,
      sig2.3,
      sig2.4,
      root.1,
      root.2,
      r.term.1,
      r.term.2,
      r.term.3,
      r.term.4,
      est.sig2.1,
      est.sig2.2,
      est.sig2.3,
      est.sig2.4,
      est.root.1,
      est.root.2,
      est.r.term.1,
      est.r.term.2,
      est.r.term.3,
      est.r.term.4,
      log.likelihood,
      convergence
    )
    
    eval(
      parse(
        text=paste(
          "setwd('likelihood_results/DDexp_pos_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    write.csv(DDexp_pos_data, file = "DDexp_pos_data.csv")
    
  } else if (model == "DDlin_neg"){
    setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/complete_DDlin_neg')
    #setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/complete_DDlin_neg')
    load('DDlin_root_sig2_values.RData')
    load('DDlin_slope_term_matrices.RData')
    load('tree_list.RData')
    
    if (!dir.exists('likelihood_results')){
      dir.create('likelihood_results')
    }
    
    rand_string_create = function(n=1,length=12){
      randomString = c(1:n) # initialize vector
      for (i in 1:n){
        randomString[i] = paste(
          sample(
            c(0:9,letters,LETTERS),
            length,
            replace=TRUE
          ),
          collapse=""
        )
      }
      return(randomString)
    }
    rand.string = rand_string_create()
    
    eval(
      parse(
        text=paste(
          "dir.create('likelihood_results/DDlin_neg_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    
    tree.size = c()
    tree.number = c()
    sig2.1 = c()
    sig2.2 = c()
    sig2.3 = c()
    sig2.4 = c()
    root.1 = c()
    root.2 = c()
    slope.term.1 = c()
    slope.term.2 = c()
    slope.term.3 = c()
    slope.term.4 = c()
    est.sig2.1 = c()
    est.sig2.2 = c()
    est.sig2.3 = c()
    est.sig2.4 = c()
    est.root.1 = c()
    est.root.2 = c()
    est.slope.term.1 = c()
    est.slope.term.2 = c()
    est.slope.term.3 = c()
    est.slope.term.4 = c()
    log.likelihood = c()
    convergence = c()
    
    for (i in Nsim){
      for (j in 1:length(sig2.matrices)){
        for (k in 1:length(pars.list[[i]][[j]])){
          x = eval(
            parse(
              text=paste(
                "load('DDlin_neg_sim_tree_",
                names(tree.list)[i],
                "_root_sig2_",
                j,
                "_tip_sig2_",
                k,
                ".RData')",
                sep = ""
              )
            )
          )
          sim.results = get(x)
          
          for (l in 1:length(sim.results)){
            result = optim(
              par = c(0.5,0,0,1,1,0.5,-0.1,-0.05,-0.05,-0.1),
              fn = log_likelihood_mv_BM_DD,
              tree = tree.list[[i]][[l]],
              sim.value = sim.results[[l]],
              model = "DDlin",
              optim = TRUE
            )
            
            #calculate the position of the values in the data frame
            data.frame.position = (i-1)*(j-1)*(k-1)*l + (j-1)*(k-1)*l + (k-1)*l + l
            
            tree.size = c(tree.size,names(tree.list[i]))
            tree.number = c(tree.number,k)
            sig2.1 = c(sig2.1,sig2.matrices[[j]][1])
            sig2.2 = c(sig2.2,sig2.matrices[[j]][2])
            sig2.3 = c(sig2.3,sig2.matrices[[j]][3])
            sig2.4 = c(sig2.4,sig2.matrices[[j]][4])
            root.1 = c(root.1,0)
            root.2 = c(root.2,0)
            r.term.1 = c(r.term.1,pars.list[[i]][[j]][[k]][1])
            r.term.2 = c(r.term.2,pars.list[[i]][[j]][[k]][2])
            r.term.3 = c(r.term.3,pars.list[[i]][[j]][[k]][3])
            r.term.4 = c(r.term.4,pars.list[[i]][[j]][[k]][4])
            est.sig2.1 = c(est.sig2.1,result[[1]][1])
            est.sig2.2 = c(est.sig2.2,result[[1]][2])
            est.sig2.3 = c(est.sig2.3,result[[1]][3])
            est.sig2.4 = c(est.sig2.4,result[[1]][4])
            est.root.1 = c(est.root.1,result[[1]][5])
            est.root.2 = c(est.root.2,result[[1]][6])
            est.r.term.1 = c(est.r.term.1,result[[1]][7])
            est.r.term.2 = c(est.r.term.2,result[[1]][8])
            est.r.term.3 = c(est.r.term.3,result[[1]][9])
            est.r.term.4 = c(est.r.term.4,result[[1]][10])
            log.likelihood = c(log.likelihood,result[[2]])
            convergence = c(convergence,result[[4]])
          }
        }
      }
    }
    
    DDlin_neg_data = data.frame(
      tree.size,
      tree.number,
      sig2.1,
      sig2.2,
      sig2.3,
      sig2.4,
      root.1,
      root.2,
      slope.term.1,
      slope.term.2,
      slope.term.3,
      slope.term.4,
      est.sig2.1,
      est.sig2.2,
      est.sig2.3,
      est.sig2.4,
      est.root.1,
      est.root.2,
      est.slope.term.1,
      est.slope.term.2,
      est.slope.term.3,
      est.slope.term.4,
      log.likelihood,
      convergence
    )
    
    eval(
      parse(
        text=paste(
          "setwd('likelihood_results/DDlin_neg_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    write.csv(DDlin_neg_data, file = "DDlin_neg_data.csv")
    
  } else if (model == "DDlin_pos"){
    setwd('/ddn/home/ckkr89/Project Data/complete_sim_results/complete_DDlin_pos')
    #setwd('//Hudson/ckkr89/My_Documents/Project/Simulation Study/Final Simulation/complete_sim_results/complete_DDlin_pos')
    load('DDlin_root_sig2_values.RData')
    load('DDlin_slope_term_matrices.RData')
    load('tree_list.RData')
    
    if (!dir.exists('likelihood_results')){
      dir.create('likelihood_results')
    }
    
    rand_string_create = function(n=1,length=12){
      randomString = c(1:n) # initialize vector
      for (i in 1:n){
        randomString[i] = paste(
          sample(
            c(0:9,letters,LETTERS),
            length,
            replace=TRUE
          ),
          collapse=""
        )
      }
      return(randomString)
    }
    rand.string = rand_string_create()
    
    eval(
      parse(
        text=paste(
          "dir.create('likelihood_results/DDlin_pos_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    
    tree.size = c()
    tree.number = c()
    sig2.1 = c()
    sig2.2 = c()
    sig2.3 = c()
    sig2.4 = c()
    root.1 = c()
    root.2 = c()
    slope.term.1 = c()
    slope.term.2 = c()
    slope.term.3 = c()
    slope.term.4 = c()
    est.sig2.1 = c()
    est.sig2.2 = c()
    est.sig2.3 = c()
    est.sig2.4 = c()
    est.root.1 = c()
    est.root.2 = c()
    est.slope.term.1 = c()
    est.slope.term.2 = c()
    est.slope.term.3 = c()
    est.slope.term.4 = c()
    log.likelihood = c()
    convergence = c()
    
    for (i in 1:length(tree.list)){
      for (j in 1:length(sig2.matrices)){
        for (k in 1:length(pars.list[[i]][[j]])){
          x = eval(
            parse(
              text=paste(
                "load('DDlin_pos_sim_tree_",
                names(tree.list)[i],
                "_root_sig2_",
                j,
                "_tip_sig2_",
                k,
                ".RData')",
                sep = ""
              )
            )
          )
          sim.results = get(x)
          
          for (l in Nsim){
            result = optim(
              par = c(0.5,0,0,1,1,0.5,-0.1,-0.05,-0.05,-0.1),
              fn = log_likelihood_mv_BM_DD,
              tree = tree.list[[i]][[l]],
              sim.value = sim.results[[l]],
              model = "DDlin",
              optim = TRUE
            )
            
            #calculate the position of the values in the data frame
            data.frame.position = (i-1)*(j-1)*(k-1)*l + (j-1)*(k-1)*l + (k-1)*l + l
            
            tree.size = c(tree.size,names(tree.list[i]))
            tree.number = c(tree.number,k)
            sig2.1 = c(sig2.1,DDpos.sig2.matrices[[j]][1])
            sig2.2 = c(sig2.2,DDpos.sig2.matrices[[j]][2])
            sig2.3 = c(sig2.3,DDpos.sig2.matrices[[j]][3])
            sig2.4 = c(sig2.4,DDpos.sig2.matrices[[j]][4])
            root.1 = c(root.1,0)
            root.2 = c(root.2,0)
            r.term.1 = c(r.term.1,pars.list[[i]][[j]][[k]][1])
            r.term.2 = c(r.term.2,pars.list[[i]][[j]][[k]][2])
            r.term.3 = c(r.term.3,pars.list[[i]][[j]][[k]][3])
            r.term.4 = c(r.term.4,pars.list[[i]][[j]][[k]][4])
            est.sig2.1 = c(est.sig2.1,result[[1]][1])
            est.sig2.2 = c(est.sig2.2,result[[1]][2])
            est.sig2.3 = c(est.sig2.3,result[[1]][3])
            est.sig2.4 = c(est.sig2.4,result[[1]][4])
            est.root.1 = c(est.root.1,result[[1]][5])
            est.root.2 = c(est.root.2,result[[1]][6])
            est.r.term.1 = c(est.r.term.1,result[[1]][7])
            est.r.term.2 = c(est.r.term.2,result[[1]][8])
            est.r.term.3 = c(est.r.term.3,result[[1]][9])
            est.r.term.4 = c(est.r.term.4,result[[1]][10])
            log.likelihood = c(log.likelihood,result[[2]])
            convergence = c(convergence,result[[4]])
          }
        }
      }
    }
    
    DDlin_pos_data = data.frame(
      tree.size,
      tree.number,
      sig2.1,
      sig2.2,
      sig2.3,
      sig2.4,
      root.1,
      root.2,
      slope.term.1,
      slope.term.2,
      slope.term.3,
      slope.term.4,
      est.sig2.1,
      est.sig2.2,
      est.sig2.3,
      est.sig2.4,
      est.root.1,
      est.root.2,
      est.slope.term.1,
      est.slope.term.2,
      est.slope.term.3,
      est.slope.term.4,
      log.likelihood,
      convergence
    )
    
    eval(
      parse(
        text=paste(
          "setwd('likelihood_results/DDlin_pos_likelihood_",
          rand.string,
          "')",
          sep=""
        )
      )
    )
    write.csv(DDlin_pos_data, file = "DDlin_pos_data.csv")
    
  } else {
    stop("Enter an applicable model")
  }
}

#max_likelihood_calculator_BM_OU_DD("OU")

no.cores = detectCores()
cl = makeCluster(no.cores, type=getClusterOption("type"),outfile='')
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cl, ls())
clusterEvalQ(cl,eval(parse(file = 'PhenotypicModel.R', n = 1)))
clusterEvalQ(cl,eval(parse(file = 'PhenotypicADiag.R', n = 1)))

clusterApply(
  cl,
  1:100,
  parallel_max_likelihood_calculator_BM_OU_DD,
  model = "DDexp_neg"
)

clusterApply(
  cl,
  1:100,
  parallel_max_likelihood_calculator_BM_OU_DD,
  model = "DDexp_pos"
)

clusterApply(
  cl,
  1:100,
  parallel_max_likelihood_calculator_BM_OU_DD,
  model = "DDlin_neg"
)

clusterApply(
  cl,
  1:100,
  parallel_max_likelihood_calculator_BM_OU_DD,
  model = "DDlin_pos"
)

stopCluster(cl)







































