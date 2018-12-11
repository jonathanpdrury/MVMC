###This function takes sim_t_comp bivariate and simulates it over multiple iterations with different parameter values.
###tree.list must take the form of a list the length of how many tree sizes you want to simulate, with values being a list of Nsim trees the same size.
###sig2.matrices is a list of each sig2 matrix to be simulated
###root is the initial value of each simulation (for now is constant among all simulations)
###pars.format and half lives are needed if using the OU or MC model
###half lives is a vector for the half life of pars.1 in the alpha or S term if OU or MC respectively.
###pars.format is used to build the rest of the matrix, 'pars.2' is half of 'pars.1' and 'cov.pars' is 3/4 of 'pars.2', a minus sign can be added in front to make the parameter negative
###If the model is OU, OU.theta is a vector defining the theta values, currently is fixed and cannot be varied.
###When using the DD model, define DD.root.rate and DD.tip.rate for the sig2 matrix at the root or tip and the R and slope term will be calculated for the DDexp and DDlin model respectively by using the equivalent sig2.matrices value as the rate at the root/tip not defined
###DD.tip.rate and DD.root.rate should be a list the length of sig2.matrices, with each value being a list containing matrices for each variation of sig2
###Nsegments is the default number of segments when simulating traits using sim_t_comp_bivariate
###Nsim is the amount of times each variant will be simulated
###model is a string defining the model used and can be 'BM', 'OU', 'MC', 'DDexp' or 'DDlin'
###If return.values is TRUE, this function will return a 'masterlist' of each simulation result, the trees used and the parameters calculated.
###If save.values is TRUE, this function will save each simulation to the working directory as well as the calculated parameter values and trees used.
###Both return.values and save.values can be TRUE in the same run of the simulation
mv_sim_multiple = function(
  tree.list,
  sig2.matrices,
  root,
  pars.format = NULL,
  half.lives = NULL,
  OU.theta = NULL,
  DD.root.rate = NULL,
  DD.tip.rate = NULL,
  Nsegments = 10000,
  Nsim = 100,
  model,
  return.values = FALSE,
  save.values = FALSE
){
  require(phytools)
  require(matlib)
  require(expm)
  ##check parameters are inputted correctly
  
  if (!return.values && !save.values){
    stop("You must select whether the values are saved or returned to the console.")
  }
  for (i in 1:length(tree.list)){
    if (length(tree.list[[i]])!=Nsim){
      stop("tree.list must be mulitple lists of length Nsim")
    }
  }
  if (!is.list(sig2.matrices)){
    stop("sig2.matrices must be a list of only matrices")
  }
  for (i in 1:length(sig2.matrices)){
    if (class(sig2.matrices[[i]])!="matrix"){
      stop("sig2.matrices must be a list of only matrices")
    }
  }
  if (length(root)!=ncol(sig2.matrices[[1]])){
    stop("root length and sig2 matrices do not match")
  }
  
  if (!is.element(model,c("BM","OU","MC","DDexp","DDlin"))){
    stop("the only accepted models are BM, OU, MC, DDexp and DDlin")
  }
  
  if (model=="OU"||model=="MC"){
    if (!is.numeric(half.lives)){
      stop("half lives must be a vector of each half life to test for the secondary term")
    }
    if (!is.list(pars.format)){
      stop("pars.format must be a list of only matrices")
    }
    for (i in 1:length(pars.format)){
      if (class(sig2.matrices[[i]])!="matrix"){
        stop("pars.format must be a list of only matrices")
      }
    }
    if (is.element("OU",model)){
      if (!is.numeric(OU.theta)||length(OU.theta)!=ncol(sig2.matrices[[1]])){
        stop("when using the OU model, a vector of theta values must be given.")
      }
    }
  } else if (model=="DDexp"||model=="DDlin"){
    if (!is.null(DD.tip.rate)){
      if (!is.list(DD.tip.rate)||length(DD.tip.rate)!=length(sig2.matrices)){
        stop("DD.tip.rate must be a list the length of sig2.matrices")
      }
      for (i in 1:length(DD.tip.rate)){
        if (!is.list(DD.tip.rate[[i]])){
          stop("each item of DD.tip.rate must be a list of matrices")
        }
        for (j in 1:length(DD.tip.rate[[i]])){
          if (class(DD.tip.rate[[i]][[j]])!="matrix"){
            stop("each item of DD.tip.rate must be a list of matrices")
          }
        }
      }
    } else if (!is.null(DD.root.rate)){
      if (!is.list(DD.root.rate)||length(DD.root.rate)!=length(sig2.matrices)){
        stop("DD.root.rate must be a list the length of sig2.matrices")
      }
      for (i in 1:length(DD.root.rate)){
        if (!is.list(DD.root.rate[[i]])){
          stop("each item of DD.root.rate must be a list of matrices")
        }
        for (j in 1:length(DD.root.rate[[i]])){
          if (class(DD.root.rate[[i]][[j]])!="matrix"){
            stop("each item of DD.root.rate must be a list of matrices")
          }
        }
      }
    } else {
      stop("you must either give the variants of sig2 matrices at the roots (if the change is positive) or at the tips (if the change is negative)")
    }
    if (!is.null(DD.tip.rate)&&!is.null(DD.root.rate)){
      stop("you must only give one of either DD.tip.rate or DD.root.rate.")
    }
  }
  
  ##if values are to be returned, make a masterlist
  if (return.values){
    masterlist = list()
    masterlist$simulated.values = list()
    masterlist$trees = tree.list
    masterlist$parameters = list()
    masterlist$parameters$sig2 = sig2.matrices
  }
  
  ##check if wanted and simulate the BM model
  if (model=="BM"){
    for (i in 1:length(tree.list)){
      ##if values are to be returned, construct data storage in masterlist
      if (return.values){
        eval(
          parse(
            text=paste(
              "masterlist$simulated.values$tree.",
              names(tree.list[i]),
              " = list()",
              sep = ""
            )
          )
        )
      }
      
      for (j in 1:length(sig2.matrices)){
        temp.data = list()
        
        ##simulate data
        for (k in 1:Nsim){
          sim_data = sim_t_comp(
            tree.list[[i]][[k]],
            list(sig2.matrices[[j]]),
            root,
            Nsegments,
            "BM"
          )
          temp.data[[k]] = sim_data
        }
        
        ##if values are to be returned, add to the masterlist simulations
        if (return.values){
          eval(
            parse(
              text=paste(
                "masterlist$simulated.values$tree.",
                names(tree.list[i]),
                "$sig2.",
                j,
                " = temp.data",
                sep = ""
              )
            )
          )
        }
        
        ##if values are to be saved, rename and save to the current directory
        if (save.values){
          eval(
            parse(
              text=paste(
                "assign('BM_sim_tree_",
                names(tree.list[i]),
                "_sig2_",
                j,
                "',temp.data)",
                sep = ""
              )
            )
          )
          
          eval(
            parse(
              text=paste(
                "save(BM_sim_tree_",
                names(tree.list[i]),
                "_sig2_",
                j,
                ",file='BM_sim_tree_",
                names(tree.list[i]),
                "_sig2_",
                j,
                ".RData')",
                sep = ""
              )
            )
          )
        }
      }
    }
  } else if (model=="OU"||model=="MC"){ ##check whether MC or OU model is wanted and simulate
    ##construct and save secondary parameters from the half lives
    pars.list = list()
    for (i in 1:length(half.lives)){
      pars.list[[i]] = list()
      pars.1 = log(2)/half.lives[i]
      pars.2 = 0.5*pars.1
      cov.pars = 0.75*pars.2
      
      for (j in 1:length(pars.format)){
        pars = pars.format[[j]]
        for (m in 1:length(pars)){
          if (pars[m]!=0){
            if (substring(pars[m],1,1)=="-"){
              if (substring(pars[m],2)=="pars.1"){
                pars[m]=-pars.1
              } else if (substring(pars[m],2)=="pars.2"){
                pars[m]=-pars.2
              } else if (substring(pars[m],2)=="cov.pars"){
                pars[m]=-cov.pars
              }
            } else if (pars[m]=="pars.1"){
              pars[m]=pars.1
            } else if (pars[m]=="pars.2"){
              pars[m]=pars.2
            } else if (pars[m]=="cov.pars"){
              pars[m]=cov.pars
            }
          } 
        }
        
        pars = matrix(as.numeric(pars),ncol=2)
        pars.list[[i]][[j]] = pars
      }
    }
    
    names(pars.list) = half.lives
    
    ##add calculated parameters to the masterlist
    if (model == "OU"){
      masterlist$parameters$alpha.parameters = pars.list
      masterlist$parameters$theta = OU.theta
    } else {
      masterlist$parameters$S.matrices = pars.list
    }
    
    ##simulate data using models
    for (i in 1:length(tree.list)){
      ##if values are to be returned, construct data storage in masterlist
      if (return.values){
        eval(
          parse(
            text=paste(
              "masterlist$simulated.values$tree.",
              names(tree.list[i]),
              " = list()",
              sep = ""
            )
          )
        )
      }
      
      for (j in 1:length(sig2.matrices)){
        if (return.values){
          eval(
            parse(
              text=paste(
                "masterlist$simulated.values$tree.",
                names(tree.list[i]),
                "$sig2.",
                j,
                "= list()",
                sep = ""
              )
            )
          )
        }
        
        for (k in 1:length(half.lives)){
          if (return.values){
            eval(
              parse(
                text=paste(
                  "masterlist$simulated.values$tree.",
                  names(tree.list[i]),
                  "$sig2.",
                  j,
                  "$half.life.",
                  half.lives[k],
                  "= list()",
                  sep = ""
                )
              )
            )
          }
          
          for (l in 1:length(pars.list[[k]])){
            temp.data = list()
            if (model=="OU"){
              sim.pars = list(sig2.matrices[[j]],pars.list[[k]][[l]],OU.theta)
            } else {
              sim.pars = list(sig2.matrices[[j]],pars.list[[k]][[l]])
            }
            
            for (m in 1:Nsim){
              sim_data = sim_t_comp(
                tree.list[[i]][[m]],
                sim.pars,
                root,
                Nsegments,
                model
              )
              temp.data[[m]] = sim_data
            }
            
            ##if values are to be returned, add to the masterlist simulations
            if (return.values){
              eval(
                parse(
                  text=paste(
                    "masterlist$simulated.values$tree.",
                    names(tree.list[i]),
                    "$sig2.",
                    j,
                    "$half.life.",
                    half.lives[k],
                    "$pars.format.",
                    l,
                    "= temp.data",
                    sep = ""
                  )
                )
              )
            }
            
            ##if values are to be saved, rename and save to the current directory
            if (save.values){
              eval(
                parse(
                  text=paste(
                    "assign('",
                    model,
                    "_sim_tree_",
                    names(tree.list[i]),
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    "',temp.data)",
                    sep = ""
                  )
                )
              )
              
              eval(
                parse(
                  text=paste(
                    "save(",
                    model,
                    "_sim_tree_",
                    names(tree.list[i]),
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    ",file=',",
                    model,
                    "_sim_tree",
                    names(tree.list[i]),
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    ".RData')",
                    sep = ""
                  )
                )
              )
            }
          }
        }
      }
    } 
  } else if (model=="DDexp"||model=="DDlin"){
    ##construct list to store parameters
    pars.list = list()
    
    ##make DD.variable to store DD.tip.rate or DD.root.rate for nicer code when constructing r.term
    if (!is.null(DD.tip.rate)){
      DD.variable = DD.tip.rate
    } else {
      DD.variable = DD.root.rate
    }
    
    ##check if wanted and simulate DDexp model
    for (i in 1:length(tree.list)){
      pars.list[[i]] = list()
      tree = tree.list[[i]]
      
      for (j in 1:length(sig2.matrices)){
        pars.list[[i]][[j]] = list()
        for (k in 1:(length(DD.variable))){
          temp.data = list()
          ##construct r.term/slope depending on the model used
          if (model=="DDexp"){
            if (!is.null(DD.tip.rate)){
              pars.list[[i]][[j]][[k]] = logm((inv(sig2.matrices[[j]])%*%DD.tip.rate[[j]][[k]]))/length(tree.list[[i]][[k]]$tip.label)
            } else {
              pars.list[[i]][[j]][[k]] = logm((inv(DD.root.rate[[j]][[k]])%*%sig2.matrices[[j]]))/length(tree.list[[i]][[k]]$tip.label)
            }
          } else if (model=="DDlin"){
            if (!is.null(DD.tip.rate)){
              pars.list[[i]][[j]][[k]] = (DD.tip.rate[[j]][[k]]-sig2.matrices[[j]])/length(tree.list[[i]][[k]]$tip.label)
            } else {
              pars.list[[i]][[j]][[k]] = (sig2.matrices[[j]]-DD.root.rate[[j]][[k]])/length(tree.list[[i]][[k]]$tip.label)
            }
          }
          
          for (m in 1:Nsim){
            sim_data = sim_t_comp(
              tree.list[[i]][[k]],
              list(sig2.matrices[[j]],pars.list[[i]][[j]][[k]]),
              root,
              Nsegments,
              "DDexp"
            )
            temp.data[[m]] = sim_data
          }
          
          if (!is.null(DD.tip.rate)){
            ##if values are to be returned, add to the masterlist simulations
            if (return.values){
              eval(
                parse(
                  text=paste(
                    "masterlist$simulated.values$tree.",
                    names(tree.list[i]),
                    "$root.sig2.",
                    j,
                    "$tip.sig2.",
                    k,
                    " = temp.data",
                    sep = ""
                  )
                )
              )
            }
            
            if (save.values){
              ##rename data
              eval(
                parse(
                  text=paste(
                    "assign('",
                    model,
                    "_sim_tree_",
                    names(tree.list[i]),
                    "_root_sig2_",
                    j,
                    "_tip_sig2_",
                    k,
                    "',temp.data)",
                    sep = ""
                  )
                )
              )
              
              ##save data
              eval(
                parse(
                  text=paste(
                    "save(",
                    model,
                    "_sim_tree_",
                    names(tree.list[i]),
                    "_root_sig2_",
                    j,
                    "_tip_sig2_",
                    k,
                    ",file='",
                    model,
                    "_sim_tree",
                    names(tree.list[i]),
                    "_root_sig2_",
                    j,
                    "_tip_sig2_",
                    k,
                    ".RData')",
                    sep = ""
                  )
                )
              )
            }
          } else {
            ##if values are to be returned, add to the masterlist simulations
            if (return.values){
              eval(
                parse(
                  text=paste(
                    "masterlist$simulated.values$tree.",
                    names(tree.list[i]),
                    "$tip.sig2.",
                    j,
                    "$root.sig2.",
                    k,
                    " = temp.data",
                    sep = ""
                  )
                )
              )
            }
            
            if (save.values){
              ##rename data
              eval(
                parse(
                  text=paste(
                    "assign('",
                    model,
                    "_sim_tree_",
                    names(tree.list[i]),
                    "_tip_sig2_",
                    j,
                    "_root_sig2_",
                    k,
                    "',temp.data)",
                    sep = ""
                  )
                )
              )
              
              ##save data
              eval(
                parse(
                  text=paste(
                    "save(",
                    model,
                    "_sim_tree_",
                    names(tree.list[i]),
                    "_tip_sig2_",
                    j,
                    "_root_sig2_",
                    k,
                    ",file='",
                    model,
                    "_sim_tree",
                    names(tree.list[i]),
                    "_tip_sig2_",
                    j,
                    "_root_sig2_",
                    k,
                    ".RData')",
                    sep = ""
                  )
                )
              )
            }
          }
        }
      }
    }
    #add names to pars.list
    names(pars.list) = tree.list
    
    if (return.values){
      if (model=="DDexp"){
        masterlist$parameters$r.matrix = pars.list
      } else {
        masterlist$parameters$slope.matrix = pars.list
      }
    }
  }
  
  ##if saving values, save parameters used.
  if (save.values){
    if (model=="BM"||model=="OU"||model=="MC"){
      eval(
        parse(
          text=paste(
            "save(sig2.matrices, file='",
            model,
            "_sig2_values.RData')",
            sep = ""
          )
        )
      )
      
      if (model=="OU"){
        eval(
          parse(
            text=paste(
              "save(pars.list,file='",
              model,
              "_alpha_matrices')",
              sep = ""
            )
          )
        )
        
        eval(
          parse(
            text=paste(
              "save(OU.theta,file='",
              model,
              "_theta_values.RData')",
              sep = ""
            )
          )
        )
      } else if (model=="MC"){
        eval(
          parse(
            text=paste(
              "save(pars.list,file='",
              model,
              "_S_matrices.RData')",
              sep = ""
            )
          )
        )
      }
    } else {
      eval(
        parse(
          text=paste(
            "save(sig2.matrices, file='",
            model,
            "_root_sig2_values.RData')",
            sep = ""
          )
        )
      )
      
      eval(
        parse(
          text=paste(
            "save(DD.tip.rate, file='",
            model,
            "_root_sig2_values.RData')",
            sep = ""
          )
        )
      )
      
      if (model=="DDexp"){
        eval(
          parse(
            text=paste(
              "save(pars.list, file='",
              model,
              "_r_term_matrices.RData')",
              sep = ""
            )
          )
        )
      } else {
        eval(
          parse(
            text=paste(
              "save(pars.list, file='",
              model,
              "_slope_matrices.RData')",
              sep = ""
            )
          )
        )
      }
    }
  }
  
  if (return.values){
    print("simulation complete")
    return(masterlist)
  } else {
    return(print("simulation complete"))
  }
}



















