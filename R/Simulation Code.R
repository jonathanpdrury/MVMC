mv_sim_multiple = function(
  tree.sizes,
  sig2.matrices,
  root,
  pars.format = NULL,
  half.lives = NULL,
  OU.theta = NULL,
  Nsim = 100,
  models = c("BM","OU","MC","DD")
){
  require(phytools)
  ##check parameters are inputted correctly
  if (!is.numeric(tree.sizes)){
    stop("tree.sizes must be a vector of the number of tips for the trees in each test.")
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
  
  for (i in c("OU","MC","DD")){
    if (is.element(i,models)){
      models.present = TRUE
      break
    }
  }
  if (models.present){
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
    if (is.element("OU",models)){
      if (!is.numeric(OU.theta)||length(OU.theta)!=ncol(sig2.matrices[[1]])){
        stop("when using the OU model, a vector of theta values must be given.")
      }
    }
  }
  accepted.models = c("BM","OU","MC","DD")
  for (i in 1:length(models)){
    if (!is.element(models[i],accepted.models)){
      stop("the only accepted models are BM, OU, MC and DD")
    }
  }
  
  ##check if wanted and simulate the BM model
  if (is.element("BM",models)){
    for (i in 1:length(tree.sizes)){
      
      for (j in 1:length(sig2.matrices)){
        temp.data = list()
        
        for (k in 1:Nsim){
          tree = pbtree(n=tree.sizes[i])
          sim_data = sim_t_comp(
            tree,
            list(sig2.matrices[[k]]),
            root,
            10000,
            "BM"
          )
          temp.data[[k]] = sim_data
        }
        
        ##save data
        eval(
          parse(
            text=paste(
              "save(temp.data,file='BM_sim_tree_",
              tree.sizes[i],
              "_sig2_",
              j,
              ".RData')",
              sep=""
            )
          )
        )
      }
    }
  }
  
  ## if other models are used, construct a list of secondary parameters
  if (models.present){
    pars.list = list()
    for (i in 1:length(half.lives)){
      pars.list[[i]] = list()
      pars.1 = log(2)/half.lives[i]
      pars.2 = 0.5*pars.1
      cov.pars = 0.75*pars.2
      
      ##construct and save in a list an alpha matrix for each format
      
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
  }
  
  ##check if wanted and simulate OU model
  if (is.element("OU",models)){
    for (i in 1:length(tree.sizes)){
      
      for (j in 1:length(sig2.matrices)){
        
        for (k in 1:length(half.lives)){
          
          for (l in 1:length(pars.list[[k]])){
            temp.data = list()
          
            for (m in 1:Nsim){
              tree = pbtree(n=tree.sizes[i])
              sim_data = sim_t_comp(
                tree,
                list(sig2.matrices[[j]],pars.list[[k]][[l]],OU.theta),
                root,
                10000,
                "OU"
              )
              temp.data[[m]] = sim_data
            }
            
            ##save data
            eval(
              parse(
                text=paste(
                  "save(temp.data,file='OU_sim_tree",
                  tree.sizes[i],
                  "_sig2_",
                  j,
                  "_half_life_",
                  half.lives[k],
                  "_pars_",
                  l,
                  ".RData')",
                  sep=""
                )
              )
            )
          }
        }
      }
    }
  }
  
  ##check if wanted and simulate MC model
  if (is.element("MC",models)){
    for (i in 1:length(tree.sizes)){
      
      for (j in 1:length(sig2.matrices)){
        
        for (k in 1:length(half.lives)){
          
          for (l in 1:length(pars.list[[k]])){
            temp.data = list()
            
            for (m in 1:Nsim){
              tree = pbtree(n=tree.sizes[i])
              sim_data = sim_t_comp(
                tree,
                list(sig2.matrices[[j]],pars.list[[k]][[l]]),
                root,
                10000,
                "MC"
              )
              temp.data[[m]] = sim_data
            }
            
            ##save data
            eval(
              parse(
                text=paste(
                  "save(temp.data,file='MC_sim_tree",
                  tree.sizes[i],
                  "_sig2_",
                  j,
                  "_half_life_",
                  half.lives[k],
                  "_pars_",
                  l,
                  ".RData')",
                  sep=""
                )
              )
            )
          }
        }
      }
    }
  }
  
  ##check if wanted and simulate DD model
  if (is.element("DD",models)){
    for (i in 1:length(tree.sizes)){
      
      for (j in 1:length(sig2.matrices)){
        
        for (k in 1:length(half.lives)){
          
          for (l in 1:length(pars.list[[k]])){
            temp.data = list()
            
            for (m in 1:Nsim){
              tree = pbtree(n=tree.sizes[i])
              sim_data = sim_t_comp(
                tree,
                list(sig2.matrices[[j]],pars.list[[k]][[l]]),
                root,
                5000,
                "DD"
              )
              temp.data[[m]] = sim_data
            }
            
            ##save data
            eval(
              parse(
                text=paste(
                  "save(temp.data,file='DD_sim_tree",
                  tree.sizes[i],
                  "_sig2_",
                  j,
                  "_half_life_",
                  half.lives[k],
                  "_pars_",
                  l,
                  ".RData')",
                  sep=""
                )
              )
            )
          }
        }
      }
    }
  }
}



















