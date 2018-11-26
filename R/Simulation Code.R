mv_sim_multiple = function(
  tree.list,
  tree.sizes,
  sig2.matrices,
  root,
  pars.format = NULL,
  half.lives = NULL,
  OU.theta = NULL,
  Nsim = 100,
  model
){
  require(phytools)
  ##check parameters are inputted correctly
  
  if (length(tree.sizes)!=length(tree.list)){
    stop("tree.sizes must be a numeric vector with each number being equal to the equivalent position on tree.list")
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
  
  if (!is.element(model,c("BM","OU","MC","DD"))){
    stop("the only accepted models are BM, OU, MC and DD")
  }
  models.present = FALSE
  for (i in c("OU","MC","DD")){
    if (is.element(i,model)){
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
    if (is.element("OU",model)){
      if (!is.numeric(OU.theta)||length(OU.theta)!=ncol(sig2.matrices[[1]])){
        stop("when using the OU model, a vector of theta values must be given.")
      }
    }
  }
  
  
  ##check if wanted and simulate the BM model
  if (model=="BM"){
    for (i in 1:length(tree.sizes)){
      
      for (j in 1:length(sig2.matrices)){
        temp.data = list()
        
        for (k in 1:Nsim){
          sim_data = sim_t_comp(
            tree.list[[i]][[k]],
            list(sig2.matrices[[k]]),
            root,
            10000,
            "BM"
          )
          temp.data[[k]] = sim_data
        }
        
        ##rename data
        eval(
          parse(
            text=paste(
              "assign('BM_sim_tree_",
              tree.sizes[i],
              "_sig2_",
              j,
              "',temp.data)",
              sep=""
            )
          )
        )
        
        ##save data
        eval(
          parse(
            text=paste(
              "save(BM_sim_tree_",
              tree.sizes[i],
              "_sig2_",
              j,
              ",file='BM_sim_tree_",
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
  } else {
    ## if other models are used, construct a list of secondary parameters
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
    
    if (model=="OU"){
      ##check if wanted and simulate OU model
      for (i in 1:length(tree.sizes)){
        
        for (j in 1:length(sig2.matrices)){
          
          for (k in 1:length(half.lives)){
            
            for (l in 1:length(pars.list[[k]])){
              temp.data = list()
              
              for (m in 1:Nsim){
                sim_data = sim_t_comp(
                  tree.list[[i]][[m]],
                  list(sig2.matrices[[j]],pars.list[[k]][[l]],OU.theta),
                  root,
                  10000,
                  "OU"
                )
                temp.data[[m]] = sim_data
              }
              
              ##rename data
              eval(
                parse(
                  text=paste(
                    "assign('OU_sim_tree_",
                    tree.sizes[i],
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    "',temp.data)",
                    sep=""
                  )
                )
              )
              
              ##save data
              eval(
                parse(
                  text=paste(
                    "save(OU_sim_tree_",
                    tree.sizes[i],
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    ",file='OU_sim_tree",
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
    } else if (model=="MC"){
      ##check if wanted and simulate MC model
      for (i in 1:length(tree.sizes)){
        
        for (j in 1:length(sig2.matrices)){
          
          for (k in 1:length(half.lives)){
            
            for (l in 1:length(pars.list[[k]])){
              temp.data = list()
              
              for (m in 1:Nsim){
                sim_data = sim_t_comp(
                  tree.list[[i]][[k]],
                  list(sig2.matrices[[j]],pars.list[[k]][[l]]),
                  root,
                  10000,
                  "MC"
                )
                temp.data[[m]] = sim_data
              }
              
              ##rename data
              eval(
                parse(
                  text=paste(
                    "assign('MC_sim_tree_",
                    tree.sizes[i],
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    "',temp.data)",
                    sep=""
                  )
                )
              )
              
              ##save data
              eval(
                parse(
                  text=paste(
                    "save(MC_sim_tree_",
                    tree.sizes[i],
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    ",file='MC_sim_tree",
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
    } else if (model=="DD"){
      ##check if wanted and simulate DD model
      for (i in 1:length(tree.sizes)){
        
        for (j in 1:length(sig2.matrices)){
          
          for (k in 1:length(half.lives)){
            
            for (l in 1:length(pars.list[[k]])){
              temp.data = list()
              
              for (m in 1:Nsim){
                sim_data = sim_t_comp(
                  tree.list[[i]][[k]],
                  list(sig2.matrices[[j]],pars.list[[k]][[l]]),
                  root,
                  5000,
                  "DD"
                )
                temp.data[[m]] = sim_data
              }
              
              ##rename data
              eval(
                parse(
                  text=paste(
                    "assign('DD_sim_tree_",
                    tree.sizes[i],
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    "',temp.data)",
                    sep=""
                  )
                )
              )
              
              ##save data
              eval(
                parse(
                  text=paste(
                    "save(DD_sim_tree_",
                    tree.sizes[i],
                    "_sig2_",
                    j,
                    "_half_life_",
                    half.lives[k],
                    "_pars_",
                    l,
                    ",file='DD_sim_tree",
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
}



















