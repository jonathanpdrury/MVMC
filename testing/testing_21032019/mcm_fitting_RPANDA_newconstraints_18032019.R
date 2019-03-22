load('~/MVMC/testing/testing_18032019/tree_list (MC).RData')
load('~/MVMC/testing/testing_18032019/MC_S_matrices.RData')
load('~/MVMC/testing/testing_18032019/MC_sig2_values.RData')
load('~/MVMC/testing/testing_18032019/MC_sim_tree_20_sig2_2_half_life_2_pars_2.RData')

source('~/MVMC/R/PhenotypicModel_MV_BETA.R', chdir = TRUE)
source('~/MVMC/R/PhenotypicADiag_MV_BETA.R', chdir = TRUE)

source('~/MVMC/R/createModel_BM_MV.R', chdir = TRUE)
source('~/MVMC/R/createModel_MC_MV.R', chdir = TRUE)

require(phytools)
require(corpcor)

require(phytools)
require(corpcor)
require(mvMORPH)

#For MC, the trees are found in tree.list[[2]], the sig2 matrix is found in sig2.matrices[[2]] and the S term in pars.list[[2]][[2]].

res.mat<-matrix(nrow=20, ncol=39)
colnames(res.mat)<-c("tree","tree.size","sig2.1","sig2.2","sig2.cov","root.1", "root.2", "S.term.1", "S.term.2", "S.term.12cov","S.term.21cov", "MC.est.sig2.1", "MC.est.sig2.2","MC.est.sig2.cov","MC.est.root.1","MC.est.root.2","MC.est.S.term.1","MC.est.S.term.2","MC.est.S.term.12cov","MC.est.S.term.21cov","MC.lnL","MC.convergence","BM.est.sig2.1","BM.est.sig2.2","BM.est.sig2.cov","BM.est.root.1","BM.est.root.2","BM.lnL","BM.convergence","OU.est.sig2.1","OU.est.sig2.2","OU.est.sig2.cov","OU.est.root.1","OU.est.root.2","OU.est.alpha.1","OU.est.alpha.2","OU.est.alpha.cov","OU.lnL","OU.convergence")
counter=3
for(i in 3:20){ #for 20 trees

	for(j in 1:1){ #for 3 sig2 matrices
			
		tree<-tree.list[[1]][[i]]
		
		data<-MC_sim_tree_20_sig2_2_half_life_2_pars_2_complete[[i]]
		
		traits = matrix(nrow=length(tree$tip.label),ncol=2)
		for (l in 1:length(tree$tip.label)){
 		 	traits[l,1] = data[[l]][1]
  			traits[l,2] = data[[l]][2]
		}
		colnames(traits) = c("trait.1","trait.2")
		rownames(traits) = names(data)

		bmm<-createModel_BM_MV(tree)
		vcv_2<-getTipDistribution(bmm,params=c(0,0,0,0,0))$Sigma
		data.sorted<-c(traits[rownames(vcv_2)[1:20],1],traits[rownames(vcv_2)[1:20],2])

		
		bmm.fit<-mvBM(tree,traits,model="BM1")
		if(bmm.fit$convergence!=0){
			bmm.fit<-mvBM(tree,traits,model="BM1",method="pic")
		}
		if(bmm.fit$convergence!=0){
			bmm.fit<-mvBM(tree,traits,model="BM1",method="sparse")
		}
		
		bmm.sig2_1<-bmm.fit$sigma[1,1]
		bmm.sig2_2<-bmm.fit$sigma[2,2]
		bmm.sig2_cov<-bmm.fit$sigma[2,1]
		bmm.m0_1<-(bmm.fit$theta[1])
		bmm.m0_2<-(bmm.fit$theta[2])
		
		mcm<-createModel_MC_MV(tree)
		params0=c(0,0,log(sqrt(bmm.sig2_1)),log(sqrt(bmm.sig2_2)),bmm.sig2_cov,-0.1,-0.1,0,0)
		mcm.fit<-fitTipData(mcm,data.sorted,params0=params0,GLSstyle=TRUE)
###FIGURE OUT CRASH ON i==2
		
		mcm.sig2_1<-exp(mcm.fit$inferredParams[3])^2 + (mcm.fit$inferredParams[7]^2)
		mcm.sig2_2<-exp(mcm.fit$inferredParams[4])^2 + (mcm.fit$inferredParams[7]^2)
		mcm.sig2_cov<-(mcm.fit$inferredParams[7]) * (exp(mcm.fit$inferredParams[3])+exp(mcm.fit$inferredParams[4]))
		mcm.S_1<-(mcm.fit$inferredParams[5])
		mcm.S_2<-(mcm.fit$inferredParams[6])
		mcm.S_12cov<-(mcm.fit$inferredParams[8])
		mcm.S_21cov<-(mcm.fit$inferredParams[9])
		mcm.m0_1<-(mcm.fit$inferredParams[1])
		mcm.m0_2<-(mcm.fit$inferredParams[2])
		
		oum.fit<-mvOU(tree,traits,model="OU1")
		if(oum.fit$convergence!=0){
			oum.fit<-mvOU(tree,traits,model="OU1",method="sparse")
		}
		if(oum.fit$convergence!=0){
			oum.fit<-mvOU(tree,traits,model="OU1",method="inverse")
		}
		if(oum.fit$convergence!=0){
			oum.fit<-mvOU(tree,traits,model="OU1",method="pseudoinverse")
		}
		oum.sig2_1<-oum.fit$sigma[1,1]
		oum.sig2_2<-oum.fit$sigma[2,2]
		oum.sig2_cov<-oum.fit$sigma[2,1]
		oum.m0_1<-(oum.fit$theta[1])
		oum.m0_2<-(oum.fit$theta[2])
		oum.alpha_1<-oum.fit$alpha[1,1]
		oum.alpha_2<-oum.fit$alpha[2,2]
		oum.alpha_cov<-oum.fit$alpha[2,1]
		k=1
		int<-c(i,20,sig2.matrices[[j]][1,1],sig2.matrices[[j]][2,2],sig2.matrices[[j]][2,1],0,0,pars.list[[2]][[j]][[k]][1,1],pars.list[[2]][[j]][[k]][2,2],pars.list[[2]][[j]][[k]][1,2],pars.list[[2]][[j]][[k]][2,1],mcm.sig2_1,mcm.sig2_2,mcm.sig2_cov,mcm.m0_1,mcm.m0_2,mcm.S_1,mcm.S_2,mcm.S_12cov,mcm.S_21cov,-mcm.fit$value, mcm.fit$convergence,bmm.sig2_1,bmm.sig2_2,bmm.sig2_cov,bmm.m0_1,bmm.m0_2,bmm.fit$LogLik, bmm.fit$convergence,oum.sig2_1,oum.sig2_2,oum.sig2_cov,oum.m0_1,oum.m0_2,oum.alpha_1,oum.alpha_2,oum.alpha_cov,oum.fit$LogLik, oum.fit$convergence)
		
		res.mat[counter,]<-int
		print(int)
		write.csv(res.mat,file="~/Dropbox/mcm_rpanda_fits_int_21032019_newconstraints.csv")
		counter=counter+1
		}
	
	}	

}	
