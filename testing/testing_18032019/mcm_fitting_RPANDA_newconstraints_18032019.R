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

res.mat<-matrix(nrow=20, ncol=27)
colnames(res.mat)<-c("tree","tree.size","sig2.1","sig2.2","sig2.cov","root.1", "root.2", "r.term.1", "r.term.2", "r.term.cov", "DDexp.est.sig2.1", "DDexp.est.sig2.2","DDexp.est.sig2.cov","DDexp.est.root.1","DDexp.est.root.2","DDexp.est.r.term.1","DDexp.est.r.term.2","DDexp.est.r.term.cov","DDexp.lnL","DDexp.convergence","BM.est.sig2.1","BM.est.sig2.2","BM.est.sig2.cov","BM.est.root.1","BM.est.root.2","BM.lnL","BM.convergence")
counter=1
for(i in 1:20){ #for 5 trees

	for(j in 2:2){ #for 3 sig2 matrices
	
		for(k in 1:1){ #for 3 slope matrices
		
		tree<-tree.list[[2]][[i]]
		
		eval(parse(text=paste0("data<-DDexp_neg_sim_tree_50_root_sig2_",j,"_tip_sig2_",k,"_complete[[",i,"]]")))
		
		traits = matrix(nrow=length(tree$tip.label),ncol=2)
		for (l in 1:length(tree$tip.label)){
 		 	traits[l,1] = data[[l]][1]
  			traits[l,2] = data[[l]][2]
		}
		colnames(traits) = c("trait.1","trait.2")
		rownames(traits) = names(data)

		bmm<-createModel_BM_MV(tree)
		vcv_2<-getTipDistribution(bmm,params=c(0,0,0,0,0))$Sigma
		data.sorted<-c(traits[rownames(vcv_2)[1:50],1],traits[rownames(vcv_2)[1:50],2])

		ddm<-createModel_DDexp_MV(tree)
		ddm.fit<-fitTipData(ddm,data.sorted,GLSstyle=TRUE)
		
		ddm.sig2_1<-exp(ddm.fit$inferredParams[3])
		ddm.sig2_2<-exp(ddm.fit$inferredParams[4])
		ddm.sig2_cov<-(ddm.fit$inferredParams[7])
		ddm.r_1<-(ddm.fit$inferredParams[5])
		ddm.r_2<-(ddm.fit$inferredParams[6])
		ddm.r_cov<-(ddm.fit$inferredParams[8])
		ddm.m0_1<-(ddm.fit$inferredParams[1])
		ddm.m0_2<-(ddm.fit$inferredParams[2])
		
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

		int<-c(i,50,sig2.matrices[[j]][1,1],sig2.matrices[[j]][2,2],sig2.matrices[[j]][2,1],0,0,pars.list[[2]][[j]][[k]][1,1],pars.list[[2]][[j]][[k]][2,2],pars.list[[2]][[j]][[k]][1,2],ddm.sig2_1,ddm.sig2_2,ddm.sig2_cov,ddm.m0_1,ddm.m0_2,ddm.r_1,ddm.r_2,ddm.r_cov,-ddm.fit$value, ddm.fit$convergence,bmm.sig2_1,bmm.sig2_2,bmm.sig2_cov,bmm.m0_1,bmm.m0_2,bmm.fit$LogLik, bmm.fit$convergence)
		
		res.mat[counter,]<-int
		print(int)
		write.csv(res.mat,file="~/Dropbox/ddm_rpanda_fits_BETA_int_18032019_newconstraints.csv")
		counter=counter+1
		}
	
	}	

}	
