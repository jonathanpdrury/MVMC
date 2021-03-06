load('~/MVMC/DDexp_neg_data/tree_list.RData')
load('~/MVMC/DDexp_neg_data/DDexp_root_sig2_values.RData')
load('~/MVMC/DDexp_neg_data/DDexp_r_term_matrices.RData')

load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_3_tip_sig2_3.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_3_tip_sig2_2.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_3_tip_sig2_1.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_3.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_2.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_2_tip_sig2_1.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_1_tip_sig2_3.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_1_tip_sig2_2.RData')
load('~/MVMC/DDexp_neg_data/DDexp_neg_sim_tree_50_root_sig2_1_tip_sig2_1.RData')

source('~/MVMC/R/PhenotypicModel_MV_BETA.R', chdir = TRUE)
source('~/MVMC/R/PhenotypicADiag_MV_BETA.R', chdir = TRUE)

source('~/MVMC/R/createModel_BM_MV.R', chdir = TRUE)
source('~/MVMC/R/createModel_DDexp_MV.R', chdir = TRUE)

require(phytools)
require(corpcor)

res.mat<-matrix(nrow=5*3*3, ncol=20)
colnames(res.mat)<-c("tree","tree.size","sig2.1","sig2.2","sig2.cov","root.1", "root.2", "r.term.1", "r.term.2", "r.term.cov", "est.sig2.1", "est.sig2.2","est.sig2.cov","est.root.1","est.root.2","est.r.term.1","est.r.term.2","est.r.term.cov","lnL","convergence")
counter=1
for(i in 1:5){ #for 5 trees

	for(j in 1:3){ #for 3 sig2 matrices
	
		for(k in 1:3){ #for 3 slope matrices
		
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
		
		int<-c(i,50,sig2.matrices[[j]][1,1],sig2.matrices[[j]][2,2],sig2.matrices[[j]][2,1],0,0,pars.list[[2]][[j]][[k]][1,1],pars.list[[2]][[j]][[k]][2,2],pars.list[[2]][[j]][[k]][1,2],ddm.sig2_1,ddm.sig2_2,ddm.sig2_cov,ddm.m0_1,ddm.m0_2,ddm.r_1,ddm.r_2,ddm.r_cov,-ddm.fit$value, ddm.fit$convergence)
		
		res.mat[counter,]<-int
		print(int)
		write.csv(res.mat,file="ddm_rpanda_fits_int.csv")
		counter=counter+1
		}
	
	}	

}	

require(mvMORPH)

#now fitting BM.

res.mat<-matrix(nrow=5*3*3, ncol=17)
colnames(res.mat)<-c("tree","tree.size","sig2.1","sig2.2","sig2.cov","root.1", "root.2", "r.term.1", "r.term.2", "r.term.cov", "est.sig2.1", "est.sig2.2","est.sig2.cov","est.root.1","est.root.2","lnL","convergence")
counter=1
for(i in 1:5){ #for 5 trees

	for(j in 1:3){ #for 3 sig2 matrices
	
		for(k in 1:3){ #for 3 slope matrices
		
		tree<-tree.list[[2]][[i]]
		
		eval(parse(text=paste0("data<-DDexp_neg_sim_tree_50_root_sig2_",j,"_tip_sig2_",k,"_complete[[",i,"]]")))
		
		traits = matrix(nrow=length(tree$tip.label),ncol=2)
		for (l in 1:length(tree$tip.label)){
 		 	traits[l,1] = data[[l]][1]
  			traits[l,2] = data[[l]][2]
		}
		colnames(traits) = c("trait.1","trait.2")
		rownames(traits) = names(data)
		
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
		
		int<-c(i,50,sig2.matrices[[j]][1,1],sig2.matrices[[j]][2,2],sig2.matrices[[j]][2,1],0,0,pars.list[[2]][[j]][[k]][1,1],pars.list[[2]][[j]][[k]][2,2],pars.list[[2]][[j]][[k]][1,2],bmm.sig2_1,bmm.sig2_2,bmm.sig2_cov,bmm.m0_1,bmm.m0_2,bmm.fit$LogLik, bmm.fit$convergence)
		
		res.mat[counter,]<-int
		print(int)
		write.csv(res.mat,file="bmm_rpanda_fits_int.csv")
		counter=counter+1
		}
	
	}	

}	