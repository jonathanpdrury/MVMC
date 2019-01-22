source('PhenotypicModel.R') #from RPANDA (Marc Manceau)
source('PhenotypicADiagonal.R') #from RPANDA (Marc Manceau)
source('DDexp_nogeo_ADiag.R') #from RPANDA 
source('DDlin_nogeo_ADiag.R') #from RPANDA 

require(phytools)

mv.vcv.rescale.DD<-function(tree, sig2.matrix, slope.matrix,model=c("DDexp","DDlin")){
	
	#fit model object
	if(model="DDexp"){	
	dd.ob<-createModel_DDexp(tree)
	}else if(model="DDlin"){
	dd.ob<-createModel_DDlin(tree)
	}
	else{ stop("model must be either 'DDexp' or 'DDlin'")}
	
	block1<-getTipDistribution(dd.ob,params=c(0,log(sig2.matrix[1,1]),slope.matrix[1,1]))$Sigma
	block2<-getTipDistribution(dd.ob,params=c(0,log(sig2.matrix[1,2]),slope.matrix[1,2]))$Sigma
	block3<-getTipDistribution(dd.ob,params=c(0,log(sig2.matrix[2,1]),slope.matrix[2,1]))$Sigma
	block4<-getTipDistribution(dd.ob,params=c(0,log(sig2.matrix[2,2]),slope.matrix[2,2]))$Sigma
	
	#compile blocks
	
	mv.vcv<-rbind(cbind(block1,block2),cbind(block3,block4))
	return(mv.vcv)
}