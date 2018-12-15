source('DDexp.rescale.3.R')

mv.vcv.rescale.DDexp<-function(tree, sig2.matrix, slope.matrix){
	
	#transform rate matrix?
	
	block1<-vcv.rescale.DDexp()
	block2<-vcv.rescale.DDexp()
	block3<-vcv.rescale.DDexp()
	block4<-vcv.rescale.DDexp()
	
	#compile blocks
}