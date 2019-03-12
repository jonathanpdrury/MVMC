##Checking existing fitting tools

fits<-read.csv("~/MVMC/testing/testing_09032019/ddm_rpanda_fits_int_08032019_newconstraints_testing.csv")

#it doesn't appear that DDexp is chosen more with more extreme slopes:
table(fits$Model.chosen ,fits$level,fits$scenario)


#sig2 1 estimated ok in unconstrained version:
hist(fits$DDexp.est.sig2.1_t,breaks=1000,xlim=c(0,5))
abline(v=1,col="red",lty=2,lwd=2)

#sig2 1 not estimated ok in constrained version:
hist(fits$DDexp_constrained.est.sig2.1_t,breaks=100)
abline(v=1,col="red",lty=2,lwd=2)

#sig2 2 estimated more or less ok in unconstrained version:
boxplot(DDexp.est.sig2.2_t~ sig2.2,fits)
abline(h=unique(fits$sig2.2),lwd=2,lty=2,col="red")

#sig2 2 slightly underestimated in constrained version:
boxplot(DDexp_constrained.est.sig2.2_t~ sig2.2,fits)
abline(h=unique(fits$sig2.2),lwd=2,lty=2,col="red")

#sig2 cov estimated more or less ok in unconstrained version:
boxplot(DDexp.est.sig2.cov_t~ sig2.cov,fits,ylim=c(-3,3))
abline(h=unique(fits$sig2.cov),lwd=2,lty=2,col="red")

#sig2 cov estimated more or less ok in constrained version:
boxplot(DDexp_constrained.est.sig2.cov_t~ sig2.cov,fits,ylim=c(-2,3))
abline(h=unique(fits$sig2.cov),lwd=2,lty=2,col="red")

#r term 1 severely overestimated
boxplot(DDexp.est.r.term.1~ r.term.1,fits)
abline(h=unique(fits$r.term.1),lwd=2,lty=2,col="red")

boxplot(DDexp_constrained.est.r.term.1~ r.term.1,fits)
abline(h=unique(fits$r.term.1),lwd=2,lty=2,col="red")

#r term 2 severely overestimated
boxplot(DDexp.est.r.term.2~ r.term.2,fits)
abline(h=unique(fits$r.term.2),lwd=2,lty=2,col="red")

boxplot(DDexp_constrained.est.r.term.2~ r.term.2,fits)
abline(h=unique(fits$r.term.2),lwd=2,lty=2,col="red")

boxplot(log(abs(DDexp.est.r.term.1)/abs(DDexp.est.r.term.2))~ log(abs(r.term.1)/abs(r.term.2)),fits)
abline(h=unique(log(abs(fits$r.term.1)/abs(fits$r.term.2))),lwd=2,lty=2,col="red")

boxplot(log(abs(DDexp_constrained.est.r.term.1)/abs(DDexp_constrained.est.r.term.2))~ log(abs(r.term.1)/abs(r.term.2)),fits)
abline(h=unique(log(abs(fits$r.term.1)/abs(fits$r.term.2))),lwd=2,lty=2,col="red")

#r term cov less terrible?

boxplot(DDexp.est.r.term.cov~ r.term.cov,fits)
abline(h=unique(fits$r.term.cov),lwd=2,lty=2,col="red")

#no cov in slope for constrained version