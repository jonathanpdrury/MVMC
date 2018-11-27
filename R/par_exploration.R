tree_020<-rcoal(20)
require(sim_t_comp_plot.r)

#rescale tree size to 10
tree_020$edge.length<-10/max(nodeHeights(tree_020))*tree_020$edge.length

#10 half lives (one half half every 1 unit)
sim_t_comp_plot(tree_020,c(1,-0.69),0,model="MC",plot=TRUE)
#obviously too extreme

#1 half lifes (one half half every 10 units)
sim_t_comp_plot(tree_020,c(1,-0.069),0,model="MC",plot=TRUE)

#3 half lives
sim_t_comp_plot(tree_020,c(1,-0.2079442),0,model="MC",plot=TRUE)

#5 half lives
sim_t_comp_plot(tree_020,c(1,-0.347),0,model="MC",plot=TRUE)


#for DDexp
log(0.01/sig2_0)/length(tree_020$tip.label)
#log(rate at tips/rate at root)/treesize

#if, for all negative slope parameter values we fix the rate at the tips to be 0.01
#sig2_0=0.1
sim_t_comp_plot(tree_020,c(1,-0.1151293),0,model="DDexp",plot=TRUE)

#sig2_0=0.5
sim_t_comp_plot(tree_020,c(1,-0.1956012),0,model="DDexp",plot=TRUE)

#sig2_0=1.0
sim_t_comp_plot(tree_020,c(1,-0.2302585),0,model="DDexp",plot=TRUE)

#if, for all positive slope parameter values we fix the rate at the root to be 0.01
#sig2_tip=0.1, positive slope (same thing as switching starting and endrate)
sim_t_comp_plot(tree_020,c(1,0.1151293),0,model="DDexp",plot=TRUE)

#sig2_tip=0.5, positive slope (same thing as switching starting and endrate)
sim_t_comp_plot(tree_020,c(1,0.1956012),0,model="DDexp",plot=TRUE)

#sig2_tip=1.0, positive slope (same thing as switching starting and endrate)
sim_t_comp_plot(tree_020,c(1,0.2302585),0,model="DDexp",plot=TRUE)


#for DDlin
(0.01-sig2_0)/length(tree_020$tip.label)
(rate at tips - rate at root)/treesize

#if, for all negative slope parameter values we fix the rate at the tips to be 0.01
#sig2_0=0.1
sim_t_comp_plot(tree_020,c(1,-0.0045),0,model="DDlin",plot=TRUE)

#sig2_0=0.5
sim_t_comp_plot(tree_020,c(1,-0.0245),0,model="DDlin",plot=TRUE)

#sig2_0=1
sim_t_comp_plot(tree_020,c(1,-0.0495),0,model="DDlin",plot=TRUE)

#if, for all positive slope parameter values we fix the rate at the root to be 0.01
#sig2_tip=0.1, positive slope (same thing as switching starting and endrate)
sim_t_comp_plot(tree_020,c(1,0.0045),0,model="DDlin",plot=TRUE)

#sig2_tip=0.5, positive slope (same thing as switching starting and endrate)
sim_t_comp_plot(tree_020,c(1,0.0245),0,model="DDlin",plot=TRUE)

#sig2_tip=1.0, positive slope (same thing as switching starting and endrate)
sim_t_comp_plot(tree_020,c(1,0.0495),0,model="DDlin",plot=TRUE)
