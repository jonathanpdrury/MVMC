setClass(
    Class = "PhenotypicModel",
    representation = representation(
        name= "character",
        period = "numeric",
        aAGamma = "function",
        numbersCopy = "numeric",
        numbersPaste = "numeric",
        initialCondition = "function",
        paramsNames = "character",
        constraints = "function",
        params0 = "numeric",
        tipLabels = "character",
        comment = "character"
    ),
    prototype=prototype(
        name = "BMtest",
        period = c(0,1,2,3,4,5,6),
        aAGamma = function(i, params){
            functiona <- function(t){
                return(rep(0,i+1))
            }
            matrixA <- diag(0, i+1)
            functionGamma <- function(t){
                return(diag(params[1], i+1))
            }
            return(list(a=functiona, A=matrixA, Gamma=functionGamma))
        },
        numbersCopy = c(1, 1, 2, 1, 2, 5),
        numbersPaste = c(2, 3, 4, 5, 6, 7),
        initialCondition = function(params){
            return(list(mean=c(0,0), var=c(0,0)))
        },
        paramsNames = c("sigma"),
        constraints = function(params){
            return(params[1] > 0)
        },
        params0 = c(1),
        tipLabels = c("A", "B", "C", "D", "E", "F", "G"),
        comment = "Toy model defined by defaut"
    ),
    validity=function(object){
        if( length(object@numbersCopy) != length(object@numbersPaste) ){
            stop("[PhenotypicModel : validation] The sequence of positions of branching lineages and the sequence of new positions for the traits in the newly born lineages should have the same length.")
        }
        if( length(object@numbersCopy) != length(object@period) ){
            stop("[PhenotypicModel : validation] The sequence of positions of branching lineages and the sequence of time periods should have the same length.")
        }
        if( length(object@params0) != length(object@paramsNames) ){
            stop("[PhenotypicModel : validation] There should be the same number of defaut parameters and parameter names.")
        }
        return(TRUE)
    }
)


###################################
#    Getters and setters
###################################

setMethod(
    f="[",
    signature="PhenotypicModel",
    definition=function(x,i,j,drop){
        switch( EXPR=i,
                "name"={return(x@name)},
                "period"={return(x@period)},
                "aAGamma"={return(x@aAGamma)},
                "numbersCopy"={return(x@numbersCopy)},
                "numbersPaste"={return(x@numberPaste)},
                "initialCondition"={return(x@initialCondition)},
                "paramsNames"={return(x@paramsNames)},
                "constraints"={return(x@constraints)},
                "params0"={return(x@params0)},
                "tipLabels"={return(x@tipLabels)},
                "comment"={return(x@comment)},
                stop("This variable name does not exist !")
        )
    }
)

setReplaceMethod(
    f="[",
    signature="PhenotypicModel",
    definition=function(x,i,j,value){
        switch( EXPR=i,
                "name"={x@name <- value},
                "period"={x@period <- value},
                "aAGamma"={x@aAGamma <- value},
                "numbersCopy"={x@numbersCopy <- value},
                "numbersPaste"={x@numberPaste <- value},
                "initialCondition"={x@initialCondition <- value},
                "paramsNames"={x@paramsNames <- value},
                "constraints"={x@constraints <- value},
                "params0"={x@params0 <- value},
                "tipLabels"={x@tipLabels <- value},
                "comment"={x@comment <- value},
                stop("This variable name does not exist !")
        )
        validObject(x)
        return(x)
    }
)

###################################
#    Affichage
###################################

setMethod(
    f="print",
    signature="PhenotypicModel",
    definition=function(x, ...){
        cat("****************************************************************\n")
        cat("*** Object of Class PhenotypicModel *** \n")
        cat("*** Name of the model : ")
        print(x@name)
        cat("*** Parameters of the model : ")
        print(x@paramsNames)
        cat("*** Description : ")
        cat(x@comment)
        cat(paste("\n*** Periods : the model is cut into ", length(x@period), " parts. \n"))
        print(x@period)
        cat("*** Lineages branching (to be copied at the end of the corresponding period) :\n")
        print(x@numbersCopy)
        cat("*** Positions of the new trait at the end of each period :\n")
        print(x@numbersPaste)
        cat("*** Initial condition :\n")
        print(x@initialCondition)
        cat("*** Vectors a_i, A_i, Gamma_i on each period i : \n")
        print(x@aAGamma)
        cat("*** Constraints on the parameters : \n")
        print(x@constraints)
        cat("*** Defaut parameter values : ")
        print(x@params0)
        cat("*** Tip labels : \n")
        print(x@tipLabels)
        cat("****************************************************************\n")
    }
)

setMethod(
    f="show",
    signature="PhenotypicModel",
    definition=function(object){
        cat("****************************************************************\n")
        cat("*** Object of Class PhenotypicModel *** \n")
        cat("*** Name of the model : ")
        print(object@name)
        cat("*** Parameters of the model : ")
        print(object@paramsNames)
        cat("*** Description : ")
        cat(object@comment)
        cat(paste("\n*** Periods : the model is cut into ", length(object@period), " parts. \n"))
        cat("For more details on the model, call : print(PhenotypicModel)\n")
        cat("****************************************************************\n")
    }
)

###################################
#    Distribution
###################################

updateBranchingMatrixSigma <- function(Sigma, copy, paste){
    # copy of a branching lineage in the matrix of covariances 'Sigma'
    n = length(Sigma[1,])/2
    newSigma1 <- diag(0, n+1)
    newSigma2 <- diag(0, n+1)
    newSigma3 <- diag(0, n+1)
    newSigma4 <- diag(0, n+1)

    Sigma1=as.matrix(Sigma[1:n,1:n])
    Sigma2=as.matrix(Sigma[1:n,(n+1):length(Sigma[1,])])
    Sigma3=as.matrix(Sigma[(n+1):length(Sigma[1,]),1:n])
    Sigma4=as.matrix(Sigma[(n+1):length(Sigma[1,]),(n+1):length(Sigma[1,])])
    
    
    newSigma1[1:(paste-1),1:(paste-1)] <- Sigma1[1:(paste-1),1:(paste-1)]
    newSigma1[paste,1:(paste-1)] <- Sigma1[copy,1:(paste-1)]
    newSigma1[1:(paste-1),paste] <- Sigma1[1:(paste-1),copy]
    newSigma1[paste,paste] <- Sigma1[copy,copy]

    newSigma2[1:(paste-1),1:(paste-1)] <- Sigma2[1:(paste-1),1:(paste-1)]
    newSigma2[paste,1:(paste-1)] <- Sigma2[copy,1:(paste-1)]
    newSigma2[1:(paste-1),paste] <- Sigma2[1:(paste-1),copy]
    newSigma2[paste,paste] <- Sigma2[copy,copy]

    newSigma3[1:(paste-1),1:(paste-1)] <- Sigma3[1:(paste-1),1:(paste-1)]
    newSigma3[paste,1:(paste-1)] <- Sigma3[copy,1:(paste-1)]
    newSigma3[1:(paste-1),paste] <- Sigma3[1:(paste-1),copy]
    newSigma3[paste,paste] <- Sigma3[copy,copy]

    newSigma4[1:(paste-1),1:(paste-1)] <- Sigma4[1:(paste-1),1:(paste-1)]
    newSigma4[paste,1:(paste-1)] <- Sigma4[copy,1:(paste-1)]
    newSigma4[1:(paste-1),paste] <- Sigma4[1:(paste-1),copy]
    newSigma4[paste,paste] <- Sigma4[copy,copy]

    if(paste < n+1){
        newSigma1[(paste+1):(n+1),1:(paste-1)] <- Sigma1[paste:n,1:(paste-1)]
        newSigma1[(paste+1):(n+1),paste] <- Sigma1[paste:n,copy]
        newSigma1[1:(paste-1),(paste+1):(n+1)] <- Sigma1[1:(paste-1),paste:n]
        newSigma1[paste,(paste+1):(n+1)] <- Sigma1[copy,paste:n]
        newSigma1[(paste+1):(n+1),(paste+1):(n+1)] <- Sigma1[paste:n, paste:n]
        
        newSigma2[(paste+1):(n+1),1:(paste-1)] <- Sigma2[paste:n,1:(paste-1)]
        newSigma2[(paste+1):(n+1),paste] <- Sigma2[paste:n,copy]
        newSigma2[1:(paste-1),(paste+1):(n+1)] <- Sigma2[1:(paste-1),paste:n]
        newSigma2[paste,(paste+1):(n+1)] <- Sigma2[copy,paste:n]
        newSigma2[(paste+1):(n+1),(paste+1):(n+1)] <- Sigma2[paste:n, paste:n]
        
        newSigma3[(paste+1):(n+1),1:(paste-1)] <- Sigma3[paste:n,1:(paste-1)]
        newSigma3[(paste+1):(n+1),paste] <- Sigma3[paste:n,copy]
        newSigma3[1:(paste-1),(paste+1):(n+1)] <- Sigma3[1:(paste-1),paste:n]
        newSigma3[paste,(paste+1):(n+1)] <- Sigma3[copy,paste:n]
        newSigma3[(paste+1):(n+1),(paste+1):(n+1)] <- Sigma3[paste:n, paste:n]
        
        newSigma4[(paste+1):(n+1),1:(paste-1)] <- Sigma4[paste:n,1:(paste-1)]
        newSigma4[(paste+1):(n+1),paste] <- Sigma4[paste:n,copy]
        newSigma4[1:(paste-1),(paste+1):(n+1)] <- Sigma4[1:(paste-1),paste:n]
        newSigma4[paste,(paste+1):(n+1)] <- Sigma4[copy,paste:n]
        newSigma4[(paste+1):(n+1),(paste+1):(n+1)] <- Sigma4[paste:n, paste:n]
    }
	
	newSigma=rbind(cbind(newSigma1,newSigma2),cbind(newSigma3,newSigma4))
    return(newSigma)
}


setGeneric(
    name="getTipDistribution",
    def=function(object="PhenotypicModel", params="numeric", v="boolean"){standardGeneric("getTipDistribution")}
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicModel",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through ODE resolution (Riccati equation: dX = AX +XA +B) ***\n(Method working for any model)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){
            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                mean1<-mean[1:(length(mean)/2)]
                mean2<-mean[((length(mean)/2)+1):length(mean)]
                if( object@numbersPaste[i] <= length(mean1) ){
      #              mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
	  #					 mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)],mean[((length(mean)/2)+1):((length(mean)/2)+(object@numbersPaste[i]-1))], mean[((length(mean)/2)+object@numbersCopy[i])], mean[((length(mean)/2)+object@numbersPaste[i]):length(mean)] )
	  					mean1 <- c( mean1[1:(object@numbersPaste[i]-1)], mean1[object@numbersCopy[i]], mean1[object@numbersPaste[i]:length(mean1)] )
	  					mean2 <- c( mean2[1:(object@numbersPaste[i]-1)], mean2[object@numbersCopy[i]], mean2[object@numbersPaste[i]:length(mean2)] )
	  					mean<-c(mean1,mean2)
	  					
                }else{
      #             mean <- c( mean, mean[object@numbersCopy[i]] )
      #              mean <- c( mean[1:(length(mean)/2)], mean[object@numbersCopy[i]], mean[((length(mean)/2)+1):length(mean)],mean[object@numbersCopy[i]] )
                   mean1 <- c( mean1, mean1[object@numbersCopy[i]] )
                   mean2 <- c( mean2, mean2[object@numbersCopy[i]] )
	  			   mean<-c(mean1,mean2)
                }
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params)
            ai <- aAGammai$a
            Ai <- aAGammai$A
            Gammai <- aAGammai$Gamma
            n = length(mean)
            # We now need to build the ODE system such that dSigma/dt = -A Sigma - Sigma A + Gamma
            derivativeSigma <- function(t,y,params){
                Sigma = matrix(y,nrow=n)
                dSigma <- -Ai %*% Sigma - t(Sigma) %*% t(Ai) + Gammai(t) %*% t(Gammai(t))
                return(list(dSigma))
            }
            # And we build a second ODE system such that dm/dt = -Ai m + ai
            derivativemean <- function(t,y,params){
                return(list(-Ai %*% y + ai(t)))
            }

            # We update the vectors of means and covariances through their ODE system resolution
            times <- c(object@period[i], object@period[i+1])
#            print(Gammai(times[1]))
#            print(Ai)
#            print(ai)
#			if((object@period[i+1]-object@period[i])> 1e-15 ){
            if((object@period[i+1]-object@period[i])> 1e-14 ){
                mean  <- ode(mean, times, derivativemean)[2, 2:(n+1)]
                sigma <- ode(as.vector(Sigma), times, derivativeSigma)[2, 2:(n*n+1)]
            }
            Sigma = matrix(sigma,nrow=n)
        }
        
        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- c(object@tipLabels,object@tipLabels)
        rownames(Sigma) <- c(object@tipLabels,object@tipLabels)
        colnames(Sigma) <- c(object@tipLabels,object@tipLabels)

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)


setGeneric(
    name="getDataLikelihood",
    def=function(object="PhenotypicModel", data="numeric", params="numeric", v="boolean"){standardGeneric("getDataLikelihood")}
)

setMethod(
    f="getDataLikelihood",
    signature="PhenotypicModel",
    definition=function(object, data, params, v=FALSE){
        if(v){
            cat("*** Computing -log( likelihood ) of tip trait data under a given set of parameters ***\n")
        }

        if(object@constraints(params)){
            n <- length(data)
            tipdistribution <- getTipDistribution(object, params)
            V<-tipdistribution$Sigma
            data1<-data[1:(n/2)]
            data2<-data[((n/2)+1):n]
			data1<-data1[rownames(V)[1:(n/2)]]
			data2<-data2[rownames(V)[((n/2)+1):n]]
			data<-c(data1,data2)
			
  			op <- getOption("show.error.messages")
  			options(show.error.messages=FALSE)
			IV=try(solve(V))
  			options(show.error.messages=op)
  			if(class(IV)=="try-error"){
    			IV=pseudoinverse(V) 
  				if(max(IV)==0){return(Inf)}
  			}

            dataminusXT <- matrix(data - tipdistribution$mean, nrow=1)
            dataminusX <- matrix(data - tipdistribution$mean, ncol=1)

            ProdVectoriel = dataminusXT %*% IV %*% dataminusX
            
            calcul <-  (ProdVectoriel + determinant(V)$modulus + n*log(2*pi)) /2
			if(is.na(calcul) | is.infinite(calcul)){calcul=-1000000}
        }else{
            calcul <- -Inf
        }
        return(as.numeric(calcul))
    }
)


###################################
#    Parameter inferences
###################################

setGeneric(
    name="fitTipData",
    def=function(object="PhenotypicModel", data="numeric", params0="numeric", GLSstyle="logical"){standardGeneric("fitTipData")}
)

setMethod(
    f="fitTipData",
    signature="PhenotypicModel",
    definition=function(object, data, params0=NULL, GLSstyle=FALSE){
        cat("*** Fit of tip trait data ***\n")
        cat("Finding the maximum likelihood estimator of the parameters, before returning the likelihood and the inferred parameters...\n")
        beginning <- Sys.time()

        n <- length(data)

        # If params0 is not given, we use the 'params0' value contained in the model
        if(is.null(params0)){
            params0 <- object@params0
        }
        # In "GLS-style" mode, there is an analytical expression for the first two parameters, namely m0_1 and m0_2
        if(GLSstyle){
            params0 <- params0[3:length(params0)]
        }

        # computing the mean vector and variance matrix for the model, returns -log(likelihood) (a real number)
        toBeOptimized <- function(params){

            if(GLSstyle){paramsPrVerif <- c(0,0, params)}else{paramsPrVerif <- params}
            if(object@constraints(paramsPrVerif)){

                if(GLSstyle){
                    tipdistribution <- getTipDistribution(object, c(0,0,params))
                    
		            V<-tipdistribution$Sigma
		            #data<-data[rownames(V)]
		  			op <- getOption("show.error.messages")
		  			options(show.error.messages=FALSE)
					IV=try(solve(V))
		  			options(show.error.messages=op)
		  			if(class(IV)=="try-error"){
		    			IV=pseudoinverse(V) 
		  				if(max(IV)==0){return(Inf)}
		  			}
		
#                    I<-matrix(rep(1,n))
#					
#                    m0 <-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%as.matrix(data)[,1]
                    I=kronecker(diag(2),matrix(1,ncol=1,nrow=(n/2)))
					
					m0 = try(solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data)
					if(class(m0)=="try-error"){
					m0 = pseudoinverse(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data
					}
					m0_1 = m0[1]
					m0_2 = m0[2]
					
                    dataminusXT <- matrix(data - (I%*%m0), nrow=1)
                    dataminusX <- matrix(data - (I%*%m0), ncol=1)
#                    dataminusXT <- matrix(data - rep(m0, times=n), nrow=1)
#                    dataminusX <- matrix(data - rep(m0, times=n), ncol=1)
                    ProdVectoriel = dataminusXT %*% IV %*% dataminusX

                    calcul <-  (ProdVectoriel + determinant(V)$modulus+ n*log(2*pi)) /2
                    params <- c(m0_1,m0_2, params)
      
                }else{
                    calcul <- getDataLikelihood(object, data, params)
                }

            }else{
                calcul <- -Inf
            }

            return(calcul)
        }

        # looking for the argmin of -log(likelihood) (i.e. argmax of likelihood)
        optimisation <- optim(params0, toBeOptimized,control= list(maxit=20000))
        inferredParams <- optimisation$par
        # In GLS-style, we got all parameters except the first two, 'm0_1' and 'm0_2' that we compute through a last call to getTipDistribution
        if(GLSstyle){
            tipdistribution <- getTipDistribution(object, c(0,0,inferredParams))
		    V<-tipdistribution$Sigma
		  	op <- getOption("show.error.messages")
		  	options(show.error.messages=FALSE)
			IV=try(solve(V))
		  	options(show.error.messages=op)
		  	if(class(IV)=="try-error"){
		    	IV=pseudoinverse(V) 
		  		if(max(IV)==0){return(Inf)}
		  	}
		    #data<-data[rownames(V)]
            I=kronecker(diag(2),matrix(1,ncol=1,nrow=(n/2)))
					
			m0 = try(solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data)
			if(class(m0)=="try-error"){
					m0 = pseudoinverse(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data
			}
			m0_1 = m0[1]
			m0_2 = m0[2]
            
            inferredParams <- c(m0_1,m0_2, inferredParams)
        }
        names(inferredParams) <- object@paramsNames

        end <- Sys.time()
        cat("Computation time :", format(end-beginning), "\n")

        return(list(value = optimisation$value, inferredParams = inferredParams,convergence=optimisation$convergence))
    }
)

setGeneric(
    name="modelSelection",
    def=function(object="PhenotypicModel", data="numeric"){standardGeneric("modelSelection")}
)

setMethod(
    f="modelSelection",
    signature="PhenotypicModel",
    definition=function(object, data){
        cat("*** Model selection with tip trait data ***\n")
        cat("For each model in \"object\", fits the model and returns its AIC value in a recap table...\n")

        aic <- c()
        names <- c()
        for(model in object){
            fit <- fitTipData(model, data)
            aic <- c(aic, 2*length(model@params0)+2*fit$value )
            names <- c(names, model@name)
        }
        names(aic) <- names

        return(sort(aic))
    }
)

####################################
##    Simulation
####################################
#
#setGeneric(
#    name="simulateTipData",
#    def=function(object="PhenotypicModel", params="numeric", method="numeric"){standardGeneric("simulateTipData")}
#)
#
#setMethod(
#    f="simulateTipData",
#    signature="PhenotypicModel",
#    definition=function(object, params, method=3){
#        cat("*** Simulation of tip trait values ***\n")
#        if( method == 1 ){
#            cat("Computing first the tip distribution, and returning a simulated dataset drawn in this distribution...\n")
#            tipdistribution <- getTipDistribution(object, params)
#            X <- rmvnorm(1, tipdistribution$mean, tipdistribution$Sigma)
#
#        }else if( method == 2 ){        
#            cat("Simulating step-by-step the whole trajectory of a realization of the model and plotting the whole trajectory, before returning the tip data...\n")
#
#        }else{
#            cat("Simulating step-by-step the whole trajectory of a realization of the model, before returning only the tip data...\n")
#
#            initialCondition <- object@initialCondition(params)
#            X <- rnorm(length(initialCondition$mean), initialCondition$mean, initialCondition$var)
#            dt <- 0.001
#            sqrtdt <- sqrt(dt)
#        
#            for(i in 1:(length(object@period)-1)){
#                # If there is a branching event at the beginning of the period
#                if(object@numbersPaste[i] != 0){
#                    if( object@numbersCopy[i] < length(X) ){
#                        X <- c( X[1:object@numbersCopy[i]], X[object@numbersCopy[i]], X[(object@numbersCopy[i]+1):length(X)] )
#                    }else{
#                        X <- c( X, X[object@numbersCopy[i]] )
#                    }
#
#                }
#                # The time period is sliced
#                time <- seq( from = object@period[i], to = object@period[i+1], by = dt )
#                for(t in time){
#                    aAGammai <- object@aAGamma(i, params)  
#                    X <- X + (aAGammai$a(t) - aAGammai$A %*% X)*dt + sqrtdt* aAGammai$Gamma(t) %*% rnorm(length(X), 0, 1)
#                }
#            }
#        }
#	X <- matrix(data=X, ncol=1)
#        return(X)
#    }
#)
