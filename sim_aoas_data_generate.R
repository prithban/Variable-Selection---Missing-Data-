##############################################################

####################################
##### Function for generating simulated data
####################################

############
##input are
##number of MC datasets (gen.size), data length(n), correlation between variables (rho), sd of error (sig)
##coefficients for variables (beta.coef), proportion of missing (prop.miss), missing type (miss.type), and correlation structure (corr.type)
##linear means all nondiagonal are rho and auto means nondiagonal = rho^|i-j|
##missing types are Non missing, MCAR, and MAR
############

############
##return is list of datasets
############

#############################################################

sim.data = function(gen.size = 100, n = 200, rho = 0.5, sig = 1, beta.coef, prop.miss = 50, miss.type = c("Non", "MCAR", "MAR"), corr.type = c("linear", "auto"))
{
	#set.seed(100)
  
  ### Inittial setup
  p = length(beta.coef)
	mu.vec = rep(0,p)
	times = 1:p
	
	### Creating the sigma matrix for different correlation structure
	if(corr.type == "linear"){
	  H = matrix(1, p, p)
	  diag(H) = 0
	}
	if(corr.type == "auto"){
	  H = abs(outer(times, times, "-"))
	}

	sig.mat = rho^H

	### drawing design matrix X
	dsgn.mtrx = mvrnorm(n, mu.vec, sig.mat)
	
	### Complete variables
	gen.var = c(4, 5, 8) 
	
	### probability of missingness under different conditions
	if(prop.miss == 50){
	  if(miss.type == "MCAR"){
	    if(p ==10){prob.miss.mcar = 0.9}
	    if(p ==40){prob.miss.mcar = 0.9}
	  }
	  if(miss.type == "MAR"){
	    if(p ==10){prob.miss.mar = 0.96352}
	    if(p ==40){prob.miss.mar = 0.9961}
	  }
	}
	if(prop.miss == 75){
	  if(miss.type == "MCAR"){
	    if(p ==10){prob.miss.mcar = 0.9}
	    if(p ==40){prob.miss.mcar = 0.9}
	  }
	  if(miss.type == "MAR"){
	    if(p ==10){prob.miss.mar = 0.76308}
	    if(p ==40){prob.miss.mar = 0.97451}
	  }
	}
	
	### Initialize
	dat = list(); prop = 0;

	### Data with no missing values
	if(miss.type == "Non"){
	  for(i in 1:gen.size){
	    res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
	    final.dat=cbind(res,dsgn.mtrx)   ## full final data
	    colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
	    dat[[i]] = final.dat
	  }
	  return(data = dat)
	}
	
	### Data with missing type MCAR 
  if(miss.type == "MCAR"){
    
    for(i in 1:gen.size){
      res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
      n.miss = n*(prop.miss/100)  ##number of missing
      miss.index = sort(sample(1:n, n.miss))  ##missing index
	    comp.index = setdiff(1:n, miss.index) ##complete index
	    miss.var = setdiff(1:p, gen.var)  ##variables with missing values
      
	    ### Missing data generation
	    miss.mat = matrix(0, nrow = n.miss, ncol = p)
      miss.mat[,gen.var] = 1 
	    a1 = matrix(0,ncol = length(miss.var))
	    ii = 1
      while(ii <= n.miss){
	      samp1 = rbinom(length(miss.var), 1, prob.miss.mcar) ##creating missing data
	      if(sum(samp1) != length(miss.var)){
	        a1 = rbind(a1, samp1)
	        ii = ii + 1
	      }
      }
	    a1 = a1[-1,]
	    miss.mat[,miss.var] = a1  
	    rm(list = c("ii","a1"))

	    ### fianl design matrix with missing values
	    dsgn.mtrx.miss = matrix(0, nrow = n, ncol = p)
	    dsgn.mtrx.miss[comp.index,] = dsgn.mtrx[comp.index,]

	    for(ii in 1:n.miss){
	      for(jj in 1:p){
	        ifelse(miss.mat[ii,jj]==0,dsgn.mtrx.miss[miss.index[ii],jj]<-NA,dsgn.mtrx.miss[miss.index[ii],jj]<-dsgn.mtrx[miss.index[ii],jj])
	      }
	    }
	    rm(list=c("ii","miss.mat"))
      
	    ### Final data
	    final.dat=cbind(res,dsgn.mtrx.miss)   ## final data with missing covariates and response
	    colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
	    prop = prop + length(which(apply(is.na(final.dat) , 1, any)==TRUE))/n  ## checking proportion of missing observations

	    dat[[i]] = final.dat
  	}
    output = list(missing.proportion = prop/gen.size, data = dat)
    return(output)
  }
	
	### Data with missing type MAR
  if(miss.type == "MAR"){
    for(i in 1:gen.size){
      res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
	    miss.var = c(1, 2, 6, 7)  ## missing indexes
	    comp.index = setdiff(1:p, miss.var)
	    x = dsgn.mtrx[,gen.var]  ## design matrix for MAR
	    
	    ### Missing probabilities
      pr1 = 1/(1 + exp(-4.5 + x[,2] + 2*x[,3])) ##probability for missing for index 1
      pr2 = 1/(1 + exp(-3.5 - x[,1] - 2*x[,2] - x[,3] + res)) ##probability for missing for index 2
      pr3 = 1/(1 + exp(-2.15 + x[,1] - 1.5*x[,2] + 0.5*x[,3] + res))  ##probability for missing for index 6
      pr4 = 1/(1 + exp(-2 + 0.25*x[,1] + x[,2] + 0.5*x[,3]))  ##probability for missing for index 7

	    prob = cbind(pr1, pr2, pr3, pr4)
	    rm(list= c("pr1", "pr2", "pr3", "pr4"))

	    ### Missing data generation
	    miss.mat = matrix(0, nrow = n, ncol = length(miss.var))
      for(ii in 1:n){
        for(jj in 1:length(miss.var)){
          miss.mat[ii,jj] = rbinom(1, 1, prob[ii,jj]) ##generating missing
          #print(prob[ii,jj])
        }
      }
      rm(ii)
      miss.mat[miss.mat==0] <- NA
      miss.mat[which(miss.mat==1)] <- dsgn.mtrx[,miss.var][which(miss.mat==1)]
      #length(which(apply(is.na(miss.mat) , 1, any)==TRUE))
      
      res.var = setdiff(1:p, c(gen.var, miss.var))
      miss.mat1 = matrix(0, nrow = n, ncol = length(res.var))
      for(ii in 1:length(res.var)){
        miss.mat1[,ii] = rbinom(n, 1, prob.miss.mar)  
      }
      rm(ii)
      miss.mat1[miss.mat1==0] <- NA
      miss.mat1[which(miss.mat1==1)] <- dsgn.mtrx[,res.var][which(miss.mat1==1)]
      #length(which(apply(is.na(miss.mat1) , 1, any)==TRUE))
      
      ### Design matrix with missing data
      dsgn.mtrx.miss = matrix(0, nrow = n, ncol = p)
      dsgn.mtrx.miss[,gen.var] = dsgn.mtrx[,gen.var]
      dsgn.mtrx.miss[,miss.var] = miss.mat
      dsgn.mtrx.miss[,res.var] = miss.mat1

      ### Fianl data
      final.dat=cbind(res,dsgn.mtrx.miss)   ## final data with missing covariates and response
      colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
      prop = prop + length(which(apply(is.na(final.dat) , 1, any)==TRUE))/n  ## checking proportion of missing observations

      dat[[i]] = final.dat
    }
    output = list(missing.proportion = prop/gen.size, data = dat)
    return(output)
  }

}
