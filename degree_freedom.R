##############################################################

####################################
#### Function for Degrees of freedom
####################################

############
##input are
##data of variables (X), response (y), weights considered (weight), 
##values of variable coefficients (beta.coef), additional parameter for penalty (a)
##value of regularization parameter (lambda), penalty types from LASSO, Bridge, SCAD, and MCP
############

############
##return are degrees of freedom and corresponding BIC
############

#############################################################

deg.frdm = function(X, y, weight, beta.coef, a = 0, lambda, type = c("lasso", "bridge", "scad", "mcp")){
	## index of nonzero beta coefficients
  ind = which(beta.coef!=0)
  
  ## initial setup
	X = as.matrix(X)
	n = length(y)
	w = as.matrix(weight)
	#w = diag(weight)
	
	## calculation of degrees of freedom and BIC
	if(length(ind) == 0){df = 0; bic.lam = log(sum(y^2)/n);
	}else{
		beta.nz = beta.coef[ind]
		x.beta = X[,ind]
		
		## LASSO penalty
		if(type == "lasso"){
			df = length(ind)
		}
		
		## Bridge penalty
		if(type == "bridge"){
		  p.lam = matrix(0, length(ind))
		  for(i in 1:length(beta.nz)){
		    if(abs(beta.nz[i]) > (lambda)^(1/a)){p.lam[i,] = a*lambda*(abs(beta.nz[i]))^(a-1)}
		    if(abs(beta.nz[i]) < (lambda)^(1/a)){p.lam[i,] = 0}
		  }
		  if(length(ind) == 1){ sig.lam = p.lam/abs(beta.nz)
		  }else{ sig.lam = diag(as.numeric(p.lam/abs(beta.nz)))}
		  df = sum(diag(x.beta%*%solve((t(x.beta)%*%w%*%x.beta + a*sig.lam), t(x.beta)%*%w)))
		}
		
		## SCAD penalty
		if(type == "scad"){
			p.lam = matrix(0, length(ind))
			for(i in 1:length(beta.nz)){
				if(abs(beta.nz[i]) < lambda){p.lam[i,] = lambda}
				if(abs(beta.nz[i]) >= lambda && abs(beta.nz[i]) < a*lambda){p.lam[i,] = (a*lambda - abs(beta.nz[i]))/(a-1)}
				if(abs(beta.nz[i]) >= a*lambda){p.lam[i,] = 0}
			}
			if(length(ind) == 1){ sig.lam = p.lam/abs(beta.nz)
			}else{ sig.lam = diag(as.numeric(p.lam/abs(beta.nz)))}
			df = sum(diag(x.beta%*%solve((t(x.beta)%*%w%*%x.beta + n*sig.lam), t(x.beta)%*%w)))
		}
		
		## MCP penalty
		if(type == "mcp"){
			p.lam = matrix(0, length(ind))
			for(i in 1:length(beta.nz)){
				if(abs(beta.nz[i]) <= a*lambda){p.lam[i,] = (a*lambda - abs(beta.nz[i]))/ a}
				if(abs(beta.nz[i]) >= a*lambda){p.lam[i,] = 0}
			}
			if(length(ind) == 1){ sig.lam = p.lam/abs(beta.nz)
			}else{ sig.lam = diag(as.numeric(p.lam/abs(beta.nz)))}
			sig.beta = (t(x.beta)%*%w%*%x.beta)/n
			sig = (t(X)%*%w%*%X)/n
			p.hat = matrix(0, length(ind), length(beta.coef));
			for(j in 1:length(ind)){p.hat[j,ind[j]] = 1}
			df = sum(diag(solve((sig.beta + sig.lam), p.hat%*%sig%*%t(p.hat))))
		}
		bic.lam = log(sum((y - X%*%beta.coef)^2)/n) + df*log(n)/n;
	}
	return(c(df,bic.lam))
}
