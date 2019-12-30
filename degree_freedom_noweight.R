##############################################################

####################################
#### Function for Degrees of freedom
####################################

############
##input are
##data of variables (X), response (y),
##values of variable coefficients (beta.coef), additional parameter for penalty (a)
##value of regularization parameter (lambda), penalty types from LASSO, Bridge, SCAD, and MCP
############

############
##return are degrees of freedom and corresponding BIC
############

#############################################################


deg.frdm.nw = function(X, y, beta, a = 0, lambda, type = c("lasso", "bridge", "scad", "mcp")){
  ind = which(beta!=0)
  X = as.matrix(X)
  n = length(y)
    if(length(ind) == 0){df = 0; bic.lam = log(sum(y^2)/n);
  }else{
    beta.nz = beta[ind]
    x.beta = X[,ind]
    
    ## LASSO
    if(type == "lasso"){
      df = length(ind)
    }
    
    ## Bridge
    if(type == "bridge"){
      p.lam = matrix(0, length(ind))
      for(i in 1:length(beta.nz)){
        if(abs(beta.nz[i]) > (lambda)^(1/a)){p.lam[i,] = a*lambda*(abs(beta.nz[i]))^(a-1)}
        if(abs(beta.nz[i]) < (lambda)^(1/a)){p.lam[i,] = 0}
      }
      if(length(ind) == 1){ sig.lam = p.lam/abs(beta.nz)
      }else{ sig.lam = diag(as.numeric(p.lam/abs(beta.nz)))}
      df = sum(diag(x.beta%*%solve((t(x.beta)%*%x.beta + a*sig.lam), t(x.beta))))
    }
    
    ## SCAD
    if(type == "scad"){
      p.lam = matrix(0, length(ind))
      for(i in 1:length(beta.nz)){
        if(abs(beta.nz[i]) < lambda){p.lam[i,] = lambda}
        if(abs(beta.nz[i]) >= lambda && abs(beta.nz[i]) < a*lambda){p.lam[i,] = (a*lambda - abs(beta.nz[i]))/(a-1)}
        if(abs(beta.nz[i]) >= a*lambda){p.lam[i,] = 0}
      }
      if(length(ind) == 1){ sig.lam = p.lam/abs(beta.nz)
      }else{ sig.lam = diag(as.numeric(p.lam/abs(beta.nz)))}
      df = sum(diag(x.beta%*%solve((t(x.beta)%*%x.beta + n*sig.lam), t(x.beta))))
    }
    
    ## MCP
    if(type == "mcp"){
      p.lam = matrix(0, length(ind))
      for(i in 1:length(beta.nz)){
        if(abs(beta.nz[i]) <= a*lambda){p.lam[i,] = (a*lambda - abs(beta.nz[i]))/ a}
        if(abs(beta.nz[i]) >= a*lambda){p.lam[i,] = 0}
      }
      if(length(ind) == 1){ sig.lam = p.lam/abs(beta.nz)
      }else{ sig.lam = diag(as.numeric(p.lam/abs(beta.nz)))}
      sig.beta = (t(x.beta)%*%x.beta)/n
      sig = (t(X)%*%X)/n
      p.hat = matrix(0, length(ind), length(beta));
      for(j in 1:length(ind)){p.hat[j,ind[j]] = 1}
      df = sum(diag(solve((sig.beta + sig.lam), p.hat%*%sig%*%t(p.hat))))
    }
    bic.lam = log(sum((y - X%*%beta)^2)/n) + df*log(n)/n;
  }
  return(c(df,bic.lam))
}
