##############################################################

####################################
##### Function for performing variable selection
##### on full simulated data corresponding to each penalty
####################################

############
##input are
##number of MC datasets (gen.size), data length(n), correlation between variables (rho), sd of error (sig)
##coefficients for variables (beta.coef), proportion of missing (prop.miss), and correlation structure (corr.type)
##linear means all nondiagonal are rho and auto means nondiagonal = rho^|i-j|
############

############
##return is a dataframe with MSE, Sensitivity, and Specificity corresponding to each penalty
############

#############################################################

sim.imput = function(gen.size = 100, n = 200, rho = 0.5, sig = 1, beta.coef, corr.type = c("linear", "auto"))
{
  non_miss_dat = sim.data(gen.size = gen.size, n = n, rho = rho, sig = sig, beta.coef = beta.coef, miss.type = "Non", corr.type = corr.type)
  pp = length(beta.coef)
  
  tab = matrix(0,ncol=12)
  
  for(i in 1:gen.size){ 		##for smaller number of iteration change
    
    tab1 = NULL
    
    ### extracting ith dataset from list
    
    dat = non_miss_dat[[i]]
    yy = dat[,1]
    true.index = 1:5
    wts = diag(rep(1,length(yy)))
    dsgn.mtrx.org = dat[,-1]
    
    
    ### LASSO fit ###
    #fit.lasso = glmnet(dsgn.mtrx.org, yy, family = "gaussian", alpha = 1)
    fit.lasso = glmreg(yy ~ ., data = data.frame(yy,dsgn.mtrx.org), family = "gaussian", penalty = "enet", alpha = 1, standardize = F)
    
    lambda.lasso = fit.lasso$lambda
    coef.mat = as.matrix(coef(fit.lasso))[-1,]
    res.lasso = matrix(0, nrow=length(lambda.lasso), ncol=3)
    
    res.lasso[,1] = lambda.lasso
    for(ii in 1:nrow(res.lasso)){
      res.lasso[ii,2:3] = deg.frdm(dsgn.mtrx.org, yy, wts, coef.mat[,ii], a=0, res.lasso[ii,1], type = "lasso")
    }
    
    min.ind = which.min(res.lasso[,3])
    #opt.lambda = res.lasso[min.ind,1]
    
    hat.beta.lasso= coef.mat[,min.ind]
    var.nonzero.cv.mcar = which(hat.beta.lasso!=0)
    
    (lasso.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA")))  ## true nonzero variable detection
    (lasso.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-lasso.true.index.cv.mcar)
    
    mse.lasso = t(beta.vec - hat.beta.lasso)%*%(beta.vec - hat.beta.lasso)
    sense.lasso = lasso.true.index.cv.mcar/ length(true.index)
    spec.lasso = 1 - lasso.false.positive.cv.mcar/(pp - length(true.index))
    
    tab1 = c(tab1, mse.lasso, sense.lasso, spec.lasso)
    
    ### SCAD fit ###
    
    #fit.scad = ncvreg(dsgn.mtrx.org, yy, family = "gaussian", penalty = "SCAD", alpha = 1)
    fit.scad = glmreg(yy ~ ., data = data.frame(yy,dsgn.mtrx.org), family = "gaussian", penalty = "snet", alpha = 1, gamma = 3.7, standardize = F)
    
    lambda.scad = fit.scad$lambda
    coef.mat = as.matrix(coef(fit.scad))[-1,]
    res.scad = matrix(0, nrow=length(lambda.scad), ncol=3)
    
    res.scad[,1] = lambda.scad
    for(ii in 1:nrow(res.scad)){
      res.scad[ii,2:3] = deg.frdm(dsgn.mtrx.org, yy, wts, coef.mat[,ii], a=3.7, res.scad[ii,1], type = "scad")
    }
    
    min.ind = which.min(res.scad[,3])
    #opt.lambda = res.scad[min.ind,1]
    
    hat.beta.scad= coef.mat[,min.ind]
    var.nonzero.cv.mcar = which(hat.beta.scad!=0)
    
    (scad.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA")))  ## true nonzero variable detection
    (scad.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-scad.true.index.cv.mcar)
    
    mse.scad = t(beta.vec - hat.beta.scad)%*%(beta.vec - hat.beta.scad)
    sense.scad = scad.true.index.cv.mcar/ length(true.index)
    spec.scad = 1 - scad.false.positive.cv.mcar/(pp - length(true.index))
    
    tab1 = c(tab1, mse.scad, sense.scad, spec.scad)
    
    ### MCP fit ###
    
    #fit.mcp = ncvreg(dsgn.mtrx.org, yy, family = "gaussian", penalty = "MCP", alpha = 1)
    fit.mcp = glmreg(yy ~ ., data = data.frame(yy,dsgn.mtrx.org), family = "gaussian", penalty = "mnet", alpha = 1, gamma = 3, standardize = F)
    
    lambda.mcp = fit.mcp$lambda
    coef.mat = as.matrix(coef(fit.mcp))[-1,]
    res.mcp = matrix(0, nrow=length(lambda.mcp), ncol=3)
    
    res.mcp[,1] = lambda.mcp
    for(ii in 1:nrow(res.mcp)){
      res.mcp[ii,2:3] = deg.frdm(dsgn.mtrx.org, yy, wts, coef.mat[,ii], a=3, res.mcp[ii,1], type = "mcp")
    }
    
    min.ind = which.min(res.mcp[,3])
    #opt.lambda = res.mcp[min.ind,1]
    
    hat.beta.mcp= coef.mat[,min.ind]
    var.nonzero.cv.mcar = which(hat.beta.mcp!=0)
    
    (mcp.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA")))  ## true nonzero variable detection
    (mcp.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-mcp.true.index.cv.mcar)
    
    mse.mcp = t(beta.vec - hat.beta.mcp)%*%(beta.vec - hat.beta.mcp)
    sense.mcp = mcp.true.index.cv.mcar/ length(true.index)
    spec.mcp = 1 - mcp.false.positive.cv.mcar/(pp - length(true.index))
    
    tab1 = c(tab1, mse.mcp, sense.mcp, spec.mcp)
    
    ### Bridge fit ###
    
    fit.bridge = gBridge(dsgn.mtrx.org, yy, family = "gaussian", alpha = 1, gamma = 0.5)
    
    lambda.bridge = fit.bridge$lambda
    coef.mat = as.matrix(coef(fit.bridge))[-1,]
    res.bridge = matrix(0, nrow=length(lambda.bridge), ncol=3)
    
    res.bridge[,1] = lambda.bridge
    for(ii in 1:nrow(res.bridge)){
      res.bridge[ii,2:3] = deg.frdm(dsgn.mtrx.org, yy, wts, coef.mat[,ii], a = 0.5, res.bridge[ii,1], type = "bridge")
    }
    
    min.ind = which.min(res.bridge[,3])
    #opt.lambda = res.bridge[min.ind,1]
    
    hat.beta.bridge= coef.mat[,min.ind]
    var.nonzero.cv.mcar = which(hat.beta.bridge!=0)
    
    (bridge.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA")))  ## true nonzero variable detection
    (bridge.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-bridge.true.index.cv.mcar)
    
    mse.bridge = t(beta.vec - hat.beta.bridge)%*%(beta.vec - hat.beta.bridge)
    sense.bridge = bridge.true.index.cv.mcar/ length(true.index)
    spec.bridge = 1 - bridge.false.positive.cv.mcar/(pp - length(true.index))
    
    tab1 = c(tab1, mse.bridge, sense.bridge, spec.bridge)
    
    tab = rbind(tab,tab1)
    
  }
  tab = tab[-1,]
  colnames(tab) = c("mse.lasso","sensitivity.lasso","specificity.lasso","mse.scad","sensitivity.scad","specificity.scad","mse.mcp","sensitivity.mcp","specificity.mcp","mse.bridge","sensitivity.bridge","specificity.bridge")
  
  return(tab)
}