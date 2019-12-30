sim.imput = function(gen.size = 100, n = 200, rho = 0.5, sig = 1, beta.coef, prop.miss = 50, miss.type = c("Non", "MCAR", "MAR"), corr.type = c("linear", "auto"), wts.type = c("sc1", "sc2", "sc3"))
{
  miss_dat = sim.data(gen.size , n, rho, sig, beta.coef, prop.miss, miss.type, corr.type)
  pp = length(beta.coef)
  
  ### Check proportion of missing ###
  
  #miss_dat$missing.proportion
  
  tab = matrix(0,ncol = 12)
  
  for(i in 1:gen.size){ 		##for smaller number of iteration change
    
    tab1 = NULL
    final.dat = miss_dat$data[[i]]
    yy = final.dat[,1]
    true.index = 1:5
    dsgn.mtrx.org = final.dat[,-1]
    
    ### finding rows with missing data and computing number of missing in each row ###
    frac = apply(final.dat, 1, function(x){ length(which(is.na(x)))})
    ind = which(frac != 0); ind1 = which(frac == 0);
    
    final.dat<-final.dat[,-1]   ## drop the response from the data before imputation
    
    final.data=data.frame(final.dat)
    imp.dat=mice(final.data,m=5,method="pmm",seed=500)  ## "m" may be 8 or 10; discuss w/ Samiran-da
    
    ## the five data sets can be obtained by complete(imp.dat,1),
    ## complete(imp.dat,2),...,complete(imp.dat,5)
    
    ### appending with only missing rows ###
    final.dat.stack=rbind(complete(imp.dat,1)[ind1,],complete(imp.dat,1)[ind,],complete(imp.dat,2)[ind,],complete(imp.dat,3)[ind,],
                          complete(imp.dat,4)[ind,],complete(imp.dat,5)[ind,])
    
    
    ### Weights construction ###
    if(wts.type == "sc1"){
      n = length(frac)
      n1 = length(ind1); n2 = dim(final.dat.stack)[1] - n1;
      wts = c(rep(1/n, n1), rep(1/(n*5), n2))
    }
    if(wts.type == "sc2"){
      n = length(frac)
      n1 = length(ind1); n2 = dim(final.dat.stack)[1] - n1;
      wts = c(rep(1/n, n1), rep((1 - frac[ind]/pp)/(n*5), 5))
    }
    if(wts.type == "sc3"){
      
      prob.wts = function(x){
        mean.x = apply(x, 2, mean)
        dist.x = rowSums((sweep(x, 2, mean.x))^2)
        wts = (1/dist.x)/(sum(1/dist.x))
      }
      
      wts = rep(0, dim(final.dat.stack)[1])
      wts[1:length(ind1)] = 1
      ind.miss = length(ind1) + 1; ind.miss.len = length(ind);
      for(i in 1:ind.miss.len){
        u = ind.miss + (i-1)
        indx = c(u + (0:4)*ind.miss.len)
        dat.parse = final.dat.stack[indx,]
        wts[indx] = prob.wts(dat.parse)
      }
      wts = wts/n
    }
    
    ###Cholesky
    wt = chol(diag(wts))
    
    yy.new = c(yy[ind1], rep(yy[ind],5))
    yy.new1 = wt%*%yy.new
    x.new = wt%*%as.matrix(final.dat.stack)
    
    
    ### Lasso fit ###
    
    fit.lasso = glmreg(yy.new ~ ., data = data.frame(yy.new,final.dat.stack), weights = wts, family = "gaussian", penalty = "enet", alpha = 1)
    
    lambda.lasso = fit.lasso$lambda
    coef.mat = as.matrix(coef(fit.lasso))[-1,]
    res.lasso = matrix(0, nrow=length(lambda.lasso), ncol=3)
    
    res.lasso[,1] = lambda.lasso
    for(ii in 1:nrow(res.lasso)){
      res.lasso[ii,2:3] = deg.frdm(final.dat.stack, yy.new, weight = diag(wts), coef.mat[,ii], a=0, res.lasso[ii,1], type = "lasso")
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
    
    fit.scad = glmreg(yy.new ~ ., data = data.frame(yy.new,final.dat.stack), weights = wts, family = "gaussian", penalty = "snet", alpha = 1)
    
    lambda.scad = fit.scad$lambda
    coef.mat = as.matrix(coef(fit.scad))[-1,]
    res.scad = matrix(0, nrow=length(lambda.scad), ncol=3)
    
    res.scad[,1] = lambda.scad
    for(ii in 1:nrow(res.scad)){
      res.scad[ii,2:3] = deg.frdm(final.dat.stack, yy.new, weight = diag(wts), coef.mat[,ii], a=3.7, res.scad[ii,1], type = "scad")
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
    
    fit.mcp = glmreg(yy.new ~ ., data = data.frame(yy.new,final.dat.stack), weights = wts, family = "gaussian", penalty = "mnet", alpha = 1)
    
    lambda.mcp = fit.mcp$lambda
    coef.mat = as.matrix(coef(fit.mcp))[-1,]
    res.mcp = matrix(0, nrow=length(lambda.mcp), ncol=3)
    
    res.mcp[,1] = lambda.mcp
    for(ii in 1:nrow(res.mcp)){
      res.mcp[ii,2:3] = deg.frdm(final.dat.stack, yy.new, weight = diag(wts), coef.mat[,ii], a=3, res.mcp[ii,1], type = "mcp")
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
    
    fit.bridge = gBridge(x.new, yy.new1, family = "gaussian", alpha = 1, gamma = 0.5)
    
    lambda.bridge = fit.bridge$lambda
    coef.mat = as.matrix(coef(fit.bridge))[-1,]
    res.bridge = matrix(0, nrow=length(lambda.bridge), ncol=3)
    
    res.bridge[,1] = lambda.bridge
    for(ii in 1:nrow(res.bridge)){
      res.bridge[ii,2:3] = deg.frdm(final.dat.stack, yy.new, weight = diag(wts), coef.mat[,ii], a = 0.5, res.bridge[ii,1], type = "bridge")
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