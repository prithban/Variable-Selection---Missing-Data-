######################################################
### Code for running the parallel simulations on MI data in R
######################################################


############################
### Example 1:
### p = 10 , 50% missing ###
############################


######################
### Simulation Study
######################
rm(list=ls())
setwd("C:/Users/pbanerjee/Documents/MI_Simulation/Variable Selection/AOAS/Parallel_Results")


### Initial setup
library(parallel)
cl = makeCluster(detectCores())
clusterSetRNGStream(cl)
clusterEvalQ(cl, c(library(MASS), library(mice), library(mpath), library(grpreg), source("sim_aoas_data_generate.R"), source("degree_freedom.R"), source("degree_freedom_noweight.R"), source("sim_imput_combined.R")))

beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 5))

prop.miss = 50
dt1 = matrix(0,ncol=60)
#################
### sigma = 1 ###
sig = 1

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
t1 = proc.time()
run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))
(proc.time() - t1)
tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");



#################
### sigma = 2 ###
sig = 2

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


stopCluster(cl)




############################
### Example 2:
### p = 10 , 75% missing ###
############################


######################
### Simulation Study
######################
rm(list=ls())
setwd("C:/Users/pbanerjee/Documents/MI_Simulation/Variable Selection/AOAS/Parallel_Results")

library(parallel)
cl = makeCluster(detectCores())
clusterSetRNGStream(cl)
clusterEvalQ(cl, c(library(MASS), library(mice), library(mpath), library(grpreg), source("sim_aoas_data_generate.R"), source("degree_freedom.R"), source("degree_freedom_noweight.R"), source("sim_imput_combined.R")))

beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 5))

prop.miss = 75

#################
### sigma = 1 ###
sig = 1

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");



#################
### sigma = 2 ###
sig = 2

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


stopCluster(cl)




############################
### Example 3:
### p = 40 , 50% missing ###
############################


######################
### Simulation Study
######################
rm(list=ls())
setwd("C:/Users/pbanerjee/Documents/MI_Simulation/Variable Selection/AOAS/Parallel_Results")

library(parallel)
cl = makeCluster(detectCores())
clusterSetRNGStream(cl)
clusterEvalQ(cl, c(library(MASS), library(mice), library(mpath), library(grpreg), source("sim_aoas_data_generate.R"), source("degree_freedom.R"), source("degree_freedom_noweight.R"), source("sim_imput_combined.R")))

beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 35))

prop.miss = 50

#################
### sigma = 1 ###
sig = 1

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
t1 = proc.time()
run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))
(proc.time() - t1)
tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)
run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)
run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");



#################
### sigma = 2 ###
sig = 2

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


stopCluster(cl)




############################
### Example 4:
### p = 40 , 75% missing ###
############################


######################
### Simulation Study
######################
rm(list=ls())
setwd("C:/Users/pbanerjee/Documents/MI_Simulation/Variable Selection/AOAS/Parallel_Results")

library(parallel)
cl = makeCluster(detectCores())
clusterSetRNGStream(cl)
clusterEvalQ(cl, c(library(MASS), library(mice), library(mpath), library(grpreg), source("sim_aoas_data_generate.R"), source("degree_freedom.R"), source("degree_freedom_noweight.R"), source("sim_imput_combined.R")))

beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 35))

prop.miss = 75

#################
### sigma = 1 ###
sig = 1

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");



#################
### sigma = 2 ###
sig = 2

clusterExport(cl, list("beta.vec", "prop.miss", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 MCAR###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.2 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.5 MCAR###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 4 : rho = 0.5 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 5 : rho = 0.8 MCAR###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MCAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MCAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");


### Case 6 : rho = 0.8 MAR###
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, prop.miss = prop.miss, miss.type = "MAR", corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1"); dt1 = rbind(dt1,apply(tab1,2,mean));

write.csv(tab1, paste0("MAR_",corr,"_",length(beta.vec),"_",prop.miss,"_",rho,"_",sig,".csv")); rm("tab1");

stopCluster(cl)

dt1 = dt1[-1,]
write.csv(dt1,"combined results.csv")