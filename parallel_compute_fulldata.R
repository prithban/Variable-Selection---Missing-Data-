######################################################
### Code for running the parallel simulations on full data in R
######################################################

rm(list=ls())
setwd("C:/Users/pbanerjee/Documents/MI_Simulation/AOAS/Parallel_Results")

library(parallel)
cl = makeCluster(detectCores())
clusterSetRNGStream(cl)
clusterEvalQ(cl, c(library(MASS), library(mice), library(mpath), library(grpreg), source("sim_aoas_data_generate.R"), source("degree_freedom.R"), source("sim_imput_fulldata.R")))

beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 5))

#################
### sigma = 1 ###
sig = 1

clusterExport(cl, list("beta.vec", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
t1 = proc.time()
run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, corr.type = corr))
(proc.time() - t1)
tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");



#################
### sigma = 2 ###
sig = 2

clusterExport(cl, list("beta.vec",  "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


stopCluster(cl)


######################
### Simulation Study on Full Data
######################
rm(list=ls())
setwd("C:/Users/pbanerjee/Documents/MI_Simulation/AOAS/Parallel_Results")

library(parallel)
cl = makeCluster(detectCores())
clusterSetRNGStream(cl)
clusterEvalQ(cl, c(library(MASS), library(mice), library(mpath), library(grpreg), source("sim_aoas_data_generate.R"), source("degree_freedom.R"), source("sim_imput_fulldata.R")))

beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 35))

#################
### sigma = 1 ###
sig = 1

clusterExport(cl, list("beta.vec", "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
t1 = proc.time()
run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec, corr.type = corr))
(proc.time() - t1)
tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");



#################
### sigma = 2 ###
sig = 2

clusterExport(cl, list("beta.vec",  "sig"))

#######################
## correlation type ###
#######################
corr = "linear"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


#######################
## correlation type ###
#######################
corr = "auto"

### Case 1 : rho = 0.2 Non###
rho = 0.2

clusterExport(cl, list("corr", "rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 2 : rho = 0.5 Non###
rho = 0.5

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


### Case 3 : rho = 0.8 Non###
rho = 0.8

clusterExport(cl, list("rho"))
clusterSetRNGStream(cl)

run1 = clusterApplyLB(cl, rep(13,8), function(nstart) sim.imput(gen.size = nstart, n = 200, rho = rho, sig = sig, beta.coef = beta.vec,  corr.type = corr))

tab1 = Reduce(rbind, run1)
rm("run1")

write.csv(tab1, paste0("Non_","fulldata",corr,"_",length(beta.vec),"_",rho,"_",sig,".csv")); rm("tab1");


stopCluster(cl)