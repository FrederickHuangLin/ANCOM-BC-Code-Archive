library(tidyverse)
source("ancom_bc_robust2.R")
source("sim_data_poi_gam_two_grp.R")

# The number of taxa, sampling depth, and sample size
n.taxa=500; samp.frac.var="large"; n.samp=c("20_30", "50_50")

# The proportion of differentially abundant taxa
prop.diff=c(0.05, 0.15, 0.25, 0.55, 0.65, 0.75)

# Set seeds
iterNum=100
abn.seed=seq(iterNum)

# Define the simulation parameters combinations
simparams=expand.grid(n.taxa, n.samp, prop.diff, abn.seed, samp.frac.var)
colnames(simparams)=c("n.taxa", "n.samp", "prop.diff", "abn.seed", "samp.frac.var")
simparams=simparams%>%mutate(obs.seed=abn.seed+1)
simparams=simparams%>%separate(col = n.samp, into = c("n.samp.grp1", "n.samp.grp2"), sep = "_")
simparams$n.samp.grp1=as.numeric(simparams$n.samp.grp1)
simparams$n.samp.grp2=as.numeric(simparams$n.samp.grp2)
simparams=simparams%>%arrange(n.taxa, n.samp.grp1, prop.diff, abn.seed, obs.seed)
simparams.list=apply(simparams, 1, paste0, collapse="_")

simparamslabels=c("n.taxa", "n.samp.grp1", "n.samp.grp2","prop.diff", 
                  "abn.seed", "samp.frac.var", "obs.seed")

library(doParallel)
library(foreach)

detectCores()
myCluster=makeCluster(8, type = "FORK")
registerDoParallel(myCluster)

start_time <- Sys.time()
simlist=foreach(i = simparams.list, .combine = 'cbind') %dopar% {
  # i = simparams.list[[1]]
  print(i)
  params = strsplit(i, "_")[[1]]
  names(params) <- simparamslabels
  
  # Paras for data generation
  n.taxa=as.numeric(params["n.taxa"])
  n.samp.grp1=as.numeric(params["n.samp.grp1"])
  n.samp.grp2=as.numeric(params["n.samp.grp2"])
  prop.diff=as.numeric(params["prop.diff"])
  abn.seed=as.numeric(params["abn.seed"])
  obs.seed=as.numeric(params["obs.seed"])
  samp.frac.var=params["samp.frac.var"]
  
  # Data generation
  test.dat=abn.tab.gen(n.taxa, n.samp.grp1, n.samp.grp2, low.abn=50, med.abn=200, high.abn=10000, prop.diff, 
                       abn.seed, obs.seed, struc.zero.prop=0.2, out.zero.prop=0.05, samp.frac.var)
  obs.abn=test.dat$obs.abn
  meta.data=cbind(Sample.ID=paste0("sub", seq(n.samp.grp1+n.samp.grp2)), 
                  group=rep(c(1, 2), c(n.samp.grp1, n.samp.grp2)))
  
  # Pre-processing
  feature.table=obs.abn; sample.var="Sample.ID"; group.var="group"; zero.cut=0.90; lib.cut=1000; neg.lb=FALSE
  pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, 
                                        group.var, zero.cut, lib.cut, neg.lb)
  feature.table=pre.process$feature.table
  group.name=pre.process$group.name
  group.ind=pre.process$group.ind
  struc.zero=pre.process$structure.zeros
  
  # Paras for ANCOM-BC
  grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
  tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05; bootNum = 100; boot.seed = 123
  
  # Run ANCOM-BC
  suppressWarnings(out <- try(ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero,
                                       adj.method, tol.EM, max.iterNum, perNum, alpha, bootNum, boot.seed), 
                              silent = TRUE))
  if (inherits(out, "try-error")) {
    FDR=NA; power=NA
  }else{
    res=cbind(out$res, diff.ind=test.dat$diff.taxa[rownames(out$feature.table)])
    
    # FDR
    FDR=ifelse(sum(res$diff.abn, na.rm = T)==0, 0, 
               sum(ifelse(res$diff.ind==0&res$diff.abn, 1, 0), na.rm = T)/
                 sum(res$diff.abn, na.rm = T))
    
    # Power
    power=sum(ifelse(res$diff.ind!=0&res$diff.abn, 1, 0), na.rm = T)/
      sum(res$diff.ind!=0, na.rm = T)
  }
  c(FDR, power)
}
end_time <- Sys.time()
end_time - start_time

stopCluster(myCluster)
write_csv(data.frame(simlist), "fdr_power_high_ancom_bc3.csv")