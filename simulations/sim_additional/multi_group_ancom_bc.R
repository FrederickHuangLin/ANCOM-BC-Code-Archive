library(tidyverse)
source("ancom_bc.R")
source("sim_data_poi_gam_multi_grp.R")

# The number of taxa, sampling depth, and sample size
n.taxa=1000; n.samp.grp=c("20 20 30 30", "50 50 50 50")

# The proportion of differentially abundant taxa
prop.diff=c(0.05, 0.15, 0.25)

# Set seeds
iterNum=100
abn.seed=seq(iterNum)

# Define the simulation parameters combinations
simparams=expand.grid(n.taxa, n.samp.grp, prop.diff, abn.seed)
colnames(simparams)=c("n.taxa", "n.samp.grp", "prop.diff", "abn.seed")
simparams=simparams%>%mutate(obs.seed=abn.seed+1)
simparams=simparams%>%arrange(n.taxa, n.samp.grp, prop.diff, abn.seed, obs.seed)
simparams.list=apply(simparams, 1, paste0, collapse="_")

simparamslabels=c("n.taxa", "n.samp.grp", "prop.diff", "abn.seed", "obs.seed")

library(doParallel)
library(foreach)

detectCores()
myCluster=makeCluster(10, type = "FORK")
registerDoParallel(myCluster)

start_time <- Sys.time()
simlist=foreach(i = simparams.list, .combine = 'cbind') %dopar% {
  # i = simparams.list[[1]]
  print(i)
  params = strsplit(i, "_")[[1]]
  names(params) <- simparamslabels
  
  # Paras for data generation
  n.taxa=as.numeric(params["n.taxa"])
  n.samp.grp=as.numeric(str_split(params["n.samp.grp"], " ")[[1]])
  prop.diff=as.numeric(params["prop.diff"])
  abn.seed=as.numeric(params["abn.seed"])
  obs.seed=as.numeric(params["obs.seed"])
  
  # Data generation
  low.abn=50; med.abn=200; high.abn=10000; struc.zero.prop=0.2; out.zero.prop=0.05
  test.dat=abn.tab.gen2(n.taxa, n.samp.grp, low.abn, med.abn, high.abn, prop.diff, 
                        abn.seed, obs.seed, struc.zero.prop, out.zero.prop)
  obs.abn=test.dat$obs.abn
  meta.data=cbind(Sample.ID=paste0("sub", seq(sum(n.samp.grp))), 
                  group=rep(c(1:length(n.samp.grp)), n.samp.grp))
  
  # Pre-processing
  feature.table=obs.abn; sample.var="Sample.ID"; group.var="group"
  zero.cut=0.90; lib.cut=1000; neg.lb=FALSE
  pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, 
                                        group.var, zero.cut, lib.cut, neg.lb)
  feature.table=pre.process$feature.table
  group.name=pre.process$group.name
  group.ind=pre.process$group.ind
  struc.zero=pre.process$structure.zeros
  
  # Paras for ANCOM-BC
  grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
  tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05
  
  # Run ANCOM-BC
  suppressWarnings(out <- try(ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero,
                                       adj.method, tol.EM, max.iterNum, perNum, alpha), 
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
write_csv(data.frame(simlist), "fdr_power_multi_group_ancom_bc.csv")



