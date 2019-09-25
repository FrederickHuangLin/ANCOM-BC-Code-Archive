library(tidyverse)
source("sim_data_poi_gam_two_grp.R")

# The number of taxa, library size, and sample size
n.taxa=1000; samp.frac.var="large"; n.samp=c("100_100")

# The proportion of differentially abundant taxa
prop.diff=0.75

# Set seeds
iterNum=100
abn.seed=seq(iterNum)

# Define the simulation parameters
simparams=expand.grid(n.taxa, n.samp, prop.diff, abn.seed, samp.frac.var)
colnames(simparams)=c("n.taxa", "n.samp", "prop.diff", "abn.seed", "samp.frac.var")
simparams=simparams%>%mutate(obs.seed=abn.seed+1)
simparams=simparams%>%separate(col = n.samp, into = c("n.samp.grp1", "n.samp.grp2"), sep = "_")
simparams=simparams%>%arrange(n.taxa, n.samp.grp1, prop.diff, abn.seed, obs.seed)
simparams.list=apply(simparams, 1, paste0, collapse="_")

simparamslabels=c("n.taxa", "n.samp.grp1", "n.samp.grp2","prop.diff", 
                  "abn.seed", "samp.frac.var", "obs.seed")

library(metagenomeSeq)
library(doParallel)
library(foreach)
simlist=foreach(i = simparams.list, .combine = 'cbind') %do% {
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
  low.abn=50; med.abn=200; high.abn=10000; struc.zero.prop=0.20; out.zero.prop=0.05
  test.dat=abn.tab.gen1(n.taxa, n.samp.grp1, n.samp.grp2, low.abn, med.abn, high.abn,
                        prop.diff, abn.seed, obs.seed, struc.zero.prop, out.zero.prop,
                        samp.frac.var)
  
  # Prepare data for metagenomeSeq
  meta.data=data.frame(group=rep(c(1, 2), c(n.samp.grp1, n.samp.grp2)))
  rownames(meta.data)=paste0("sub", seq(n.samp.grp1+n.samp.grp2))
  
  countdata=test.dat$obs.abn
  
  zero.threshold=0.90
  taxa.info.ind=apply(countdata, 1, function(x) sum(x==0)/(n.samp.grp1+n.samp.grp2))
  countdata=countdata[which(taxa.info.ind<zero.threshold), ]+1L
  
  # Run metagenomeSeq
  phenotypeData = AnnotatedDataFrame(meta.data)
  obj = newMRexperiment(countdata, phenoData=phenotypeData, featureData=NULL)
  
  # Calculating normalization factors
  obj = cumNorm(obj)
  
  # Zero-inflated Log-Normal mixture model
  pd = pData(obj)
  mod = model.matrix(~group, data = pd)
  
  suppressWarnings(fit <- try(fitFeatureModel(obj, mod)))
  if (inherits(fit, "try-error")) {
    power=NA; FDR=NA
  } else{
    out=MRcoefs(fit, number = nrow(countdata))
    out=data.frame(taxa=rownames(out), FDR=out$adjPvalues)
    out=out[match(rownames(countdata), as.character(out$taxa)), ]
    out$FDR[is.na(out$FDR)]=1
    
    res=data.frame(taxa=out$taxa,
                   diff.test=ifelse(out$FDR<0.05, 1, 0), 
                   diff.ind=test.dat$diff.taxa[which(taxa.info.ind<zero.threshold)])
    FDR=ifelse(sum(res$diff.test==1, na.rm = T)==0, 0,
               sum(ifelse(res$diff.ind==0&res$diff.test==1, 1, 0), na.rm = T)/
                 sum(res$diff.test==1, na.rm = T))
    power=sum(ifelse(res$diff.ind!=0&res$diff.test==1, 1, 0), na.rm = T)/
      sum(res$diff.ind!=0, na.rm = T)
  }
  
  c(FDR, power)
}
start_time <- Sys.time()

end_time <- Sys.time()
end_time - start_time

write_csv(data.frame(simlist), "fdr_power_zilg_large.csv")