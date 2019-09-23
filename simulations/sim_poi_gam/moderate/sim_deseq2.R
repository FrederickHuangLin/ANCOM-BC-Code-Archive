library(tidyverse)
source("sim_data_poi_gam_two_grp.R")

# The number of taxa, library size, and sample size
n.taxa=1000; samp.frac.var="moderate"; n.samp=c("20_30", "50_50")

# The proportion of differentially abundant taxa
prop.diff=c(0.05, 0.25, 0.50, 0.75)

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

library(DESeq2)
library(doParallel)
library(foreach)

start_time <- Sys.time()
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
  test.dat=abn.tab.gen(n.taxa, n.samp.grp1, n.samp.grp2, low.abn, med.abn, high.abn,
                       prop.diff, abn.seed, obs.seed, struc.zero.prop, out.zero.prop,
                       samp.frac.var)
  
  # Prepare data for DESeq2
  countdata=test.dat$obs.abn # Format for DESeq2: taxa are rows
  coldata=data.frame(group=as.factor(rep(c(1, 2), c(n.samp.grp1, n.samp.grp2))))
  rownames(coldata)=paste0("sub", seq(n.samp.grp1+n.samp.grp2))
  
  zero.threshold=0.90
  taxa.info.ind=apply(countdata, 1, function(x) sum(x==0)/(n.samp.grp1+n.samp.grp2))
  feature_table=round(countdata[which(taxa.info.ind<zero.threshold), ])+1L
  
  count.table=DESeqDataSetFromMatrix(
    countData = feature_table, colData = coldata, design = ~ group)
  
  # Run DESeq2
  suppressWarnings(dds <- try(DESeq(count.table, quiet = TRUE), silent = TRUE))
  if (inherits(dds, "try-error")) {
    # If the parametric fit failed, try the local.
    suppressWarnings(dds <- try(DESeq(count.table, fitType = "local", quiet = TRUE), silent = TRUE))
    if (inherits(dds, "try-error")) {
      # If local fails, try the mean
      suppressWarnings(dds <- try(DESeq(count.table, fitType = "mean", quiet = TRUE), silent = TRUE))
    }
  }
  if (inherits(dds, "try-error")) {
    # If still bad
    FDR=NA; power=NA
  }else{
    res = results(dds)
    res$id = rownames(res)
    # Some DESeq2 results (for example) had NA adjusted p-values, so replace them with 1
    res[is.na(res[, "padj"]), "padj"] = 1
    
    res=data.frame(taxa=res$id,
                   diff.test=ifelse(res$padj<0.05, 1, 0), 
                   diff.ind=test.dat$diff.taxa[which(taxa.info.ind<zero.threshold)],
                   effec.size=test.dat$effect.size[which(taxa.info.ind<zero.threshold)])
    FDR=ifelse(sum(res$diff.test==1, na.rm = T)==0, 0,
               sum(ifelse(res$diff.ind==0&res$diff.test==1, 1, 0), na.rm = T)/
                 sum(res$diff.test==1, na.rm = T))
    power=sum(ifelse(res$diff.ind!=0&res$diff.test==1, 1, 0), na.rm = T)/
      sum(res$diff.ind!=0, na.rm = T)
    
  }
  
  c(FDR, power)
}
end_time <- Sys.time()
end_time - start_time

write_csv(data.frame(simlist), "fdr_power_deseq2_large.csv")
