library(tidyverse)

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

## Read in original data
dat.ancom_bc=read_csv("archive/check/fdr_power_ancom_bc_large.csv")
dat.deseq2=read_csv("archive/check/fdr_power_deseq2_large.csv")
dat.edger=read_csv("archive/check/fdr_power_edger_large.csv")
dat.zilg=read_csv("archive/check/fdr_power_zilg_large.csv")
dat.zig=read_csv("archive/check/fdr_power_zig_large.csv")
dat.wilcox_un=read_csv("archive/check/fdr_power_wilcox_un_large.csv")
dat.wilcox_tss=read_csv("archive/check/fdr_power_wilcox_tss_large.csv")

## Reshaping data
simpattern=distinct(simparams, n.taxa, n.samp.grp1, n.samp.grp2, prop.diff)

data_summary = function(eval_data, method){
  FDR=tapply(as.numeric(eval_data[1, ]), 
             rep(seq(nrow(simpattern)), each=iterNum), function(x) mean(x, na.rm = T))
  FDRSD=tapply(as.numeric(eval_data[1, ]), 
               rep(seq(nrow(simpattern)), each=iterNum), function(x) sd(x, na.rm = T))
  power=tapply(as.numeric(eval_data[2, ]), 
               rep(seq(nrow(simpattern)), each=iterNum), function(x) mean(x, na.rm = T))
  powerSD=tapply(as.numeric(eval_data[2, ]), 
                 rep(seq(nrow(simpattern)), each=iterNum), function(x) sd(x, na.rm = T))
  data_sum = data.frame(FDR, FDRSD, power, powerSD, simpattern, method)
  data_sum = data_sum%>%unite(n.samp.grp, n.samp.grp1, n.samp.grp2, sep = ", ")
  
  return(data_sum)
}

eval.dat.list = list(dat.ancom_bc, dat.deseq2, dat.edger, 
                     dat.zilg, dat.zig, dat.wilcox_un, dat.wilcox_tss)
method.list = list("ANCOM-BC", "DESeq2", "edgeR", "ZILG", "ZIG", "Wilcoxon", "Wilcoxon + TSS")

dat.fig.list = vector(mode = "list", length = length(eval.dat.list))
for (i in 1:length(eval.dat.list)) {
  dat.fig.list[[i]] = data_summary(eval.dat.list[[i]], method.list[[i]])
}

## Merge data
dat.fig=Reduce('rbind', dat.fig.list)
dat.fig$n.samp.grp=factor(dat.fig$n.samp.grp)
dat.fig$method=factor(dat.fig$method)
dat.fig$prop.diff=factor(dat.fig$prop.diff)



