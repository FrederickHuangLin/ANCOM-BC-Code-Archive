abn.tab.gen2=function(n.taxa, n.samp.grp, low.abn, med.abn, high.abn,
                      prop.diff, abn.seed, obs.seed, struc.zero.prop, out.zero.prop){
  # Total number of samples
  n.samp=sum(n.samp.grp)
  n.grp=length(n.samp.grp)
  
  set.seed(abn.seed) # This seed is used to control whether you would like to have the same population
  low.prop=0.6 # Proportion of low abundance 
  med.prop=0.3 # Proportion of medium abundance
  hi.prop=0.1  # Proportion of high abundance
  # Indices for taxa abundance 
  index=sample(c(1, 2, 3), n.taxa, replace = T, prob = c(low.prop, med.prop, hi.prop)) 
  
  # Poisson parameters
  lambda=rep(NA, n.taxa)
  lambda[which(index==1)]=rgamma(length(which(index==1)), shape=low.abn, rate=1)
  lambda[which(index==2)]=rgamma(length(which(index==2)), shape=med.abn, rate=1)
  lambda[which(index==3)]=rgamma(length(which(index==3)), shape=high.abn, rate=1)
  
  # Which taxa are differentially abundant
  diff.ind=rep(0, n.taxa) # 0 = nondifferentially abundant
  diff.taxa=sample(c(1:n.taxa), floor(n.taxa*prop.diff), replace=F)
  # k = # groups that are differentially abundant with the 1st group
  diff.ind[diff.taxa]=sapply(diff.ind[diff.taxa], function(x) x=sample(1:(n.grp-1), size = 1))
  # Structural zeros
  struc.zero.ind=sample(which(diff.ind==0), struc.zero.prop*length(which(diff.ind==0)), replace = F)
  # |k| = # groups that are differentially abundant with the 1st group
  diff.ind[struc.zero.ind]=sapply(diff.ind[struc.zero.ind], function(x) 
    x=-sample(1:(n.grp-1), size = 1))
  
  # Signal size
  effect.size=matrix(1, nrow = n.taxa, ncol = n.grp)
  for (i in diff.taxa) {
    k = unlist(diff.ind[i])
    effect.size[i, ]=c(1,
                       sample(c(runif(k, 1, 10), runif(k, 0.1, 1)), k),
                       rep(1, n.grp-k-1))
  }
  for (i in struc.zero.ind) {
    k = abs(diff.ind[i])
    effect.size[i, ]=c(1,
                       rep(0, k),
                       rep(1, n.grp-k-1))
  }
  effect.size=apply(effect.size, 2, as.numeric)
  rownames(effect.size)=paste0("taxon", seq(n.taxa))
  colnames(effect.size)=paste0("group", seq(n.grp))
  
  # Mean absolute abundance in the ecosystem
  temp.mat=apply(effect.size, 2, function(x) x*lambda)
  
  # Absolute abundance in the ecosystem
  abn.list=vector(mode = "list", length = n.grp)
  for (i in 1:n.grp) {
    abn.list[[i]]=rpois(n.taxa*n.samp.grp[i], rep(temp.mat[, i], n.samp.grp[i]))
    abn.list[[i]]=matrix(abn.list[[i]], ncol=n.samp.grp[i])
  }
  abn.mat=Reduce('cbind', abn.list)
  # Outlier zeros
  out.ind=rep(0, n.taxa); out.ind[sample(seq(n.taxa), out.zero.prop*n.taxa, replace = F)]=1
  names(out.ind)=paste0("taxon", seq(n.taxa))
  abn.mat[which(out.ind==1), sample(seq(n.samp), out.zero.prop*n.samp, replace = F)]=0
  
  # Microbial load
  abn.total=colSums(abn.mat)
  names(abn.total)=paste0("sub", seq(n.samp))
  
  # Library size
  depth=1/sample(c(runif(n.samp, 10, 50), runif(n.samp, 100, 500)), n.samp, replace = T)
  obs.total=round(max(abn.total)*depth)
  names(obs.total)=paste0("sub", seq(n.samp))
  
  # Absolute abundance in the sample
  set.seed(obs.seed)
  obs.list=lapply(1:n.samp, function(i) 
    phyloseq:::rarefaction_subsample(x=abn.mat[, i], sample.size=obs.total[i]))
  obs.mat=Reduce('cbind', obs.list)
  
  # Prepare outputs
  temp.dat=data.frame(temp.mat, row.names = NULL)
  colnames(temp.dat)=paste0("group", seq(n.grp))
  rownames(temp.dat)=paste0("taxon", seq(n.taxa))
  
  abn.dat=data.frame(abn.mat, row.names = NULL)
  colnames(abn.dat)=paste0("sub", seq(n.samp))
  rownames(abn.dat)=paste0("taxon", seq(n.taxa))
  
  obs.dat=data.frame(obs.mat, row.names = NULL)
  colnames(obs.dat)=paste0("sub", seq(n.samp))
  rownames(obs.dat)=paste0("taxon", seq(n.taxa))
  
  diff.ind=unlist(diff.ind)
  names(diff.ind)=paste0("taxon", seq(n.taxa))
  
  c=obs.total/abn.total
  names(c)=paste0("sub", seq(n.samp))
  
  test.data=list(temp.dat, abn.dat, obs.dat, effect.size, 
                 diff.ind, out.ind, c, abn.total, obs.total)
  names(test.data)=c("template", "pop.abn", "obs.abn", "effect.size", "diff.taxa", 
                     "outlier.ind", "samp.depth", "abn.total", "obs.total")
  return(test.data)
}