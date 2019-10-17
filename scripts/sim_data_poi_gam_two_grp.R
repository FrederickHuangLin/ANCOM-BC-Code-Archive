library(phyloseq)

abn.tab.gen1=function(n.taxa, n.samp.grp1, n.samp.grp2, low.abn, med.abn, high.abn,
                      prop.diff, abn.seed, obs.seed, struc.zero.prop, out.zero.prop,
                      balanced.micro.load = c(TRUE, FALSE),
                      balanced.lib.size = c(TRUE, FALSE),
                      samp.frac = c("large", "small")){
  # Total number of samples
  n.samp=n.samp.grp1+n.samp.grp2
  
  set.seed(abn.seed)
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
  
  if(balanced.micro.load==TRUE){
    # Construct balanced microbial load in the ecosystem
    
    # Which taxa are differentially abundant
    diff.ind=rep(0, n.taxa)
    # Differentially abundant taxa
    diff.pos=sample(c(1:n.taxa), floor(n.taxa*prop.diff), replace=FALSE)
    diff.ind[diff.pos]=1
    # Structural zeros
    szero.pos=sample(which(diff.ind==0), struc.zero.prop*length(which(diff.ind==0)), replace = FALSE)
    diff.ind[szero.pos]=-1
    
    # Effect size
    effect.size=rep(1, n.taxa)
    effect.size[diff.pos]=runif(length(diff.pos), 1, 10)
    effect.size[szero.pos]=0
    names(effect.size)=paste0("taxon", seq(n.taxa))
    
    # Mean absolute abundance in the ecosystem
    temp.grp1=round(lambda); temp.grp2=round(lambda)
    for (i in c(diff.pos, szero.pos)) {
      rand.ind=sample(c(1, 2), 1)
      if (rand.ind==1){
        temp.grp1[i]=temp.grp1[i]*effect.size[i]
      }else{
        temp.grp2[i]=temp.grp2[i]*effect.size[i]
      }
    }
    temp.dat=data.frame(temp.grp1, temp.grp2, effect.size)
    rownames(temp.dat)=paste0("taxon", seq(n.taxa))
  }else{
    # Construct unbalanced microbial in the ecosystem
    
    # Which taxa are differentially abundant
    diff.ind=rep(0, n.taxa)
    # Group1 is higher than group2
    diff1.ind=sample(c(1:n.taxa), floor(n.taxa*prop.diff), replace=FALSE)
    diff.ind[diff1.ind]=1
    # Group2 is higher than group1
    wt=runif(1, 0, 1)
    diff2.ind=sample(diff1.ind, wt*length(diff1.ind), replace=FALSE)
    diff.ind[diff2.ind]=2
    # Structural zeros
    diff3.ind=sample(which(diff.ind==0), struc.zero.prop*length(which(diff.ind==0)), replace = FALSE)
    diff.ind[diff3.ind]=-1
    
    # Effect size
    effect.size=rep(1, n.taxa)
    effect.size[diff1.ind]=runif(length(diff1.ind), 1, 10)
    effect.size[diff2.ind]=runif(length(diff2.ind), 0.1, 1)
    effect.size[diff3.ind]=0
    names(effect.size)=paste0("taxon", seq(n.taxa))
    
    # Mean absolute abundance in the ecosystem
    temp.grp1=round(lambda*effect.size)
    temp.grp2=round(lambda)
    for (i in which(effect.size!=1)) {
      if(temp.grp1[i]==temp.grp2[i]) temp.grp1[i]=temp.grp1[i]+1
    }
    temp.dat=data.frame(temp.grp1, temp.grp2, effect.size)
    rownames(temp.dat)=paste0("taxon", seq(n.taxa))
  }
  
  # Absolute abundance in ecosystems
  abn.mat=matrix(0, ncol=n.samp, nrow=n.taxa)
  for(i in 1:n.taxa){
    abn.mat[i, ]=c(rpois(n.samp.grp1, temp.grp1[i]), rpois(n.samp.grp2, temp.grp2[i]))
  }
  # Outlier zeros
  out.ind=rep(0, n.taxa); out.ind[sample(seq(n.taxa), out.zero.prop*n.taxa, replace = F)]=1
  names(out.ind)=paste0("taxon", seq(n.taxa))
  abn.mat[which(out.ind==1), sample(seq(n.samp), out.zero.prop*n.samp, replace = F)]=0
  
  # Microbial load
  abn.total=colSums(abn.mat)
  names(abn.total)=paste0("sub", seq(n.samp))
  abn.prob.mat=t(t(abn.mat)/abn.total)
  
  if(samp.frac=="large"){
    depth=1/runif(n.samp, 5, 10)
  }else{
    depth=1/sample(c(runif(n.samp, 10, 50), runif(n.samp, 100, 500)), n.samp, replace = T)
  }
  
  # library size
  if(balanced.lib.size==TRUE){
    # Construct balanced library sizes
    obs.total=round(max(abn.total)*depth)
  }else{
    # Construct unbalanced library sizes
    obs.total=round(abn.total*depth)
  }
  names(obs.total)=paste0("sub", seq(n.samp))
  
  # Absolute abundance in samples
  set.seed(obs.seed)
  obs.list=lapply(1:n.samp, function(i)
    phyloseq:::rarefaction_subsample(x=abn.mat[, i], sample.size=obs.total[i]))
  obs.mat=Reduce('cbind', obs.list)
  
  # Prepare outputs
  abn.dat=data.frame(abn.mat, row.names = NULL)
  rownames(abn.dat)=paste0("taxon", seq(n.taxa))
  colnames(abn.dat)=paste0("sub", seq(n.samp))
  
  obs.dat=data.frame(obs.mat, row.names = NULL)
  rownames(obs.dat)=paste0("taxon", seq(n.taxa))
  colnames(obs.dat)=paste0("sub", seq(n.samp))
  
  abn.prob.dat=data.frame(abn.prob.mat, row.names = NULL)
  rownames(abn.prob.dat)=paste0("taxon", seq(n.taxa))
  colnames(abn.prob.dat)=paste0("sub", seq(n.samp))
  
  grp.ind=c(rep(1, n.samp.grp1), rep(2, n.samp.grp2))
  names(grp.ind)=paste0("sub", seq(n.samp))
  
  names(diff.ind)=paste0("taxon", seq(n.taxa))
  
  # Sampling fractions
  c.mult=obs.total/abn.total
  names(c.mult)=paste0("sub", seq(n.samp))
  
  test.data=list(temp.dat, abn.dat, obs.dat, effect.size, grp.ind, 
                 diff.ind, out.ind, c.mult, abn.total, obs.total)
  names(test.data)=c("mean.eco.abn", "eco.abn", "obs.abn", "effect.size", "grp", 
                     "diff.taxa", "outlier", "samp.frac", "micro.load", "lib.size")
  return(test.data)
}