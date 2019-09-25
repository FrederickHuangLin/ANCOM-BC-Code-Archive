library(phyloseq)

# Abudance table generating function
abn.tab.gen = function(template, n.taxa, n.samp.grp1, n.samp.grp2, prop.diff, abn.seed, obs.seed,
                       struc.zero.prop){
  # Template for absolute abundance in the ecosystem
  set.seed(abn.seed)
  abn.temp = sample(template, n.taxa)
  taxa.id = names(abn.temp); n.samp = n.samp.grp1 + n.samp.grp2

  # Set which taxa are differentially abundant
  diff.ind = rep(0, n.taxa)
  # Group1 is higher than group2
  diff1.ind = sample(1:n.taxa, floor(n.taxa*prop.diff), replace = FALSE)
  diff.ind[diff1.ind] = 1
  # Group2 is higher than group1
  wt = runif(1, 0, 1)
  diff2.ind = sample(diff1.ind, wt*length(diff1.ind), replace = FALSE)
  diff.ind[diff2.ind] = 2
  # Structural zeros
  diff3.ind = sample(which(diff.ind == 0), struc.zero.prop * length(which(diff.ind==0)), replace = FALSE)
  diff.ind[diff3.ind] = -1
  
  # effect size
  effect.size = rep(1, n.taxa)
  effect.size[diff1.ind] = runif(length(diff1.ind), 1, 10)
  effect.size[diff2.ind] = runif(length(diff2.ind), 0.1, 1)
  effect.size[diff3.ind] = 0
  names(effect.size) = taxa.id
  
  # Absolute abundance in the ecosystem of two groups
  abn.grp1 = round(abn.temp * effect.size)
  abn.grp2 = round(abn.temp)
  abn.mat = cbind(abn.grp1, abn.grp2)
  
  abn.total = colSums(abn.mat)
  names(abn.total) = c("grp1", "grp2")
  
  # Number of taxa that are sampled for each subject
  depth = 1/sample(c(runif(n.samp, 10, 50), runif(n.samp, 100, 500)), n.samp, replace  =  T)
  obs.total = round(max(abn.total)*depth)
  names(obs.total) = paste0("sub", seq(n.samp))
  
  # Specimen abundance
  set.seed(obs.seed)
  obs.mat = matrix(NA, nrow = n.taxa, ncol  =  n.samp)
  for (i in 1:n.samp.grp1) {
    obs.mat[, i] = phyloseq:::rarefaction_subsample(x = abn.mat[, 1], sample.size = obs.total[i])
  }
  for (i in (n.samp.grp1 + 1):n.samp) {
    obs.mat[, i] = phyloseq:::rarefaction_subsample(x = abn.mat[, 2], sample.size = obs.total[i])
  }
  
  # Prepare output data sets
  abn.dat = data.frame(abn.mat, row.names  =  NULL)
  rownames(abn.dat) = taxa.id
  colnames(abn.dat) = c("grp1", "grp2")
  
  obs.dat = data.frame(obs.mat, row.names  =  NULL)
  rownames(obs.dat) = taxa.id
  colnames(obs.dat) = paste0("sub", seq(n.samp))
  
  grp.ind = c(rep(1, n.samp.grp1), rep(2, n.samp.grp2))
  names(grp.ind) = paste0("sub", seq(n.samp))
  
  names(diff.ind) = taxa.id
  
  c.mult = c(obs.total[1:n.samp.grp1]/abn.total[1], 
             obs.total[(n.samp.grp1 + 1):n.samp]/abn.total[2])
  names(c.mult) = paste0("sub", seq(n.samp))
  
  test.data = list(abn.dat, obs.dat, effect.size, grp.ind, 
                   diff.ind, c.mult, abn.total, obs.total)
  names(test.data) = c("pop.abn", "obs.abn", "effect.size", "grp", 
                       "diff.taxa", "mult", "abn.total", "obs.total")
  return(test.data)
}
