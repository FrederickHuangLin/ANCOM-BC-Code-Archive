# Load dependencies
library(nloptr)
library(dplyr)

# Data Pre-Processing
feature_table_pre_process = function(feature.table, meta.data, sample.var, group.var, 
                                     zero.cut = 0.90, lib.cut = 1000, neg.lb){
  feature.table = data.frame(feature.table, check.names = FALSE)
  meta.data = data.frame(meta.data, check.names = FALSE)
  # Drop unused levels
  meta.data[] = lapply(meta.data, function(x) if(is.factor(x)) factor(x) else x)
  
  sample.ID = colnames(feature.table)
  meta.data = meta.data[match(sample.ID, meta.data[, sample.var]), ]
  
  # 1. Identify outliers within each taxon
  group = factor(meta.data[, group.var])
  group.name = levels(group)
  grp.ind.origin = lapply(1:nlevels(group), function(i) which(group == group.name[i]))
  n.grp.origin = length(grp.ind.origin)
  n.samp.grp.origin = sapply(grp.ind.origin, length)
  feature.table = feature.table[, unlist(grp.ind.origin)]
  
  z = log(feature.table + 1)
  f = z; f[f == 0] = NA; f = colMeans(f, na.rm = T)
  f.mean = unlist(tapply(f, rep(1:n.grp.origin, n.samp.grp.origin), mean))
  e = f-rep(f.mean, n.samp.grp.origin)
  y = t(t(z) - e)
  
  outlier_check = function(x){
    mu1 = quantile(x, 0.25); mu2 = quantile(x, 0.75)
    sigma1 = quantile(x, 0.75)-quantile(x, 0.25); sigma2 = sigma1
    pi = 0.75
    n = length(x)
    epsilon = 100
    tol = 1e-5
    score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1-pi)*dnorm(x, mean = mu2, sd = sigma2))
    while (epsilon > tol) {
      grp1.ind = score >= 1
      mu1.new = mean(x[grp1.ind]); mu2.new = mean(x[!grp1.ind])
      sigma1.new = sd(x[grp1.ind]); if(is.na(sigma1.new)) sigma1.new = 0
      sigma2.new = sd(x[!grp1.ind]); if(is.na(sigma2.new)) sigma2.new = 0
      pi.new = sum(grp1.ind)/n
      
      para = c(mu1.new, mu2.new, sigma1.new, sigma2.new, pi.new)
      if(any(is.na(para))) break
      
      score = pi.new*dnorm(x, mean = mu1.new, sd = sigma1.new)/
        ((1-pi.new)*dnorm(x, mean = mu2.new, sd = sigma2.new))
      
      epsilon = sqrt((mu1-mu1.new)^2 + (mu2-mu2.new)^2 + 
                       (sigma1-sigma1.new)^2 + (sigma2-sigma2.new)^2 + (pi-pi.new)^2)
      mu1 = mu1.new; mu2 = mu2.new; sigma1 = sigma1.new; sigma2 = sigma2.new; pi = pi.new
    }
    
    if(mu1 + 1.96*sigma1 < mu2 - 1.96*sigma2){
      if(pi > 0.85){
        out.ind = (!grp1.ind)
      }else if(pi < 0.15){
        out.ind = grp1.ind
      }else{
        out.ind = rep(FALSE, n)
      }
    }else{
      out.ind = rep(FALSE, n)
    }
    return(out.ind)
  }
  feature.table.out = t(apply(y, 1, function(i)
    unlist(tapply(i, rep(1:n.grp.origin, n.samp.grp.origin), function(j) outlier_check(j)))))
  feature.table[feature.table.out] = NA
  
  # 2. Discard taxa with zeros  >=  zero.cut
  taxa.zero.prop = apply(feature.table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  filter.taxa = which(taxa.zero.prop >= zero.cut)
  if(length(filter.taxa)>0){
    feature.table = feature.table[-filter.taxa, ]
  }
  
  # 3. Discard samples with library size < lib.cut
  library.size = colSums(feature.table, na.rm = T)
  sample.ID = colnames(feature.table)
  meta.data = meta.data[match(sample.ID, meta.data[, sample.var]), ]
  
  if(any(library.size<lib.cut)){
    filter.subject = which(library.size<lib.cut)
    feature.table = feature.table[, -filter.subject]
    meta.data = meta.data[-filter.subject, ]
  }
  
  # 4. Re-order the OTU table
  group = factor(meta.data[, group.var])
  group.name = levels(group)
  grp.ind = lapply(1:nlevels(group), function(i) which(group == group.name[i]))
  n.grp = length(grp.ind)
  n.samp.grp = sapply(grp.ind, length)
  
  n.taxa = nrow(feature.table)
  taxa.id = rownames(feature.table)
  n.samp = ncol(feature.table)
  
  # 5. Identify taxa with structure zeros
  present.table = as.matrix(feature.table)
  present.table[is.na(present.table)] = 0
  present.table[present.table != 0] = 1
  
  p.hat.mat = t(apply(present.table, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) mean(y, na.rm = T)))))
  sample.size = t(apply(feature.table, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) length(y[!is.na(y)])))))
  p.hat.lo.mat = p.hat.mat - 1.96 * sqrt(p.hat.mat*(1 - p.hat.mat)/sample.size)
  colnames(p.hat.mat) = levels(group)
  colnames(p.hat.lo.mat) = levels(group)
  
  struc.zero = matrix(0, nrow = n.taxa, ncol = n.grp)
  struc.zero[p.hat.mat == 0] = 1
  # Whether we need to classify a taxon into structural zero by its negative lower bound?
  if(neg.lb) struc.zero[p.hat.lo.mat <= 0] = 1
  rownames(struc.zero) = taxa.id
  colnames(struc.zero) = paste0("structural.zero (", levels(group), ")")
  
  # Entries considered to be structural zeros are set to be 0s
  ind.zero = struc.zero[, rep(1:n.grp, times = n.samp.grp)]
  feature.table = feature.table * (1 - ind.zero)
  
  # 6. Return results
  res = list(feature.table = feature.table, library.size = library.size, 
             group.name = group.name, group.ind = grp.ind, structure.zeros = struc.zero)
  return(res)
}

# ANCOM-BC main function
ANCOM_BC = function(feature.table, grp.name, grp.ind, struc.zero, adj.method = "bonferroni", 
                    tol.EM = 1e-5, max.iterNum = 100, perNum = 1000, alpha = 0.05){
  n.taxa.raw = nrow(feature.table)
  taxa.id.raw = rownames(feature.table)
  n.samp = ncol(feature.table)
  sample.id = colnames(feature.table)
  
  n.grp = length(grp.ind)
  n.samp.grp = sapply(grp.ind, length)
  
  ### 0. Discard taxa with structural zeros for the moment
  comp.taxa.pos = which(apply(struc.zero, 1, function(x) all(x == 0))) # position of complete taxa (no structural zeros)
  O = feature.table[comp.taxa.pos, ]
  n.taxa = nrow(O)
  taxa.id = rownames(O)
  n.samp = ncol(O)
  y = as.matrix(log(O + 1))
  
  ### 1. Initial estimates of sampling fractions and mean absulute abundances
  mu = t(apply(y, 1, function(i) tapply(i, rep(1:n.grp, n.samp.grp), function(j)
    mean(j, na.rm = T))))
  d = colMeans(y - mu[, rep(1:n.grp, times = n.samp.grp)], na.rm = T)
  
  ## Iteration in case of missing values of y
  iterNum = 0
  epsilon = 100
  while (epsilon > tol.EM & iterNum < max.iterNum) {
    # Updating mu
    mu.new = t(apply(t(t(y) - d), 1, function(i) tapply(i, rep(1:n.grp, n.samp.grp), function(j)
      mean(j, na.rm = T))))
    
    # Updating d
    d.new = colMeans(y - mu.new[, rep(1:ncol(mu.new), times = n.samp.grp)], na.rm = T)
    
    # Iteration
    epsilon = sqrt(sum((mu.new - mu)^2) + sum((d.new - d)^2))
    iterNum = iterNum + 1
    
    mu = mu.new
    d = d.new
  }
  
  mu.var.each = (y-t(t(mu[, rep(1:ncol(mu), times = n.samp.grp)])+d))^2
  mu.var = t(apply(mu.var.each, 1, function(x) tapply(x, rep(1:n.grp, n.samp.grp), function(y)
    mean(y, na.rm = T))))
  sample.size = t(apply(y, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) length(y[!is.na(y)])))))
  mu.var = mu.var/sample.size
  
  ### 2. Estimate the bias (between-group difference of sampling fractions) by E-M algorithm
  bias.em.vec = rep(NA, n.grp - 1)
  bias.wls.vec = rep(NA, n.grp - 1)
  bias.var.vec = rep(NA, n.grp - 1)
  for (i in 1:(n.grp-1)) {
    Delta = mu[, 1] - mu[, 1+i]
    nu = rowSums(mu.var[, c(1, 1+i)])
    
    ## 2.1 Initials
    pi0_0 = 0.75
    pi1_0 = 0.125
    pi2_0 = 0.125
    delta_0 = mean(Delta[Delta >= quantile(Delta, 0.25, na.rm = T)&
                         Delta <= quantile(Delta, 0.75, na.rm = T)], na.rm = T)
    l1_0 = mean(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    l2_0 = mean(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    kappa1_0 = var(Delta[Delta < quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
    kappa2_0 = var(Delta[Delta > quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    if(is.na(kappa2_0)|kappa2_0 == 0) kappa2_0 = 1
    
    ## 2.2 Apply E-M algorithm
    # 2.21 Store all paras in vectors/matrices
    pi0.vec = c(pi0_0); pi1.vec = c(pi1_0); pi2.vec = c(pi2_0)
    delta.vec = c(delta_0); l1.vec = c(l1_0); l2.vec = c(l2_0)
    kappa1.vec = c(kappa1_0); kappa2.vec = c(kappa2_0)
    
    # 2.22 E-M iteration
    iterNum = 0
    epsilon = 100
    while (epsilon > tol.EM & iterNum < max.iterNum) {
      # print(iterNum)
      ## Current value of paras
      pi0 = pi0.vec[length(pi0.vec)]; pi1 = pi1.vec[length(pi1.vec)]; pi2 = pi2.vec[length(pi2.vec)]
      delta = delta.vec[length(delta.vec)]; 
      l1 = l1.vec[length(l1.vec)]; l2 = l2.vec[length(l2.vec)]
      kappa1 = kappa1.vec[length(kappa1.vec)]; kappa2 = kappa2.vec[length(kappa2.vec)]
      
      ## E-step
      pdf0 = sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta, sqrt(nu[i])))
      pdf1 = sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta + l1, sqrt(nu[i] + kappa1)))
      pdf2 = sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta + l2, sqrt(nu[i] + kappa2)))
      r0i = pi0*pdf0/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r0i[is.na(r0i)] = 0
      r1i = pi1*pdf1/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r1i[is.na(r1i)] = 0
      r2i = pi2*pdf2/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2); r2i[is.na(r2i)] = 0
      
      ## M-step
      pi0_new = mean(r0i, na.rm = T); pi1_new = mean(r1i, na.rm = T); pi2_new = mean(r2i, na.rm = T)
      delta_new = sum(r0i*Delta/nu + r1i*(Delta-l1)/(nu+kappa1) + r2i*(Delta-l2)/(nu+kappa2), na.rm = T)/
        sum(r0i/nu + r1i/(nu+kappa1) + r2i/(nu+kappa2), na.rm = T)
      l1_new = min(sum(r1i*(Delta-delta)/(nu+kappa1), na.rm = T)/sum(r1i/(nu+kappa1), na.rm = T), 0)
      l2_new = max(sum(r2i*(Delta-delta)/(nu+kappa2), na.rm = T)/sum(r2i/(nu+kappa2), na.rm = T), 0)
      
      # Nelder-Mead simplex algorithm for kappa1 and kappa2
      obj.kappa1 = function(x){
        log.pdf = log(sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+l1, sqrt(nu[i]+x))))
        log.pdf[is.infinite(log.pdf)] = 0
        -sum(r1i*log.pdf, na.rm = T)
      }
      kappa1_new = neldermead(x0 = kappa1, fn = obj.kappa1, lower = 0)$par
      
      obj.kappa2 = function(x){
        log.pdf = log(sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+l2, sqrt(nu[i]+x))))
        log.pdf[is.infinite(log.pdf)] = 0
        -sum(r2i*log.pdf, na.rm = T)
      }
      kappa2_new = neldermead(x0 = kappa2, fn = obj.kappa2, lower = 0)$par
      
      ## Merge to the paras vectors/matrices
      pi0.vec = c(pi0.vec, pi0_new); pi1.vec = c(pi1.vec, pi1_new); pi2.vec = c(pi2.vec, pi2_new)
      delta.vec = c(delta.vec, delta_new)
      l1.vec = c(l1.vec, l1_new); l2.vec = c(l2.vec, l2_new)
      kappa1.vec = c(kappa1.vec, kappa1_new); kappa2.vec = c(kappa2.vec, kappa2_new)
      
      ## Calculate the new epsilon
      epsilon = sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 + (delta_new-delta)^2+
                       (l1_new-l1)^2 + (l2_new-l2)^2 + (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
      iterNum = iterNum+1
    }
    # 2.23 Estimate the bias
    bias.em.vec[i] = delta.vec[length(delta.vec)]
    
    # 2.24 The WLS estimator of bias
    # Cluster 0
    C0 = which(Delta >= quantile(Delta, pi1_new, na.rm = T) & Delta < quantile(Delta, 1 - pi2_new, na.rm = T))
    # Cluster 1
    C1 = which(Delta < quantile(Delta, pi1_new, na.rm = T))
    # Cluster 2
    C2 = which(Delta >= quantile(Delta, 1 - pi2_new, na.rm = T))
    
    nu_temp = nu
    nu_temp[C1] = nu_temp[C1] + kappa1_new
    nu_temp[C2] = nu_temp[C2] + kappa2_new
    wls.deno = sum(1 / nu_temp)
    
    wls.nume = 1 / nu_temp
    wls.nume[C0] = (wls.nume * Delta)[C0]
    wls.nume[C1] = (wls.nume * (Delta - l1_new))[C1]
    wls.nume[C2] = (wls.nume * (Delta - l2_new))[C2];   
    wls.nume = sum(wls.nume)
    
    bias.wls.vec[i] = wls.nume / wls.deno
      
    # 2.25 Estimate the variance of bias  
    bias.var.vec[i] = 1 / wls.deno
    if (is.na(bias.var.vec[i])) bias.var.vec[i] = 0
  }
  bias.em.vec = c(0, bias.em.vec)
  bias.wls.vec = c(0, bias.wls.vec)
  
  ### 3. Final estimates of mean absolute abundane and sampling fractions
  mu.adj.comp = t(t(mu) + bias.em.vec)
  colnames(mu.adj.comp) = grp.name; rownames(mu.adj.comp) = taxa.id
  
  d.adj = d - rep(bias.em.vec, sapply(grp.ind, length))
  names(d.adj) = sample.id
  
  ### 4. Hypothesis testing
  W.numerator = matrix(apply(mu.adj.comp, 1, function(x) combn(x, 2, FUN = diff)), ncol = n.taxa)
  W.numerator = t(W.numerator)
  # Variance of estimated mean difference
  W.denominator1 = matrix(apply(mu.var, 1, function(x) combn(x, 2, FUN = sum)), ncol = n.taxa)
  # Variance of delta_hat
  if (length(bias.var.vec) < 2) {
    W.denominator2 = bias.var.vec
  }else {
    W.denominator2 = c(bias.var.vec, combn(bias.var.vec, 2, FUN = sum))
  }
  W.denominator = W.denominator1 + W.denominator2 + 2 * sqrt(W.denominator1 * W.denominator2)
  W.denominator = t(sqrt(W.denominator))
  grp.pair = combn(n.grp, 2)
  colnames(W.numerator) = sapply(1:ncol(grp.pair), function(x) 
    paste0("mean.difference (", grp.name[grp.pair[2, x]], " - ", grp.name[grp.pair[1, x]], ")"))
  colnames(W.denominator) = sapply(1:ncol(grp.pair), function(x) 
    paste0("se (", grp.name[grp.pair[2, x]], " - ", grp.name[grp.pair[1, x]], ")"))
  rownames(W.numerator) = taxa.id; rownames(W.denominator) = taxa.id
  
  if (length(grp.name) == 2) {
    ## Two-group comparison
    W = W.numerator/W.denominator
    p.val = sapply(W, function(x) 2*pnorm(abs(x), mean = 0, sd = 1, lower.tail = F))
    q.val = p.adjust(p.val, method = adj.method)
    q.val[is.na(q.val)] = 1
  } else {
    ## Multi-group comparison: Permutation test
    # Test statistics
    W.each = W.numerator/W.denominator
    W.each[is.na(W.each)] = 0 # Replace missing values with 0s
    W = apply(abs(W.each), 1, max)
    
    # Test statistics under null
    W.null.list = lapply(1:perNum, function(x) {
      set.seed(x)
      mu.adj.comp.null = matrix(rnorm(n.taxa * n.grp), nrow = n.taxa, ncol = n.grp) * sqrt(mu.var)
      W.numerator.null = matrix(apply(mu.adj.comp.null, 1, function(x) combn(x, 2, FUN = diff)), ncol = n.taxa)
      W.numerator.null = t(W.numerator.null)
      
      W.each.null = W.numerator.null/W.denominator
      W.each.null[is.na(W.each.null)] = 0
      W.null = apply(abs(W.each.null), 1, max)
      return(W.null)
    })
    W.null = Reduce('cbind', W.null.list)
    
    # Test results
    p.val = apply(W.null - W, 1, function(x) sum(x > 0)/perNum)
    q.val = p.adjust(p.val, method = adj.method)
    q.val[is.na(q.val)] = 1
  }
  
  W = matrix(W, ncol = 1)
  colnames(W) = "W"
  res.comp = data.frame(W.numerator, W.denominator, W = W, p.val, q.val, check.names = FALSE)
  
  ### 5. Combine results from structural zeros
  mu.adj = matrix(NA, nrow = n.taxa.raw, ncol = n.grp)
  colnames(mu.adj) = grp.name; rownames(mu.adj) = taxa.id.raw
  mu.adj[comp.taxa.pos, ] = mu.adj.comp
  
  if (length(comp.taxa.pos) < n.taxa.raw) {
    O.incomp = feature.table[-comp.taxa.pos, ]
    ind.incomp = struc.zero[-comp.taxa.pos, rep(1:n.grp, times = n.samp.grp)]
    y.incomp = log(O.incomp + 1)
    d.incomp = t(t(1 - ind.incomp) * d) # Sampling fractions for entries considered to be structural zeros are set to be 0s
    y.adj.incomp = y.incomp - d.incomp
    mu.incomp = t(apply(y.adj.incomp, 1, function(i) 
      tapply(i, rep(1:n.grp, n.samp.grp), function(j) mean(j, na.rm = T))))
    # In case of negative values for mean absolute abundances
    mu.adj.incomp = mu.incomp
    mu.adj.incomp[mu.adj.incomp == 0] = NA
    mu.adj.incomp = t(t(mu.adj.incomp) + abs(apply(mu.incomp, 2, min)))
    mu.adj.incomp[is.na(mu.adj.incomp)] = 0
  }else{
    mu.adj.incomp = NA
  }
  mu.adj[-comp.taxa.pos, ] = mu.adj.incomp
  colnames(mu.adj) = paste0("mean.absolute.abundance (", grp.name, ")")
  rownames(mu.adj) = taxa.id.raw
  
  ### 6. Outputs
  W.numerator = matrix(apply(mu.adj, 1, function(x) combn(x, 2, FUN = diff)), ncol = n.taxa.raw)
  W.numerator = t(W.numerator)
  W.denominator = matrix(0, ncol = ncol(W.numerator), nrow = nrow(W.numerator))
  
  res = data.frame(W.numerator, W.denominator, W = Inf, p.val = 0, q.val = 0, check.names = FALSE)
  res[comp.taxa.pos, ] = res.comp
  colnames(res) = colnames(res.comp); rownames(res) = taxa.id.raw
  res = res%>%mutate(diff.abn = ifelse(q.val < alpha, TRUE, FALSE))
  
  out = list(feature.table = feature.table, res = res, d = d.adj, mu = mu.adj, bias.em = bias.em.vec, bias.wls = bias.wls.vec)
  return(out)
}
