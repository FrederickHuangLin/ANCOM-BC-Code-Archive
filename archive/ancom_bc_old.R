# Load dependencies
library(nloptr)
library(dplyr)

# Data Pre-Processing
feature_table_pre_process=function(feature.table, meta.data, sample.var, group.var, 
                                   zero.cut, lib.cut, neg.lb){
  feature.table=data.frame(feature.table, check.names=FALSE)
  meta.data=data.frame(meta.data, check.names=FALSE)
  # Drop unused levels
  meta.data[]=lapply(meta.data, function(x) if(is.factor(x)) factor(x) else x)
  
  sample.ID=colnames(feature.table)
  meta.data=meta.data[match(sample.ID, meta.data[, sample.var]), ]
  
  # 1. Identify outliers within each taxon
  group=factor(meta.data[, group.var])
  group.name=levels(group)
  grp.ind.origin=lapply(1:nlevels(group), function(i) which(group==group.name[i]))
  n.grp.origin=length(grp.ind.origin)
  n.samp.grp.origin=sapply(grp.ind.origin, length)
  feature.table=feature.table[, unlist(grp.ind.origin)]
  
  z=log(feature.table+1)
  f=z; f[f==0]=NA; f=colMeans(f, na.rm = T)
  f.mean=unlist(tapply(f, rep(1:n.grp.origin, n.samp.grp.origin), mean))
  e=f-rep(f.mean, n.samp.grp.origin)
  y=t(t(z)-e)
  
  outlier_check=function(x){
    mu1=quantile(x, 0.25); mu2=quantile(x, 0.75)
    sigma1=quantile(x, 0.75)-quantile(x, 0.25); sigma2=sigma1
    pi=0.75
    n=length(x)
    epsilon=100
    tol=1e-5
    score=pi*dnorm(x, mean = mu1, sd=sigma1)/((1-pi)*dnorm(x, mean = mu2, sd=sigma2))
    while (epsilon>tol) {
      grp1.ind=score>=1
      mu1.new=mean(x[grp1.ind]); mu2.new=mean(x[!grp1.ind])
      sigma1.new=sd(x[grp1.ind]); if(is.na(sigma1.new)) sigma1.new=0
      sigma2.new=sd(x[!grp1.ind]); if(is.na(sigma2.new)) sigma2.new=0
      pi.new=sum(grp1.ind)/n
      
      para=c(mu1.new, mu2.new, sigma1.new, sigma2.new, pi.new)
      if(any(is.na(para))) break
      
      score=pi.new*dnorm(x, mean = mu1.new, sd=sigma1.new)/
        ((1-pi.new)*dnorm(x, mean = mu2.new, sd=sigma2.new))
      
      epsilon=sqrt((mu1-mu1.new)^2+(mu2-mu2.new)^2+
                     (sigma1-sigma1.new)^2+(sigma2-sigma2.new)^2+(pi-pi.new)^2)
      mu1=mu1.new; mu2=mu2.new; sigma1=sigma1.new; sigma2=sigma2.new; pi=pi.new
    }
    
    if(mu1+1.96*sigma1<mu2-1.96*sigma2){
      if(pi>0.85){
        out.ind=(!grp1.ind)
      }else if(pi<0.15){
        out.ind=grp1.ind
      }else{
        out.ind=rep(FALSE, n)
      }
    }else{
      out.ind=rep(FALSE, n)
    }
    return(out.ind)
  }
  feature.table.out=t(apply(y, 1, function(i)
    unlist(tapply(i, rep(1:n.grp.origin, n.samp.grp.origin), function(j) outlier_check(j)))))
  feature.table[feature.table.out]=NA
  
  # 2. Discard taxa with zeros >= zero.cut
  taxa.zero.prop=apply(feature.table, 1, function(x) sum(x==0, na.rm = T)/length(x[!is.na(x)]))
  filter.taxa=which(taxa.zero.prop>=zero.cut)
  if(length(filter.taxa)>0){
    feature.table=feature.table[-filter.taxa, ]
  }
  
  # 3. Discard samples with library size < lib.cut
  library.size=colSums(feature.table, na.rm = T)
  sample.ID=colnames(feature.table)
  meta.data=meta.data[match(sample.ID, meta.data[, sample.var]), ]
  
  if(any(library.size<lib.cut)){
    filter.subject=which(library.size<lib.cut)
    feature.table=feature.table[, -filter.subject]
    meta.data=meta.data[-filter.subject, ]
  }
  
  # 4. Re-order the OTU table
  group=factor(meta.data[, group.var])
  group.name=levels(group)
  grp.ind=lapply(1:nlevels(group), function(i) which(group==group.name[i]))
  n.grp=length(grp.ind)
  n.samp.grp=sapply(grp.ind, length)
  
  n.taxa=nrow(feature.table)
  taxa.id=rownames(feature.table)
  n.samp=ncol(feature.table)
  
  # 5. Identify taxa with structure zeros
  present.table=as.matrix(feature.table)
  present.table[is.na(present.table)]=0
  present.table[present.table!=0]=1
  
  p.hat.mat=t(apply(present.table, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) mean(y, na.rm = T)))))
  sample.size=t(apply(feature.table, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) length(y[!is.na(y)])))))
  p.hat.lo.mat=p.hat.mat-1.96*sqrt(p.hat.mat*(1-p.hat.mat)/sample.size)
  colnames(p.hat.mat)=levels(group)
  colnames(p.hat.lo.mat)=levels(group)
  
  struc.zero=matrix(0, nrow = n.taxa, ncol = n.grp)
  struc.zero[p.hat.mat==0]=1
  # Whether we need to classify a taxon into structural zero by its negative lower bound?
  if(neg.lb) struc.zero[p.hat.lo.mat<=0]=1
  rownames(struc.zero)=taxa.id
  colnames(struc.zero)=paste0("Structural Zero in ", levels(group))
  
  # 6. Return results
  res=list(feature.table=feature.table, library.size=library.size, 
           group.name=group.name, group.ind=grp.ind, structure.zeros=struc.zero)
  return(res)
}

# ANCOM-BC main function
ANCOM_BC=function(feature.table, grp.name, grp.ind, struc.zero, adj.method, 
                  tol.EM, max.iterNum, perNum, alpha){
  n.taxa.origin=nrow(feature.table)
  taxa.id.origin=rownames(feature.table)
  n.samp=ncol(feature.table)
  sample.id=colnames(feature.table)
  
  n.grp=length(grp.ind)
  n.samp.grp=sapply(grp.ind, length)
  
  ### 0. Separate out taxa with no structural zeros
  info.taxa.pos=which(apply(struc.zero, 1, function(x) all(x==0)))
  O=feature.table[info.taxa.pos, ]
  n.taxa=nrow(O)
  taxa.id=rownames(O)
  n.samp=ncol(O)
  y=log(O+1)
  
  ### 1. Estimate sampling fractions and mean abundances by iteratively least squares
  d=rep(0, n.samp)
  mu=t(apply(y, 1, function(i) tapply(i, rep(1:n.grp, n.samp.grp), function(j)
    mean(j, na.rm = T))))
  iterNum=0
  epsilon=100
  while (epsilon>tol.EM&iterNum<max.iterNum) {
    # Updating mu
    mu.new=t(apply(t(t(y)-d), 1, function(i) tapply(i, rep(1:n.grp, n.samp.grp), function(j)
      mean(j, na.rm = T))))
    
    # Updating d
    d.new=colMeans(y-mu.new[, rep(1:ncol(mu.new), times = n.samp.grp)], na.rm = T)
    
    # Iteration
    epsilon=sqrt(sum((mu.new-mu)^2)+sum((d.new-d)^2))
    iterNum=iterNum+1
    
    mu=mu.new
    d=d.new
  }
  
  mu.var.raw=(y-t(t(mu[, rep(1:ncol(mu), times = n.samp.grp)])+d))^2
  mu.var=t(apply(mu.var.raw, 1, function(x) tapply(x, rep(1:n.grp, n.samp.grp), function(y)
    mean(y, na.rm = T))))
  sample.size=t(apply(y, 1, function(x)
    unlist(tapply(x, rep(1:n.grp, n.samp.grp), function(y) length(y[!is.na(y)])))))
  mu.var=mu.var/sample.size
  
  ### 2. Estimate the bias of sampling fractions by E-M algorithm
  bias.est.vec=rep(NA, n.grp-1)
  for (i in 1:(n.grp-1)) {
    Delta=mu[, 1]-mu[, 1+i]
    Delta.var.est=rowSums(mu.var[, c(1, 1+i)])
    
    ## 2.1 Initials
    pi1_0=0.125
    pi2_0=0.75
    pi3_0=0.125
    delta_0=mean(Delta[Delta>=quantile(Delta, 0.25, na.rm = T)&
                         Delta<=quantile(Delta, 0.75, na.rm = T)], na.rm = T)
    d1_0=mean(Delta[Delta<quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    d2_0=mean(Delta[Delta>quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    psi1.sq_0=var(Delta[Delta<quantile(Delta, 0.125, na.rm = T)], na.rm = T)
    if(is.na(psi1.sq_0)|psi1.sq_0==0) psi1.sq_0=1
    psi2.sq_0=var(Delta[Delta>quantile(Delta, 0.875, na.rm = T)], na.rm = T)
    if(is.na(psi2.sq_0)|psi2.sq_0==0) psi2.sq_0=1
    
    ## 2.2 Apply E-M algorithm
    # 2.21 Store all paras in vectors/matrices
    pi1.vec=c(pi1_0); pi2.vec=c(pi2_0); pi3.vec=c(pi3_0)
    delta.vec=c(delta_0); d1.vec=c(d1_0); d2.vec=c(d2_0)
    psi1.sq.vec=c(psi1.sq_0); psi2.sq.vec=c(psi2.sq_0)
    
    # 2.22 E-M iteration
    iterNum=0
    epsilon=100
    sigmai.sq=Delta.var.est
    while (epsilon>tol.EM&iterNum<max.iterNum) {
      # print(iterNum)
      ## Current value of paras
      pi1=pi1.vec[length(pi1.vec)]; pi2=pi2.vec[length(pi2.vec)]; pi3=pi3.vec[length(pi3.vec)]
      delta=delta.vec[length(delta.vec)]; 
      d1=d1.vec[length(d1.vec)]; d2=d2.vec[length(d2.vec)]
      psi1.sq=psi1.sq.vec[length(psi1.sq.vec)]; psi2.sq=psi2.sq.vec[length(psi2.sq.vec)]
      
      ## E-step
      pdf1=sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+d1, sqrt(sigmai.sq[i]+psi1.sq)))
      pdf2=sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta, sqrt(sigmai.sq[i])))
      pdf3=sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+d2, sqrt(sigmai.sq[i]+psi2.sq)))
      r1i=pi1*pdf1/(pi1*pdf1+pi2*pdf2+pi3*pdf3); r1i[is.na(r1i)]=0
      r2i=pi2*pdf2/(pi1*pdf1+pi2*pdf2+pi3*pdf3); r2i[is.na(r2i)]=0
      r3i=pi3*pdf3/(pi1*pdf1+pi2*pdf2+pi3*pdf3); r3i[is.na(r3i)]=0
      
      ## M-step
      pi1_new=mean(r1i, na.rm = T); pi2_new=mean(r2i, na.rm = T); pi3_new=mean(r3i, na.rm = T)
      delta_new=sum(r1i*(Delta-d1)/(sigmai.sq+psi1.sq)+r2i*Delta/sigmai.sq+
                      r3i*(Delta-d2)/(sigmai.sq+psi2.sq), na.rm = T)/
        sum(r1i/(sigmai.sq+psi1.sq)+r2i/sigmai.sq+r3i/(sigmai.sq+psi2.sq), na.rm = T)
      d1_new=min(sum(r1i*(Delta-delta)/(sigmai.sq+psi1.sq), na.rm = T)/
                   sum(r1i/(sigmai.sq+psi1.sq), na.rm = T), 0)
      d2_new=max(sum(r3i*(Delta-delta)/(sigmai.sq+psi2.sq), na.rm = T)/
                   sum(r3i/(sigmai.sq+psi2.sq), na.rm = T), 0)
      
      # Nelder-Mead simplex algorithm for psi1.sq, psi2.sq, and sigmai.sq
      obj.psi1.sq=function(x){
        log.pdf=log(sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+d1, sqrt(sigmai.sq[i]+x))))
        log.pdf[is.infinite(log.pdf)]=0
        -sum(r1i*log.pdf, na.rm = T)
      }
      psi1.sq_new=neldermead(x0 = psi1.sq, fn = obj.psi1.sq, lower = 0)$par
      
      obj.psi2.sq=function(x){
        log.pdf=log(sapply(seq(n.taxa), function(i) dnorm(Delta[i], delta+d2, sqrt(sigmai.sq[i]+x))))
        log.pdf[is.infinite(log.pdf)]=0
        -sum(r3i*log.pdf, na.rm = T)
      }
      psi2.sq_new=neldermead(x0 = psi2.sq, fn = obj.psi2.sq, lower = 0)$par
      
      ## Merge to the paras vectors/matrices
      pi1.vec=c(pi1.vec, pi1_new); pi2.vec=c(pi2.vec, pi2_new); pi3.vec=c(pi3.vec, pi3_new)
      delta.vec=c(delta.vec, delta_new)
      d1.vec=c(d1.vec, d1_new); d2.vec=c(d2.vec, d2_new)
      psi1.sq.vec=c(psi1.sq.vec, psi1.sq_new); psi2.sq.vec=c(psi2.sq.vec, psi2.sq_new)
      
      ## Calculate the new epsilon
      epsilon=sqrt((pi1_new-pi1)^2+(pi2_new-pi2)^2+(pi3_new-pi3)^2+(delta_new-delta)^2+
                     (d1_new-d1)^2+(d2_new-d2)^2+(psi1.sq_new-psi1.sq)^2+(psi2.sq_new-psi2.sq)^2)
      iterNum=iterNum+1
    }
    # 2.23 Estimate the bias term
    bias.est.vec[i]=delta.vec[length(delta.vec)]
  }
  bias.est.vec=c(0, bias.est.vec)
  
  ### 3. Data transformation: remove the between-group bias
  mu.nonstrc.zero.adj=t(t(mu)+bias.est.vec)
  
  ### 4. Permutation test settings
  ## 4.1 Values of the test statistics for the current data set
  mu.nonstrc.zero.adj.diff=t(apply(mu.nonstrc.zero.adj, 1, function(x) combn(x, 2, FUN=diff)))
  mu.var.diff=t(apply(mu.var, 1, function(x) combn(x, 2, FUN=sum)))
  grp.pair=combn(n.grp, 2)
  colnames(mu.nonstrc.zero.adj.diff)=sapply(1:ncol(grp.pair), function(x) 
    paste0(grp.name[grp.pair[2, x]], " - ", grp.name[grp.pair[1, x]]))
  colnames(mu.var.diff)=sapply(1:ncol(grp.pair), function(x) 
    paste0(grp.name[grp.pair[2, x]], " - ", grp.name[grp.pair[1, x]]))
  
  stat.each=mu.nonstrc.zero.adj.diff/sqrt(mu.var.diff)
  stat.each[is.na(stat.each)]=0 # Replace missing values with 0s
  stat.pool=apply(abs(stat.each), 1, max)
  
  ## 4.2 Permutation data sets
  stat.pool.sim.list=lapply(1:perNum, function(x) {
    set.seed(x)
    mu.nonstrc.zero.adj.sim=matrix(rnorm(n.taxa*n.grp), nrow = n.taxa, ncol = n.grp)*sqrt(mu.var)
    mu.nonstrc.zero.adj.diff.sim=t(apply(mu.nonstrc.zero.adj.sim, 1, function(x) combn(x, 2, FUN=diff)))
    
    stat.each.sim=mu.nonstrc.zero.adj.diff.sim/sqrt(mu.var.diff)
    stat.each.sim[is.na(stat.each.sim)]=0
    stat.pool.sim=apply(abs(stat.each.sim), 1, max)
    return(stat.pool.sim)
  })
  stat.pool.sim=Reduce('cbind', stat.pool.sim.list)
  
  ### 5. Test result
  p.val.nonstrc.zero=apply(stat.pool.sim-stat.pool, 1, function(x) sum(x>0)/perNum)
  q.val.nonstrc.zero=p.adjust(p.val.nonstrc.zero, method = adj.method)
  q.val.nonstrc.zero[is.na(q.val.nonstrc.zero)]=1
  
  p.val=rep(0, n.taxa.origin); p.val[info.taxa.pos]=p.val.nonstrc.zero
  q.val=rep(0, n.taxa.origin); q.val[info.taxa.pos]=q.val.nonstrc.zero
  
  mu.adj=matrix(NA, nrow = n.taxa.origin, ncol = n.grp)
  mu.adj[info.taxa.pos, ]=mu.nonstrc.zero.adj
  if(length(info.taxa.pos)<n.taxa.origin){
    O.strc.zero=feature.table[-info.taxa.pos, ]
    ind.strc.zero=struc.zero[-info.taxa.pos, rep(1:n.grp, times = n.samp.grp)]
    O.strc.zero.adj=O.strc.zero*(1-ind.strc.zero)
    y.strc.zero=log(O.strc.zero.adj+1)
    d.strc.zero=t(t(1-ind.strc.zero)*d)
    y.strc.zero.adj=y.strc.zero-d.strc.zero
    mu.strc.zero=t(apply(y.strc.zero.adj, 1, function(i) 
      tapply(i, rep(1:n.grp, n.samp.grp), function(j) mean(j, na.rm = T))))
    # Make it the relative mean difference (with related to the smallest value)
    mu.strc.zero.adj=mu.strc.zero
    mu.strc.zero.adj[mu.strc.zero.adj==0]=NA
    mu.strc.zero.adj=t(t(mu.strc.zero.adj)+abs(apply(mu.strc.zero, 2, min)))
    mu.strc.zero.adj[is.na(mu.strc.zero.adj)]=0
  }else{
    mu.strc.zero.adj=NA
  }
  mu.adj[-info.taxa.pos, ]=mu.strc.zero.adj
  colnames(mu.adj)=paste0("mean_abundance", grp.name)
  rownames(mu.adj)=taxa.id.origin
  
  res = data.frame(mu.adj, p.val, q.val, diff.abn=ifelse(q.val<0.05, 1, 0))
  
  out=list(res=res, d=d, bias.est=bias.est.vec)
  return(out)
}
