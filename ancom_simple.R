ancom.W = function(otu_data, var_data, main.var, sig){
  
  n_otu=dim(otu_data)[2]
  otu_ids=colnames(otu_data)
  
  n.samp=nrow(otu_data)
  group=as.factor(var_data[, grep(main.var, colnames(var_data))])
  grp1.ind=which(group==levels(group)[1])
  grp2.ind=which(group==levels(group)[2])
  n.samp1=length(grp1.ind)
  n.samp2=length(grp2.ind)
  
  # Log of the abundance
  log_otu_data = apply(otu_data, 2, function(x) log(1+x))
  # Log-ratio of the abundance
  lr = apply(log_otu_data, 2, function(x) log_otu_data[, 1]-x)
  for (i in 2:n_otu){
    lr_new = apply(log_otu_data, 2, function(x) log_otu_data[, i]-x)
    lr=cbind(lr, lr_new)
  }
  pval=apply(lr, 2, function(x) wilcox.test(x[1:n.samp1], x[n.samp1+1:n.samp])$p.value)
  pval[which(is.nan(pval))]=1
  pval.mat=matrix(pval, byrow = T, ncol = n_otu)
  
  pval.adj.mat=apply(pval.mat, 2, function(x) p.adjust(x, method = "BH"))
  W = apply(pval.adj.mat, 2, function(x) sum(x<sig))
  
  return(W)
  }



ANCOM.main = function(OTUdat, Vardat, main.var, sig, prev.cut){
  
  p.zeroes=apply(OTUdat,2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  OTUdat.thinned=OTUdat[, which(p.zeroes<prev.cut)]
  
  otu.names=colnames(OTUdat.thinned)
  
  W.detected = ancom.W(otu_data=OTUdat.thinned, var_data=Vardat, main.var, sig)
  W_stat = W.detected
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE

  return(W_frame)
}

