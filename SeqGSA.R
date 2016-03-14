# compute the maxmean statistic for a set of z values
maxmean = function(z){
  n = length(z)
  pos = sum(z[z > 0])/n
  neg = -sum(z[z < 0])/n
  max = max(pos, neg)
  list(pos=pos,neg=neg,max=max)
}

# SeqGSA function
# x, table of count
# y, sample group
# z.func, function to compute gene-level z statistics, z.func(x,y), output is a vector of z values.
# genesets, list of vectors, each vector contains the row numbers of the genes
# length, vector of gene length
# weighted = T by default. If weighted = F, the non-weighted sampling is used.
# pow = 1 by default. Tuning parameter sampling weights. pow > 1 puts more weights on genes sets of similar length, 
# pow < 1 puts less weights on gene sets of similar length.
# restand.basis, use the all genes in the data or the catalog of genesets for restandardization.
# restand.first = F by default. If restand.first = T, the maxmean statistics are restandardized by positive and negative part first and then take the maximum value. 
# nrand, number of randomly sampled sets.
SeqGSA_func = function(x, y, z.func, genesets, len = 1, weighted = T, pow = 1, restand.first = F,
                restand.basis = c("all","catalog"), nrand = 1000) {
  n = nrow(x)
  # setsizes = sapply(genesets, length)
  if(restand.basis=="all") {
    freq = rep(1/n,n)
  }
  if(restand.basis=="catalog") {
    catalog = unlist(genesets)
    freq = tabulate(catalog, n)/length(catalog)
  }
  z = z.func(x,y)
  mm = data.frame(do.call("rbind", lapply(genesets, function(set) unlist(maxmean(z[set])))))
  cdf.all = rank(len)/sum(!is.na(len))
  cdf.all[is.na(len)] = NA
  
  s.restand = rep(NA, length(genesets))
  for (i in 1:length(genesets)) {
    geneset = genesets[[i]]
    size = length(geneset)
    if (!weighted) {
      weights = freq
    } else {
      cdf.set = cdf.all[genesets[[i]]]
      dist = outer(cdf.all, cdf.set, function(l1,l2) 1 - abs(l1 - l2))
      avg = colSums(dist, na.rm = T)
      dist = t(apply(dist,1,function(x) x/avg))
      dist[is.na(dist)] = 1/n
      weights = freq*rowSums(dist)
      weights = weights/sum(weights)
      weights = weights^pow
      weights = weights/sum(weights)
    }
    
    if (restand.first) {
      s = mm[i,1:2]
      mu.dagger.pos = sum(z*(z>0)*weights)
      mu.dagger.neg = -sum(z*(z<0)*weights)
      sd.dagger.pos = sqrt(sum((z*(z>0)-mu.dagger.pos)^2*weights)/size)
      sd.dagger.neg = sqrt(sum((z*(z<0)-mu.dagger.neg)^2*weights)/size)
      s.pos.restand = (s$pos - mu.dagger.pos)/sd.dagger.pos
      s.neg.restand = (s$neg - mu.dagger.neg)/sd.dagger.neg
      s.restand[i] = ifelse(s.pos.restand > s.neg.restand, s.pos.restand, -s.neg.restand)
    } else {
      s = mm[i,3]
      rand.sets = replicate(nrand, sample(n, size, prob = weights))
      rand.s = apply(rand.sets, 2, function(set) maxmean(z[set])$max)
      mu.dagger = mean(rand.s)
      sd.dagger = sd(rand.s)
      s.restand[i] = (s - mu.dagger)/sd.dagger
    }
  }
  s.restand
}

SeqGSA = function(x, y, z.func, genesets, len, weighted = T, pow = 1, restand.first = F, restand.basis = "catalog", 
                nrand = 1000, nperm = 1000, packages = c("edgeR"), funs = c("maxmean"), ncores = 8) {
  s.set = wGSA_func(x, y, z.func, genesets, len, weighted, pow, restand.first, restand.basis, nrand)
  
  cl = makeCluster(ncores, type = "SOCK")
  registerDoSNOW(cl)
  s.perm = foreach(i = 1:nperm, .combine = "rbind", .export = c("SeqGSA_func", funs), .packages = packages) %dopar% {
    SeqGSA_func(x, sample(y), z.func, genesets, len, weighted, pow, restand.first, restand.basis, nrand)
  }
  stopCluster(cl)
  p.hi = sapply(1:length(genesets), function(i) sum(s.perm[,i] > s.set[i])/nperm)
  p.lo = sapply(1:length(genesets), function(i) sum(s.perm[,i] < s.set[i])/nperm)
  list(p.hi,p.lo)
}

