x = count; y = group; z.func = exactZ; genesets = cat2id; len = len; weighted = T; pow = 1; packages="edgeR"
restand.first = F; restand.basis = "all"; nrand = 1000; nperm = 1000; funs = c("exactZ","maxmean"); ncores = 8

s.set = SeqGSA_func(x, y, z.func, genesets, len, weighted, pow, restand.first, restand.basis, nrand)
s.perm = foreach(i = 1:nperm, .combine = "rbind", .export = c("SeqGSA_func", funs), .packages = packages) %dopar% {
  wGSA_func(x, sample(y), z.func, genesets, len, weighted, pow, restand.first, restand.basis, nrand)
}
weighted = F
s2.set = SeqGSA_func(x, y, z.func, genesets, len, weighted, pow, restand.first, restand.basis, nrand)
s2.perm = foreach(i = 1:nperm, .combine = "rbind", .export = c("SeqGSA_func", funs), .packages = packages) %dopar% {
  wGSA_func(x, sample(y), z.func, genesets, len, weighted, pow, restand.first, restand.basis, nrand)
}

hist(s.perm[,63], freq=F)
hist(s2.perm[,63], freq=F, add=T, lty=3)
s.set[63]
s2.set[63]

histRes0 <- hist(s.perm[,63],breaks=seq(min(s.perm[,63]),max(s.perm[,63]),length=20),plot=FALSE)
xvals0 <- histRes0$breaks
yvals0 <- histRes0$density
xvals0 <- c(xvals0,xvals0[length(xvals0)])
yvals0 <- c(0,yvals0,0)

histRes1 <- hist(s2.perm[,63],breaks=seq(min(s2.perm[,63]),max(s2.perm[,63]),length=20),plot=FALSE)
xvals1 <- histRes1$breaks
yvals1 <- histRes1$density
xvals1 <- c(xvals1,xvals1[length(xvals1)])
yvals1 <- c(0,yvals1,0)

pdf("PermutationDensity.pdf")
par(las=1)
plot(xvals0,yvals0,type="S", xlab = "Standardized Maxmean", ylab = "Density", xlim = c(-2.7,4), ylim = c(0,0.5), col="red", lty = 2, xaxt="n", 
     cex.lab=1.0, cex.axis=1.0)
lines(xvals1,yvals1,type="S", lty = 2)
axis(1,at=c(-2,0,2),labels=c(-2,0,2), lwd=1, cex.axis=1.0)
axis(1,at=s.set[63],labels=expression('  s'[w]^'*'), col.ticks="red", lwd=2, padj = 0.25, cex.axis=1.5)
axis(1,at=s2.set[63],labels=expression('s'^'*'), lwd=2, cex.axis=1.5)
legend("topright", c("weighted permutation","unweighted permutation"), bty="n", col = c("red", "black"), lty=c(2,2), cex = 1.2)
dev.off()
