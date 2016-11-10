install.packages("gplots")
install.packages("RColorBrewer")

require("gplots")
require("RColorBrewer")


z = exactZ(count,group)
z63 = z[cat2id[[63]]]
xy = DGEList(counts=count, group=group)
xy = calcNormFactors(xy)
mat_data <- count[cat2id[[63]],c(which(group==1),which(group==2))]
#mat_data = mat_data/xy$samples$norm.factors
tot.count = colSums(count)
TPM.factors = tot.count/exp(mean(log(tot.count)))
mat_data = mat_data/TPM.factors
names(mat_data) = c(paste0("female",1:17),paste("male",1:24))
mat_data = mat_data[order(z63, decreasing = T),]
rnames = row.names(mat_data)
mat_data = as.matrix(mat_data)
row.names(mat_data) = rnames
mat_data = log(mat_data+1)
#mat_data = (mat_data - 0.5*rowMeans(mat_data[,1:17]) - 0.5*rowMeans(mat_data[,18:41]))/apply(mat_data,1,sd)
mat_data = (mat_data - rowMeans(mat_data))/apply(mat_data,1,sd)
#mat_data[,1:17] = (mat_data[,1:17] - rowMeans(mat_data[,1:17]))/apply(mat_data[,1:17],1,sd) + rowMeans(mat_data[,1:17])
#mat_data[,18:41] = (mat_data[,18:41] - rowMeans(mat_data[,18:41]))/apply(mat_data[,18:41],1,sd) + rowMeans(mat_data[,18:41])
c1 = hclust(dist(t(mat_data[,1:17])), method = "complete", members = NULL)
ord1 = c1$order
c2 = hclust(dist(t(mat_data[,18:41])), method = "complete", members = NULL)
ord2 = c2$order + 17
#mat_data = mat_data[,c(ord1,ord2)]
mat_data = mat_data[,c(ord1,rev(ord2))]

# corr = cor(mat_data)
# diag(corr) = 0
# max.corr = apply(abs(corr),2,max)
# best.sample = sapply(1:nrow(corr), function(i) which.min(abs(corr[,i]-max.corr[i])))

quantile.range <- quantile(mat_data, probs = seq(0, 1, 0.01))
palette.breaks <- c(seq(quantile.range["5%"], 0, length = 10), seq(0, quantile.range["95%"], length = 10)[-1])
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], by=0.1)
# use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
# color.palette  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
color.palette  <- colorRampPalette(c("#C80000", "#FFFFDF", "#008000"))(length(palette.breaks) - 1)
color.palette  <- colorRampPalette(c("#008000", "#FFFFDF", "#B0171F"))(length(palette.breaks) - 1)
color.palette  <- colorRampPalette(c("#009ACD", "#BFEFFF"))(length(palette.breaks) - 1)

png("D:/Dropbox/thesis/GSA/Sim tools/heatmap2.png", width = 10*300, height = 10*300, res = 300, pointsize = 8)
heatmap.2(x = mat_data, Rowv=NULL, Colv=NULL, 
          #col = rev(rainbow(20*10, start = 0/6, end = 4/6)), 
          col = color.palette,
          breaks = palette.breaks,
          scale="none",
          margins=c(8,6), # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='none', 
          denscol="black",
          keysize=5, 
          #key.xlab="log-tranformed count (row scaled)",
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(4.0,0,4.0,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), 
          lhei=c(0.5, 4.5), lwid=c(1, 10, 1),
          cexRow = 1.9, cexCol = 0.1, labCol = NA, labRow = "ENSG00000154620"
          )
axis(1, 15.8/41, "", las=1, cex.axis= 1, lty = 1, lwd.ticks = 2, tck = -0.02)
axis(1, c(8/41,27/41), c("female","male"), las=1, cex.axis= 2, lty = 0)
dev.off()
