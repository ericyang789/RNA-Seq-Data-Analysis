library(edgeR)
infile <- "merged.23"
group <- factor(c(1,1,1,2,2,2))
outfile <-"result23_norm.out"
x <- read.table(infile, sep="	", row.names=1)
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tab <- topTags(et, nrow(x))
write.table(tab, file=outfile)
