# @author dylan sosa
setwd('/Users/Dylan/Documents/SLU/5.1/Bioinformatics 1')
library('edgeR')
countsSRR = read.delim('/Users/Dylan/Documents/SLU/5.1/Bioinformatics 1/counts_out', header = FALSE, row.names = 1)
# change row names to transcript:
head(countsSRR)
names(countsSRR)[1] = 'SRR2567786'
names(countsSRR)[2] = 'SRR2567795'

groups = c('SRR2567786','SRR2567795')
d = DGEList(counts = countsSRR, group = factor(groups))
d

d.full = d
head(d$counts)
head(cpm(d))

apply(d$counts, 2, sum)

keep = rowSums(cpm(d)>1) >= 2
# trimming 
# 1

d = d[keep,]
dim(d)
# two groups 

d$samples$lib.size = colSums(d$counts)
d$samples

d = calcNormFactors(d)
d

require(DESeq)
cds = newCountDataSet(data.frame(d$counts), d$samples$group)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds, method = 'blind')
plotDispEsts(cds)

bcv = 0.2
et = exactTest(d, dispersion=bcv^2)
de1 = decideTestsDGE(et, adjust.method = 'BH', p.value = 0.05)
# FDR is 0.05
summary(de1)

de1tags = rownames(d)[as.logical(de1)]
plotSmear(et, de.tags = de1tags)
abline(h = c(-2,2), col = 'green')
