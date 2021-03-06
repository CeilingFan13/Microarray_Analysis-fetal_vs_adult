#!/usr/bin/env Rscript
# microarray_analysis.R

library("limma")
library("affy")
library(knitr)

# get CEL files
affy.data <- ReadAffy()

# converts an instance of AffyBatch into an instance of ExpressionSet
# standardize for mas5 probe set
eset.mas5 <- mas5(affy.data)
# extracts the normalized expression matrix
exprSet.nologs <- exprs(eset.mas5)
# name columns of the matrix
colnames(exprSet.nologs) <-
  c(
    "brain.1",
    "brain.2",
    "fetal.brain.1",
    "fetal.brain.2",
    "fetal.liver.1",
    "fetal.liver.2",
    "liver.1",
    "liver.2"
  )

# binary representation (log base 2 binary dataset) represent fold change
exprSet <- log(exprSet.nologs, 2)
kable(head(exprSet))
# export tab-separated file if necessary
write.table(exprSet,
            file = "Summary_mas5_matrix.txt",
            quote = F,
            sep = "\t")

#-------------------------------------------------------------------------------

# making P/M/A calls
data.mas5calls <- mas5calls(affy.data)
# extracts expression matrix containing A and P for each gene and tissue
data.mas5calls.calls <- exprs(data.mas5calls)
# load in microarray image file
# green for brain, red for fetal brain
brain.fetalbrain.2color <-
  read.maimages(
    "brain.fetalbrain.2color.data.txt",
    columns = list(
      G = "brain.1",
      R = "fetal.brain.1",
      Gb = "bg1",
      Rb = "bg2"
    )
  )
# standardize microarray image using local regression
brain.fetalbrain.2color.loess <-
  normalizeWithinArrays(brain.fetalbrain.2color,
                        method = "loess")
# two panel graph, 1 row
par(mfrow = c(1, 2))
# plot original data and corrected data
plotMA(brain.fetalbrain.2color)
plotMA(brain.fetalbrain.2color.loess)

#-------------------------------------------------------------------------------

# calculate mean of each set
brain.mean <- apply(exprSet[, c("brain.1", "brain.2")], 1, mean)
fetal.brain.mean <-
  apply(exprSet[, c("fetal.brain.1", "fetal.brain.2")],
        1, mean)
liver.mean <- apply(exprSet[, c("liver.1", "liver.2")], 1, mean)
fetal.liver.mean <-
  apply(exprSet[, c("fetal.liver.1", "fetal.liver.2")],
        1, mean)
# calculate ratio: log(A/B) = log(A) - log(B)
brain.fetal.to.adult <- fetal.brain.mean - brain.mean
liver.fetal.to.adult <- fetal.liver.mean - liver.mean
# dataframe with all the information from this section
# mean value of each data set and ratio
all.data <- cbind(
  exprSet,
  brain.mean,
  fetal.brain.mean,
  liver.mean,
  fetal.liver.mean,
  brain.fetal.to.adult,
  liver.fetal.to.adult
)
# write tab-separated file
write.table(all.data,
            file = "Microarray_ALL.txt",
            quote = F,
            sep = "\t")

#-------------------------------------------------------------------------------

# two tailed t test betwee two brain and two liver
brain.p.value.all.genes <-
  apply(exprSet, 1, function(x) {
    t.test(x[1:2], x[3:4])$p.value
  })
liver.p.value.all.genes <-
  apply(exprSet, 1, function(x) {
    t.test(x[5:6], x[7:8])$p.value
  })
# place Absent/Present values in AP, showing status in each cell
AP <- apply(data.mas5calls.calls, 1, paste, collapse = "")
# make sure the cell has at least one presence
genes.present = names(AP[AP != "AAAAAAAA"])
# takes rows from exprSet where the gene has values in at least one chip
exprSet.present <- exprSet[genes.present, ]
# new dataset with p-vals
brain.raw.pvals.present <- brain.p.value.all.genes[genes.present]
liver.raw.pvals.present <- liver.p.value.all.genes[genes.present]
# false discovery rate correction for type I error
brain.fdr.pvals.present <-
  p.adjust(brain.raw.pvals.present, method = "fdr")
liver.fdr.pvals.present <-
  p.adjust(liver.raw.pvals.present, method = "fdr")
# sort in increasing order, most significant to least significant
brain.fdr.pvals.present.sorted <-
  brain.fdr.pvals.present[order(brain.fdr.pvals.present)]
liver.fdr.pvals.present.sorted <-
  liver.fdr.pvals.present[order(liver.fdr.pvals.present)]
# p-val less than 0.01
brain.DE.probesets <-
  names(brain.raw.pvals.present[brain.raw.pvals.present < 0.01])
liver.DE.probesets <-
  names(liver.raw.pvals.present[liver.raw.pvals.present < 0.01])
# get only ratio information of p-val less than 0.01
brain.DE.log2.ratios <-
  all.data[brain.DE.probesets, c("brain.fetal.to.adult", "liver.fetal.to.adult")]
liver.DE.log2.ratios <-
  all.data[liver.DE.probesets, c("brain.fetal.to.adult", "liver.fetal.to.adult")]
# write out the report
write.table(
  brain.DE.log2.ratios,
  "brain.DE.log2.ratios.txt",
  sep = "\t",
  quote = F
)
write.table(
  liver.DE.log2.ratios,
  "liver.DE.log2.ratios.txt",
  sep = "\t",
  quote = F
)

#-------------------------------------------------------------------------------

par(mfrow = c(1, 1))
x.data <- all.data[, "brain.mean"]
y.data <- all.data[, "fetal.brain.mean"]
# log2 expression in fetal vs adult brain
plot(
  x.data,
  y.data,
  main = "Log2 expression in fetal brain (n=2) vs adult brain (n=2)",
  xlab = "brain",
  ylab = "fetal brain",
  col = "blue",
  cex = 0.5
)
abline(0,1)
# create vocalno plot
average = (x.data + y.data) / 2
ratio = x.data - y.data
plot(
  average,
  ratio,
  main = "MA plot of fetal brain vs brain",
  pch = 19,
  cex = 0.2,
  col = "red"
)
# differential expressed genes not on line with slope = 0
abline(0,0)

# raw and FDR corrected p-values for all of the "present" genes
expression.plus.pvals = cbind(
  exprSet.present,
  brain.raw.pvals.present,
  brain.fdr.pvals.present,
  liver.raw.pvals.present,
  liver.fdr.pvals.present
)
log2.ratios = expression.plus.pvals[, "brain.1"] -  expression.plus.pvals[, "fetal.brain.1"]
p.values = expression.plus.pvals[, "brain.raw.pvals.present"]

# setting up two-column graph
par(mfrow=c(1,2), mtext("Title for Two Plots", outer = TRUE, cex = 1.5))
# p-val vs fold change
plot(log2.ratios, p.values)
# log base 10 transform of the raw p-values
plot(log2.ratios, -log(p.values, 10) )


