bins.raw <- read.table("bin-genotype.2013-08-12.txt", header = T, sep = "\t")
flank.mids <- read.table("flank-mids.tsv", header = T, sep = "\t")

bin.bounds <- read.table("sample-file/bins.tsv", header = T, sep = "\t")
length(bins$pos[bins$pos >= bin.bounds$start && bins$pos <= bin.bounds$end])

bins <- merge(flank.mids, bins.raw, by = c('chr', 'pos'))
bins <- bins[, c(1:4, order + 4)]




bin.bounds[bins$chr[1] == bin.bounds$chr & bins$pos[1] >= bin.bounds$start & bins$pos[1] <= bin.bounds$end, ]




bounds <- function(chr, mid, df = bin.bounds) {
  hit <- df[df$chr   == chr &
            df$start <= mid &
            df$end   >= mid, ]
  c(hit$start, hit$end)
}
bounds('A10', 17507149)


bins$chr[1:10]
bins$pos[1:10]



flanks <- bounds(bins$chr[1:10], bins$pos[1:10])
flanks <- bounds(bins$chr[32:35], bins$pos[32:35])
flanks <- bounds(bins$chr, bins$pos)
c( bins$start, bins$end ) <- matrix(flanks, ncol = 2)
bins$start <- matrix(flanks, ncol = 2)[,1]



library(ggplot2)
library(reshape)

bins.m <- melt(bins, id = c('chr', 'pos', 'start', 'end'), variable_name = 'BIL')
bins.m$bil.idx <- as.integer(bins.m$BIL)
head(bins.m)
qplot()


#test
bins.m <- melt(bins[1:100, 1:20], id = c('chr', 'pos'), variable_name = 'BIL')
head(bins.m)
dim(bins.m)
qplot(pos, y = 1, color = value, data = bins.m) + facet_grid(BIL ~ chr)
qplot(pos, y = 1, color = value, data = bins.m, geom = 'tile') + facet_grid(BIL ~ chr)
qplot(xmin = pos - 100, xmax = pos + 100, ymin = 0, ymax = 1, color = value, data = bins.m, geom = 'rect') + facet_grid(BIL ~ chr)
qplot(xmin = pos - 10000, xmax = pos + 10000, ymin = 0, ymax = 1, fill = value, data = bins.m, geom = 'rect') + facet_grid(BIL ~ chr)



qplot(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = value, color = value, data = bins.m, geom = 'rect') + facet_grid(BIL ~ chr, scales = "free_x")
qplot(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = value, color = value, data = bins.m[bins.m$BIL == 'RIL_1', ], geom = 'rect') + facet_grid(BIL ~ chr, scales = "free_x")
qplot(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = value, color = value, data = bins.m[1:1365, ], geom = 'rect') + facet_grid(BIL ~ chr, scales = "free_x")
qplot(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = value, color = value, data = bins.m[1:6825, ], geom = 'rect') + facet_grid(BIL ~ chr, scales = "free_x")


qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m, geom = 'rect') + facet_grid(. ~ chr, scales = "free_x") + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black"))
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m, geom = 'rect') + facet_grid(. ~ chr) + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black"))
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m, geom = 'rect') + facet_grid(. ~ chr, scales = "free_x", space = "free_x") + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black"))

#THIS ONE
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m, geom = 'rect') + facet_grid(. ~ chr, scales = "free_x", space = "free_x") + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black")) + theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())  + labs(colour = "Genotype", fill = "Genotype")
ggsave("ril-bins.png")

qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value,
      data = bins.m[bins.m$chr == 'A01', ], geom = 'rect') + facet_grid(. ~ chr, scales = "free_x", space = "free_x") + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black")) + theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())  + labs(colour = "Genotype", fill = "Genotype")

qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m[bins.m$BIL == 'RIL_103' | bins.m$BIL == 'RIL_136' , ], geom = 'rect') + facet_grid(. ~ chr, scales = "free_x", space = "free_x") + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black")) + theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())  + labs(colour = "Genotype", fill = "Genotype")
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m[bins.m$BIL %in% colnames(bins)[106:125] , ], geom = 'rect') + facet_grid(. ~ chr, scales = "free_x", space = "free_x") + scale_colour_manual(values = c("orange", "sky blue", "black")) + scale_fill_manual(values = c("orange", "sky blue", "black")) + theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())  + labs(colour = "Genotype", fill = "Genotype")
ggsave("ril-bins.zoom20.png")



qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m, geom = 'rect') + facet_grid(. ~ chr, scales = "free_x")
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m[bins.m$BIL == 'RIL_1', ], geom = 'rect') + facet_grid(. ~ chr, scales = "free_x")
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m[1:1365, ], geom = 'rect') + facet_grid(. ~ chr, scales = "free_x")
qplot(xmin = start, xmax = end, ymin = bil.idx, ymax = bil.idx + 1, fill = value, color = value, data = bins.m[1:6825, ], geom = 'rect') + facet_grid(. ~ chr, scales = "free_x")

ggplot()



bins <- read.table("bin-genotype.2013-08-12.txt", header = T, sep = "\t", as.is = T)
bins[bins == "R500"] <- 0
bins[bins == "IMB211"] <- 1


dist(t(bins[, 3:20]), upper = T, diag = F, method = "binary")
distances <- dist(t(bins[, 3:20]), method = "binary")



distances <- dist(t(bins[, 3:ncol(bins)]), method = "binary")
hc <- hclust(distances)
plot(hc)
plot(hc, hang = -1)


#hc method: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".


dist.mat <- as.matrix(dist(distances))
heatmap(dist.mat)






order <- hc$order
colnames(bins[, order + 2])












# junk
set.seed(123)
x <- rnorm(10)
hc <- hclust(dist(x))
dd <- as.dendrogram(hc)
dd.reorder <- reorder(dd, 10:1)
plot(dd, main = "random dendrogram 'dd'")


