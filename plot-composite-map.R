GetOrderOfSamplesClusteredByGenotype <- function(bin.genotypes,
                                       par1 = "par1", par2 = "par2") {
  library(plyr)

  bins.binary <- colwise(as.character)(bin.genotypes)
  bins.binary[bins.binary == par1] <- 0
  bins.binary[bins.binary == par2] <- 1
  bins.binary[bins.binary == "HET"] <- 1

  distances <- dist(t(bins.binary[, 5:ncol(bins.binary)]), method = "binary")
  hc <- hclust(distances)
  clustered.sample.order <- hc$order

  clustered.sample.order
}


ClusterAndMeltBinGenotypes <- function(bin.genotypes, order,
                                       par1 = "par1", par2 = "par2") {
  library(reshape)

  samples.sorted <- bin.genotypes[, c(1:4, order + 4)]

  bin.genotypes.m <- melt(samples.sorted,
                          id = c("chr", "bin.mid", "bin.start", "bin.end"),
                          variable_name = "sample")
  bin.genotypes.m$value <- factor(bin.genotypes.m$value,
                                  levels = c(par2, "HET", par1))
  bin.genotypes.m$sample.idx <- as.integer(bin.genotypes.m$sample)

  bin.genotypes.m
}


SubsetByChrAndHasIntrogression <- function(bin.genotypes, chromosome) {
  chr.subset <- subset(bin.genotypes, chr == chromosome)
  has.introgression <- colwise(
    function(z) {length(unique(z)) > 1})(chr.subset[5:ncol(chr.subset)])
  chr.subset[, c(rep(TRUE, 4), as.vector(has.introgression, mode="logical"))]
}


PlotCompositeMap <- function(bin.genotypes.melted, stacked.chromosomes = FALSE,
                             par1 = "par1", par2 = "par2",
                             col1 = "sky blue", colh = "black", col2 = "orange",
                             plot.file = "composite-map.png",
                             plot = TRUE, save = FALSE,
                             chr.text.size = 7, chr.text.angle = 0,
                             ggtitle = "Composite Genotype Map", ...) {
  library(ggplot2)

  composite.map <- ggplot(bin.genotypes.melted) +
    geom_rect(aes(
      xmin = bin.start,
      xmax = bin.end,
      ymin = sample.idx,
      ymax = sample.idx + 1,
      fill = value,
      color = value
    )) +
    scale_colour_manual(values = c(col2, colh, col1)) +
    scale_fill_manual(values = c(col2, colh, col1)) +
    theme(
      strip.text.x = element_text(
        size = chr.text.size,
        angle = chr.text.angle
      ),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(colour = "Genotype", fill = "Genotype") +
    ggtitle(ggtitle)

  if (stacked.chromosomes) {
    composite.map <- composite.map +
      facet_grid(chr ~ ., scales = "free_y", space = "free_y")
  } else {
    composite.map <- composite.map +
      facet_grid(. ~ chr, scales = "free_x", space = "free_x")
  }

  if (plot)
    print(composite.map)

  if (save)
    ggsave(filename = plot.file, plot = composite.map, ...)
}


setwd("/Users/mfc/tomato/dan/bil-composite-plot-genetic-distance/")

bins.physical.file <- "data/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like"
bins.genetic.file <- "data/gen.bins.info.RDS"

bins.physical <- read.table(bins.physical.file, header = T, sep = "\t")
bins.genetic <- readRDS(bins.genetic.file)

bins.physical$chr <- bins.genetic$chr

par1 <- "M82"
par2 <- "PEN"

order <- GetOrderOfSamplesClusteredByGenotype(bins.physical,
                                              par1 = par1, par2 = par2)

bins.physical.m <- ClusterAndMeltBinGenotypes(bins.physical, order,
                                              par1 = par1, par2 = par2)
bins.genetic.m <- ClusterAndMeltBinGenotypes(bins.genetic, order,
                                             par1 = par1, par2 = par2)

PlotCompositeMap(bins.physical.m, par1 = par1, par2 = par2, col1 = "magenta",
                 col2 = "green", plot.file = "plots/composite-map.physical.png",
                 save = TRUE, plot=FALSE, chr.text.size = 12, width = 10,
                 height = 7.5)
PlotCompositeMap(bins.genetic.m, par1 = par1, par2 = par2, col1 = "magenta",
                 col2 = "green", plot.file = "plots/composite-map.genetic.png",
                 save = TRUE, plot=FALSE, chr.text.size = 12, width = 10,
                 height = 7.5)


# Cluster by chromosome after removing samples without an introgression
bins.physical.m <- data.frame()
bins.genetic.m <- data.frame()

for (chromosome in unique(bins.physical$chr)) {
  bins.physical.chr <- SubsetByChrAndHasIntrogression(bins.physical, chromosome)
  bins.genetic.chr <- SubsetByChrAndHasIntrogression(bins.genetic, chromosome)

  order <- GetOrderOfSamplesClusteredByGenotype(bins.physical.chr,
                                                par1 = par1, par2 = par2)

  bins.physical.chr.m <- ClusterAndMeltBinGenotypes(bins.physical.chr, order,
                                                    par1 = par1, par2 = par2)
  bins.genetic.chr.m <- ClusterAndMeltBinGenotypes(bins.genetic.chr, order,
                                                   par1 = par1, par2 = par2)

  bins.physical.m <- rbind(bins.physical.m, bins.physical.chr.m)
  bins.genetic.m <- rbind(bins.genetic.m, bins.genetic.chr.m)
}

PlotCompositeMap(bins.physical.m, stacked.chromosomes = TRUE,
                 par1 = par1, par2 = par2, col1 = "magenta", col2 = "green",
                 plot.file = "plots/composite-map.physical.cluster-by-chr.png",
                 save = TRUE, plot=FALSE, chr.text.size = 12, width = 7.5,
                 height = 10)
PlotCompositeMap(bins.genetic.m, stacked.chromosomes = TRUE,
                 par1 = par1, par2 = par2, col1 = "magenta", col2 = "green",
                 plot.file = "plots/composite-map.genetic.cluster-by-chr.png",
                 save = TRUE, plot=FALSE, chr.text.size = 12, width = 7.5,
                 height = 10)
