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


PlotCompositeMap <- function(bin.genotypes.melted,
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
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
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

par1 <- "M82"
par2 <- "PEN"

order <- GetOrderOfSamplesClusteredByGenotype(bins.physical,
                                              par1 = par1, par2 = par2)

bins.physical.m <- ClusterAndMeltBinGenotypes(bins.physical, order,
                                              par1 = par1, par2 = par2)
bins.genetic.m <- ClusterAndMeltBinGenotypes(bins.genetic, order,
                                              par1 = par1, par2 = par2)

PlotCompositeMap(bins.physical.m, par1 = par1, par2 = par2, col1 = "magenta",
  col2 = "green", plot.file = "plots/composite-map.physical.png", save = TRUE,
  plot=FALSE, chr.text.size = 12, width = 10, height = 7.5)
PlotCompositeMap(bins.genetic.m, par1 = par1, par2 = par2, col1 = "magenta",
  col2 = "green", plot.file = "plots/composite-map.genetic.png", save = TRUE,
  plot=FALSE, chr.text.size = 12, width = 10, height = 7.5)
