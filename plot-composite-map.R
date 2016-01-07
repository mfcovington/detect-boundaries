ClusterSamplesByGenotype <- function(bins.physical, bins.genetic,
                                     par1 = "par1", par2 = "par2") {
  library(plyr)
  library(reshape)

  bins.binary <- colwise(as.character)(bins.physical)
  bins.binary[bins.binary == par1] <- 0
  bins.binary[bins.binary == par2] <- 1
  bins.binary[bins.binary == "HET"] <- 1

  distances <- dist(t(bins.binary[, 5:ncol(bins.binary)]), method = "binary")
  hc <- hclust(distances)
  order <- hc$order
  samples.sorted <- bins.physical[, c(1:4, order + 4)]

  bins.physical.m <- melt(samples.sorted,
                          id = c("chr", "bin.mid", "bin.start", "bin.end"),
                          variable_name = "sample")
  bins.physical.m$value <- factor(bins.physical.m$value,
                                  levels = c(par2, "HET", par1))
  bins.physical.m$sample.idx <- as.integer(bins.physical.m$sample)

  bins.genetic.sorted <- cbind(bins.genetic[, 1:4],
                               samples.sorted[, 5:ncol(samples.sorted)])
  bins.genetic.m <- melt(bins.genetic.sorted,
                         id = c("chr", "bin.mid", "bin.start", "bin.end"),
                         variable_name = "sample")
  bins.genetic.m$value <- factor(bins.genetic.m$value,
                                 levels = c(par2, "HET", par1))
  bins.genetic.m$sample.idx <- as.integer(bins.genetic.m$sample)

  bins.physical.m <<- bins.physical.m
  bins.genetic.m <<- bins.genetic.m
}


PlotCompositeMap <- function(bin.genotypes,
                             par1 = "par1", par2 = "par2",
                             col1 = "sky blue", colh = "black", col2 = "orange",
                             plot.file = "composite-map.png",
                             plot = TRUE, save = FALSE,
                             chr.text.size = 7, chr.text.angle = 0,
                             ggtitle = "Composite Genotype Map", ...) {
  library(ggplot2)

  composite.map <- ggplot(bin.genotypes) +
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

ClusterSamplesByGenotype(bins.physical, bins.genetic, par1 = par1, par2 = par2)

PlotCompositeMap(bins.physical.m, par1 = par1, par2 = par2, col1 = "magenta",
  col2 = "green", plot.file = "plots/composite-map.physical.png", save = TRUE,
  plot=FALSE, chr.text.size = 12, width = 10, height = 7.5)
PlotCompositeMap(bins.genetic.m, par1 = par1, par2 = par2, col1 = "magenta",
  col2 = "green", plot.file = "plots/composite-map.genetic.png", save = TRUE,
  plot=FALSE, chr.text.size = 12, width = 10, height = 7.5)
