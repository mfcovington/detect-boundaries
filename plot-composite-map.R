PlotCompositeMap <- function(bounds.dir, bins.file, plot.file="composite-map.png") {
  library(ggplot2)
  library(reshape)

  bins.raw <- read.table(bins.file, header = T, sep = "\t")

  bins <- colwise(as.character)(bins.raw)
  par1 <- "M82"
  par2 <- "PEN"
  bins[bins == par1] <- 0
  bins[bins == par2] <- 1
  bins[bins == "HET"] <- 1

  distances <- dist(t(bins[, 5:ncol(bins)]), method = "binary")
  hc <- hclust(distances)
  order <- hc$order
  bins.sorted <- bins.raw[, c(1:4,order + 4)]

  bins.m <- melt(bins.sorted, id = c('chr', 'bin.mid', 'bin.start', 'bin.end'), variable_name = 'BIL')
  bins.m$value <- factor(bins.m$value, levels = c(par2, "HET", par1))
  bins.m$bil.idx <- as.integer(bins.m$BIL)

  qplot(
    xmin = bin.start,
    xmax = bin.end,
    ymin = bil.idx,
    ymax = bil.idx + 1,
    fill = value,
    color = value,
    data = bins.m,
    geom = 'rect'
  ) +
  facet_grid(
    . ~ chr,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_colour_manual(values = c("green", "black", "magenta")) +
  scale_fill_manual(values = c("green", "black", "magenta")) +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(colour = "Genotype", fill = "Genotype")

  ggsave(plot.file)
}

bounds.dir <- "~/Dropbox/lab/tomato/bils/re-seq/hmm-results/boundaries-fixed"
bins.file <- "/Users/mfc/Dropbox/lab/tomato/bils/re-seq/hmm-results/boundaries-fixed/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like"

PlotCompositeMap(bounds.dir, bins.file)
