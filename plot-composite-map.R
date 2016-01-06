PlotCompositeMap <- function(bins.file,
                             par1="par1", par2="par2",
                             col1="sky blue", colh="black", col2="orange",
                             plot.file="composite-map.png",
                             plot=FALSE, save=TRUE,
                             chr.text.size=7, chr.text.angle=0,
                             ggtitle="Composite Genotype Map", ...) {
  library(ggplot2)
  library(reshape)

  bins <- read.table(bins.file, header = T, sep = "\t")

  bins.binary <- colwise(as.character)(bins)
  bins.binary[bins.binary == par1] <- 0
  bins.binary[bins.binary == par2] <- 1
  bins.binary[bins.binary == "HET"] <- 1

  distances <- dist(t(bins.binary[, 5:ncol(bins.binary)]), method = "binary")
  hc <- hclust(distances)
  order <- hc$order
  bins.sorted <- bins[, c(1:4,order + 4)]

  bins.m <- melt(bins.sorted, id = c('chr', 'bin.mid', 'bin.start', 'bin.end'), variable_name = 'BIL')
  bins.m$value <- factor(bins.m$value, levels = c(par2, "HET", par1))
  bins.m$bil.idx <- as.integer(bins.m$BIL)

  composite.map <- ggplot(bins.m) +
  geom_rect(aes(
    xmin = bin.start,
    xmax = bin.end,
    ymin = bil.idx,
    ymax = bil.idx + 1,
    fill = value,
    color = value,
  )) +
  facet_grid(
    . ~ chr,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_colour_manual(values = c(col2, colh, col1)) +
  scale_fill_manual(values = c(col2, colh, col1)) +
  theme(
    strip.text.x = element_text(size = chr.text.size, angle = chr.text.angle),
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

bins.file <- "/Users/mfc/Dropbox/lab/tomato/bils/re-seq/hmm-results/boundaries-fixed/bin-genotypes.BILs.2014-12-07.imputed-NAs.merged-like"

PlotCompositeMap(bins.file, par1="M82", par2="PEN", col1="magenta", col2="green")
