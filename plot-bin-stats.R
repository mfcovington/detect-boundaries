CountAndMeasureIntrogressions <- function(boundaries.dir,
                                          par1 = "par1", par2 = "par2") {
  # Accumulate counts and lengths data for introgressions

  counts.df <- data.frame(id = character(),
                          par1 = integer(),
                          het = integer(),
                          par2 = integer())
  lengths.df <- counts.df

  filelist <- list.files(boundaries.dir, pattern = ".+\\.boundaries")

  for (filename in filelist) {
    id <- sub("(^.+)\\.boundaries$", "\\1", filename)

    df <- read.table(paste0(boundaries.dir, filename),
                     colClasses = c("factor", rep("numeric", 2), "factor"))
    colnames(df) <- c("chr", "start", "end", "genotype")
    df$length <- df$end - df$start + 1

    counts <- count(df, "genotype")
    par1_counts <- max(counts$freq[counts$genotype == par1], 0)
    het_counts  <- max(counts$freq[counts$genotype == "HET"], 0)
    par2_counts <- max(counts$freq[counts$genotype == par2], 0)

    counts.df <- rbind(counts.df, data.frame(id = id,
                                             par1 = par1_counts,
                                             het = het_counts,
                                             par2 = par2_counts))

    lengths <- aggregate(length ~ genotype, data = df, sum)
    par1_lengths <- max(lengths$length[lengths$genotype == par1], 0)
    het_lengths  <- max(lengths$length[lengths$genotype == "HET"], 0)
    par2_lengths <- max(lengths$length[lengths$genotype == par2], 0)

    lengths.df <- rbind(lengths.df, data.frame(id = id,
                                               par1 = par1_lengths,
                                               het = het_lengths,
                                               par2 = par2_lengths))
  }

  counts.df  <<- counts.df
  lengths.df <<- lengths.df
}


PlotIntrogressionsPerSample <- function (
      counts.df, border.color = "black", fill.color = "skyblue",
      xlab = "# of introgressions per sample", ylab = "# of Samples",
      plot.file = "introgressions-per-sample.png",
      plot = TRUE, save = FALSE, ...) {

  # Run CountAndMeasureIntrogressions() to generate input data frame

  library(ggplot2)

  max.count.introgression.combined <- max(counts.df$het + counts.df$par2)

  introgression.histogram <- ggplot(data = counts.df) +
    geom_histogram(aes(x = het + par2),
                   binwidth = 1,
                   origin = 0.5,
                   color = border.color,
                   fill = fill.color) +
    xlim(0, max.count.introgression.combined + 0.5) +
    xlab(xlab) +
    ylab(ylab)

  if (plot)
    print(introgression.histogram)

  if (save)
    ggsave(filename = plot.file, plot = introgression.histogram, ...)
}


PlotPercentIntrogressed <- function(
      counts.df, border.color = "black", fill.color = "skyblue",
      xlab = "% of introgressed genotype in genome", ylab = "# of Samples",
      plot.file = "percent-introgressed.png",
      plot = TRUE, save = FALSE, ...) {

  # Run CountAndMeasureIntrogressions() to generate input data frame

  library(ggplot2)

  introgression.histogram <- ggplot(data = lengths.df) +
    geom_histogram(aes(x = 100 * (par2 + het) / (par1 + het + par2)),
                   binwidth = 0.15,
                   color = border.color,
                   fill = fill.color) +
    xlab(xlab) +
    ylab(ylab) +
    scale_x_log10()

  if (plot)
    print(introgression.histogram)

  if (save)
    ggsave(filename = plot.file, plot = introgression.histogram, ...)
}


PlotBinsPerChromosome <- function(
      bin.geno.df, border.color = "black", fill.color = "skyblue",
      plot.file = "bins-per-chromosome.png", plot = TRUE, save = FALSE,
      ggtitle = "Bins per Chromosome", ...) {

  library(ggplot2)

  bins.per.chr <- ggplot(bin.geno.df, aes(x = chr)) +
    geom_bar(color = border.color, fill = fill.color) +
    ggtitle(ggtitle) +
    xlab("Chromosome") +
    ylab("Number of unique bins")

  if (plot)
    print(bins.per.chr)

  if (save)
    ggsave(filename = plot.file, plot = bins.per.chr, ...)
}


PlotDistributionOfIntrogressions <- function(
      bin.geno.df, genetic.distance = FALSE,
      par1 = "par1", par2 = "par2",
      color.introgression = "orange", color.het = "black",
      plot.file = "distribution-of-introgressions.png",
      plot = TRUE, save = FALSE,
      chr.text.size = 12, chr.text.angle = 270,
      ggtitle = "Distribution of Introgressions Across Bins", ...) {

  library(ggplot2)

  bin.geno.df$par1 <- apply(bin.geno.df, 1, function(line) sum(line == par1))
  bin.geno.df$het  <- apply(bin.geno.df, 1, function(line) sum(line == "HET"))
  bin.geno.df$par2 <- apply(bin.geno.df, 1, function(line) sum(line == par2))

  if (!genetic.distance) {
    # added due to broken feature in ggplot 0.9.1:
    # is it still broken/required?
    max_pos <- max(bin.geno.df$bin.end)
    temp_max <- floor(max_pos / 10000000)
    max_pos_fixed <- temp_max * 10000000
    label_max <- temp_max * 10
  }

  offset <- 1.05
  max_count <-  max(bin.geno.df$het + bin.geno.df$par2)
  max_lab <- (max_count %/% 10) * 10

  distribution <- ggplot(bin.geno.df, aes(xmin = bin.start, xmax = bin.end)) +
    geom_rect(aes(ymin = 0, ymax = par2), fill = color.introgression) +
    geom_rect(aes(ymin = par2, ymax = het + par2), fill = color.het) +
    facet_grid(chr ~ .) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    ggtitle(ggtitle) +
    scale_y_continuous(
      '# of BILs with introgression',
      breaks = c(0, max_lab / 2, max_lab),
      labels = c(0, max_lab / 2, max_lab),
      limits = c(0, offset * max_count)
    ) +
    theme(
      strip.text.y = element_text(
        size = chr.text.size,
        angle = chr.text.angle
      )
    )

  if (genetic.distance) {
    distribution <- distribution + xlab('Bin position on chromosome (cM)')
  } else {
    distribution <- distribution +
      scale_x_continuous(
        'Bin position on chromosome (Mb)',
        breaks = seq(0, max_pos_fixed, 10000000),
        labels = seq(0, label_max,     10)
      )
  }

  if (plot)
    print(distribution)

  if (save)
    ggsave(filename = plot.file, plot = distribution, ...)
}
