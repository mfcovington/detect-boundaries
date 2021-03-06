#!/usr/bin/env Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)

id <- args[1]

if (!exists("id"))
  stop("No Sample ID provided!")

# if (exists("args[2]")) {
  par1_id <- args[2]
# } else {
  # par1_id <- 'R500'
# }

# if (exists("args[3]")) {
  par2_id <- args[3]
# } else {
#   par2_id <- 'IMB211'
# }

print("running")

# setwd("/Users/mfc/git.repos/detect-boundaries")
library(ggplot2)

snps.dir <- "filtered-snps/"
snps.files <- list.files(path = snps.dir, pattern = id)
# snps.files <- list.files(path = "filtered-snps/", pattern = paste0(id, "*.snps"))
# print(snps.files)
# print(paste0("filtered-snps/", id, "*.filtered.snps"))
# q()
snps.df <- data.frame(chr = character(0), pos = integer(0), cov = integer(0), par = character(0))

for (file in snps.files) {
  file <- paste0(snps.dir, file)
  print(paste("Reading", file))
  snps.in <- read.table(file, sep = '\t')
  colnames(snps.in) <- c('chr', 'pos', 'cov', 'par')
  snps.df <- rbind(snps.df, snps.in)
}

snps.df$score[snps.df$par == par1_id] <- -snps.df$cov[snps.df$par == par1_id]
snps.df$score[snps.df$par == par2_id] <-  snps.df$cov[snps.df$par == par2_id]
# qplot(pos, score, data = snps.df, geom = 'line')
# qplot(pos, score, data = snps.df, geom = 'area')

snps.df$binary[snps.df$par == par1_id] <- -1
snps.df$binary[snps.df$par == 'HET']   <- 0
snps.df$binary[snps.df$par == par2_id] <- 1
# qplot(pos, binary, data = snps.df, geom = 'line')
# qplot(pos, binary, data = snps.df, geom = 'area')


# ggplot(snps.df, aes(x = pos, y = score)) + facet_grid(chr ~ .) + geom_area()
boundary.plot <- ggplot(snps.df, aes(x = pos, y = binary)) +
                   facet_grid(chr ~ .) +
                   geom_area(color = 'blue', fill = 'orange')
# ggplot(snps.df, aes(x = pos, y = binary)) + facet_grid(chr ~ .) + geom_line()

# ggsave(filename = paste(id, "png", sep = "."),
#        plot     = boundary.plot,
#        width    = 10,
#        height   = 8)

plot.name <- paste0("plots/", id, ".png")
print(paste("Saving", plot.name))
ggsave(filename = plot.name,
       plot     = boundary.plot,
       width    = 10,
       height   = 8)
