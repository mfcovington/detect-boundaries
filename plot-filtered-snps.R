setwd("/Users/mfc/git.repos/detect-boundaries")
library(ggplot2)

par1_id <- 'R500'
par2_id <- 'IMB211'
snps.in <- read.table('snps.out', sep = '\t')
colnames(snps.in) <- c('chr', 'pos', 'cov', 'par')

snps.in$score[snps.in$par == par1_id] <- -snps.in$cov[snps.in$par == par1_id]
snps.in$score[snps.in$par == par2_id] <-  snps.in$cov[snps.in$par == par2_id]
qplot(pos, score, data = snps.in, geom = 'line')
qplot(pos, score, data = snps.in, geom = 'area')

snps.in$binary[snps.in$par == par1_id] <- -1
snps.in$binary[snps.in$par == par2_id] <- 1
qplot(pos, binary, data = snps.in, geom = 'line')
qplot(pos, binary, data = snps.in, geom = 'area')

