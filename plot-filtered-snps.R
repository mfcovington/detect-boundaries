setwd("/Users/mfc/git.repos/detect-boundaries")
library(ggplot2)

id <- 'RIL_41'

snps.files <- list.files(pattern = paste0(id, ".*snps"))
par1_id <- 'R500'
par2_id <- 'IMB211'

snps.df <- data.frame(chr = character(0), pos = integer(0), cov = integer(0), par = character(0))

for (file in snps.files) {
  print(file)
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
ggplot(snps.df, aes(x = pos, y = binary)) + facet_grid(chr ~ .) + geom_area(color = 'blue', fill = 'orange')
# ggplot(snps.df, aes(x = pos, y = binary)) + facet_grid(chr ~ .) + geom_line()
