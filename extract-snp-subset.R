


vcf.dir <- "/Users/mfc/git.repos/snps-from-rils/merged.EC50.minflter.vcf/summaries.het_ratio_0_1.alt_ratmin_0_05.filters/"
vcf.summary.files <- paste0(vcf.dir, list.files(vcf.dir, "merged.*.EC50.minflter.vcf.summary"))
snps <- read.table(c(vcf.summary.files), header = TRUE, sep = "\t", as.is = TRUE)[c(1:3, 11)]
snps <- snps[snps$filter == '.', 2:4]



read.multi.tables <- function(file.names, ...) {
    require(plyr)
    ldply(file.names, function(fn) data.frame(filename = fn, read.table(fn, ...)))
}

snps.multi <- read.multi.tables(vcf.summary.files, header = TRUE, sep = "\t", as.is = TRUE)[c(2:4, 12)]
snps.multi <- snps.multi[snps.multi$filter == '.', 2:4]


vcf.summary <- "/Users/mfc/git.repos/snps-from-rils/merged.EC50.minflter.vcf/summaries.het_ratio_0_1.alt_ratmin_0_05.filters/merged.A01.EC50.minflter.vcf.summary"
snps <- read.table(vcf.summary, header = TRUE, sep = "\t", as.is = TRUE)[c(1:3, 11)]
snps <- snps[snps$filter == '.', 2:4]

bins <- read.table("bins.tsv", header = TRUE, sep = "\t", as.is = TRUE)

chr    <- bins$chr[1]
start  <- bins$start[1]
end    <- bins$end[1]
length <- end - start + 1

# extract observed ratio data for bin
region <- snps[snps$chr == chr & snps$pos >= start & snps$pos <= end, ]


left  <- region$pos[1]
right <- region$pos[length(region$pos)]

flank       <- 0.2
flank_extra <- flank + 0.1

flank_adj       <- length * flank
flank_extra_adj <- length * flank_extra

obs_rat_min       <- 0.9
obs_rat_min_extra <- 0.8

left_flank        <- left  + flank_adj
right_flank       <- right - flank_adj
left_flank_extra  <- left  + flank_extra_adj
right_flank_extra <- right - flank_extra_adj




left.snp <- (region$pos[region$observed_ratio == 1 &
                    region$pos <= left_flank, ])[1] # change left_flank

if (is.na(left.snp))
  left.snp <- (region$pos[region$observed_ratio == 1 &
                      region$pos <= left_flank, ])[1] # change left_flank

if (is.na(left.snp$pos))
  left.snp <- (region$pos[region$observed_ratio == 1 &
                      region$pos <= left_flank, ])[1] # change left_flank

if (is.na(left.snp$pos))
  left.snp <- (region$pos[region$observed_ratio >= obs_rat_min &
                      region$pos <= left_flank, ])[1]

if (is.na(left.snp$pos))
  left.snp <- (region$pos[region$observed_ratio >= obs_rat_min_extra &
                      region$pos <= left_flank, ])[1]

if (is.na(left.snp$pos))
  left.snp <- (region$pos[region$observed_ratio >= obs_rat_min &
                      region$pos <= left_flank_extra, ])[1]

if (is.na(left.snp$pos))
  left.snp <- (region$pos[region$observed_ratio >= obs_rat_min_extra &
                      region$pos <= left_flank_extra, ])[1]

if (is.na(left.snp$pos))
  left.snp <- (region$pos[region$observed_ratio >= obs_rat_min_extra &
                      region$pos <= left_flank_extra, ])[1]

if (is.na(left.snp$pos))


if (is.na(left.snp$pos))


if (!exists("good"))
  print("sorry")



region$pos[which.max(region$observed_ratio[region$pos <= left_flank]), ]


if (!exists("good"))
  good <- (region[region$pos <= left_flank & region$observed_ratio == 2, ])[1, 2]
if (is.na(good$pos))
  print("sorry")


get.mid.max <- function(region) {
  max.pos <- region$pos[region$observed_ratio == max(region$observed_ratio)]
  max.length <- length(max.pos)
  if (max.length == 1) {
    max.mid <- max.pos
  } else {
    max.mid <- max.pos[round(max.length / 2)]
  }
  max.mid
}

system.time(
  for (i in 1:nrow(bins)) {
  # for (i in 1:100) {
    chr    <- bins$chr[i]
    start  <- bins$start[i]
    end    <- bins$end[i]
  print(i)
    # extract observed ratio data for bin
    region <- snps[snps$chr == chr & snps$pos >= start & snps$pos <= end, ]

    max.mid <- get.mid.max(region)
    # print(bins[i, ])
    # print(max.mid)
    # print(snps[snps$pos == max.mid, ])
    bins$max.mid[i] <- max.mid
  }
)






