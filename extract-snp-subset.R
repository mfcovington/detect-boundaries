


read.multi.tables <- function(file.names, ...) {
    require(plyr)
    ldply(file.names, function(fn) data.frame(filename = fn, read.table(fn, ...)))
}

vcf.dir <- "/Users/mfc/git.repos/snps-from-rils/merged.EC50.minflter.vcf/summaries.het_ratio_0_1.alt_ratmin_0_05.filters/"
vcf.summary.files <- paste0(vcf.dir, list.files(vcf.dir, "merged.*.EC50.minflter.vcf.summary"))
obs.rat <- read.multi.tables(vcf.summary.files, header = TRUE, sep = "\t", as.is = TRUE)[c(2:4, 12)]
obs.rat <- obs.rat[obs.rat$filter == '.', 2:4]

bins <- read.table("bins.tsv", header = TRUE, sep = "\t", as.is = TRUE)

polydb.files <- list.files("./", "polyDB.*")
snps <- read.multi.tables(polydb.files, header = TRUE, sep = "\t", as.is = TRUE)[2:3]

####### MIDDLES ######

get.mid.max <- function(region) {
  max.pos <- region$pos[region$observed_ratio == max(region$observed_ratio)]
  max.count <- length(max.pos)
  if (max.count == 1) {
    max.mid <- max.pos
  } else {
    mid   <- round(sum(range(region$pos)) / 2)
    diffs <- max.pos - mid

    min.diff.idx <- which.min(abs(diffs))

    max.mid <- mid + diffs[min.diff.idx]
  }
  max.mid
}

# bin.sizes <- 0
bin.size.min <- 50
for (i in 1:nrow(bins)) {
# for (i in 1:100) {
  chr    <- bins$chr[i]
  start  <- bins$start[i]
  end    <- bins$end[i]

  # ignore small bins
  bin.size <- end - start + 1
  # bin.sizes[i] <- bin.size
  if (bin.size == 1) {
    bins$snp.mid[i] <- NA
    next
  }

  # extract observed ratio data for bin
  region <- obs.rat[obs.rat$chr == chr & obs.rat$pos >= start & obs.rat$pos <= end, ]

  # ignore positions absent from polyDB files (i.e. insufficient coverage in parent)
  region$observed_ratio[!paste(region[, 1], region[, 2], sep = ':') %in% paste(snps[, 1], snps[, 2], sep = ':')] <- 0
  max.mid <- get.mid.max(region)
  bins$snp.mid[i] <- max.mid
}
write.table(bins, file = "bins-snp.mid", quote = FALSE, sep = "\t", row.names = FALSE)


####### FLANKS ######
chr    <- bins$chr[1]
start  <- bins$start[1]
end    <- bins$end[1]
length <- end - start + 1

# extract observed ratio data for bin
region <- obs.rat[obs.rat$chr == chr & obs.rat$pos >= start & obs.rat$pos <= end, ]


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



