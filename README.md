


`merge-replicates.pl`: Use if there are multiple replicates for a single sample and the genotype files need to be merged.

- Input: Genotype files (one per chromosome per replicates for a sample)
- Output: Merged genotype files (one per chromosome)

`filter-snps.pl`: 

- Input:
    - Sample ID
    - Genotype files (one per chromosome for a sample)
- Output: 
    - Boundaries file
    - Filtered SNP files (one per chromosome)


### Notes:

detect-boundaries:`filter-snps.pl`

- IN:
    - Genotyped files:
- OUT:
    - Boundary files: `sample-file/RIL_300.boundaries`, etc.
    - Filtered SNP files: `RIL_300.A01.filtered.snps`, etc.


### Merging Boundaries with `merge-boundaries.pl`

This script takes a collection of `*.boundaries` files and merges all of the boundaries. The bins resulting from the merge are written to a file: `bins.tsv`. Boundary and bin stats (count, min size, max size, and mean size) are printed to the screen.

Sample Usage:

```shell
./merge-boundaries.pl \
  --bam_file ~/git.repos/sample-files/bam/IMB211.good.bam \
  --chr_list A01,A02,A03,A04,A05,A06,A07,A08,A09,A10 \
  sample-file/RIL_*.boundaries
```

- IN:
    - Boundary files: `sample-file/RIL_300.boundaries`, etc.
    - Representative `.bam` file (for extracting chromosome lengths)
    - (Optional) Chromosome list (for limiting analysis to a subset of sequences)
- OUT: 
    - Comprehensive list of bins and their locations: `bins.tsv`

detect-boundaries:`extract-snp-subset.R`

- IN:
    - VCF summaries: `/Users/mfc/git.repos/snps-from-rils/merged.EC50.minflter.vcf/summaries.het_ratio_0_1.alt_ratmin_0_05.filters/merged.*.EC50.minflter.vcf.summary`
    - Bins: `bins.tsv`
- OUT: 
    - `sample-file/bins-snp.mid`


detect-boundaries:`get-genotypes-for-snp-subset.pl`

- IN:
    - `sample-file/bins-snp.mid`
    - Genotyped files:
    - polyDB...
    - list of Sample IDs
- OUT:
    - `bin-genotype`

*Version 0.6.0*
