


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


detect-boundaries:`merge-boundaries.pl`

- IN:
    - Boundary files: `sample-file/RIL_300.boundaries`, etc.
    - Chromosome lengths
- OUT: 
    - `bins.tsv`


detect-boundaries:`extract-snp-subset.R`

- IN:
    - VCF summaries: `/Users/mfc/git.repos/snps-from-rils/merged.EC50.minflter.vcf/summaries.het_ratio_0_1.alt_ratmin_0_05.filters/merged.*.EC50.minflter.vcf.summary`
    - Bins: `bins.tsv`
- OUT: 
    - `sample-file/bins-snp.mid`


detect-boundaries:`get-genotypes-for-subset.pl`

- IN:
    - `sample-file/bins-snp.mid`
    - Genotyped files:
    - polyDB...
    - list of Sample IDs
- OUT:
    - `bin-genotype`


