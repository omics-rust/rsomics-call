# rsomics-call

`rsomics-call` is the alignment-likelihood and lightweight small-variant
calling product in the rsomics portfolio.

The repository is under implementation and is not published. Its release
contract contains three commands:

```text
rsomics-call pileup
rsomics-call call
rsomics-call run
```

The current code establishes coordinate-merged SAM/BAM/CRAM input, multi-input
sample projection, bounded reference caching, and streaming SNP likelihoods
with per-input depth limits, deterministic deep-sample selection, and
allele-aligned site and sample quality evidence. Indexed region input merges
records across samples by coordinate, sorts selections by alignment-header
order, merges overlapping or adjacent intervals, and emits each selected site
once. Streaming inclusion targets require no alignment index, normalize
overlapping selections, preserve alignment order, and intersect with indexed
regions when both are present. Target files accept BED, one-based tabular, and
VCF coordinates with transparent plain, gzip, or BGZF input. The library
includes full and partial BAQ, insertion and
deletion likelihoods, sample-specific reference consensus, glocal realignment,
STR penalties, typed indel annotations, zero-copy, haploid, diploid, and
independently grouped multiallelic calls, a fused typed calling path, and strict
likelihood and called-record streams for plain VCF, BGZF VCF, raw BCF, and BGZF
BCF. Pileup records carry the bcftools default bias metrics plus strand-aware
allele depths, quality means and sums, strand bias, mismatch, and soft-clip
annotations at sample and site scope. These paths are checked against bcftools
and HTSlib 1.24. Called records retain those annotations, replace the internal
`I16` evidence with model-specific `DP4`, `MQ`, and `PV4`, and preserve
allele-indexed fields when callers trim alternates. Typed ploidy resolution
supports constant, GRCh37, GRCh38, and checked custom definitions with
sex-aware or fixed sample assignments and allocation-free repeated queries. It
does not yet expose a command-line binary. Incomplete commands stay absent
until target-exclusion behavior, gVCF behavior, complete command orchestration,
and oracle and performance gates pass.

The historical `rsomics-vcf-mpileup` and `rsomics-vcf-call` repositories are
implementation and fixture sources, not dependencies. Their single-sample
shells are not the product boundary.

The multiallelic allele-selection, genotype, and quality model follows
bcftools 1.24 `mcall.c`, Copyright Genome Research Ltd., under its MIT/Expat
license retained in `LICENSES/BCFTOOLS-MIT.txt`.

The consensus allele-frequency posterior follows bcftools 1.24 `ccall.c` and
`prob1.c`; its attribution and license are retained in
`THIRD_PARTY_LICENSES.md`.

The established indel likelihood path follows bcftools 1.24
`bam2bcf_indel.c` and `str_finder.c`, and HTSlib 1.24 `probaln.c`; their
attribution and licenses are retained in `THIRD_PARTY_LICENSES.md`,
`LICENSES/BCFTOOLS-MIT.txt`, and `LICENSES/HTSLIB-MIT.txt`.

Pileup annotations and caller-side bias tests follow bcftools 1.24 `bam2bcf.c`,
`ccall.c`, and `mcall.c`; Fisher exact and related numerical kernels follow
HTSlib 1.24 `kfunc.c`. Their MIT notices are retained in the same license files.

Ploidy presets, custom definitions, and sample binding follow bcftools 1.24
`ploidy.c` and `vcfcall.c` under the retained bcftools MIT notice.

The revised MAQ error model follows HTSlib 1.24 `errmod.c`, Copyright Broad
Institute and Genome Research Ltd., under the MIT/Expat license retained in
`LICENSES/HTSLIB-MIT.txt`. Its deterministic `drand48` sequence follows the
FreeBSD implementation carried by HTSlib; the notice is retained in
`LICENSES/FREEBSD-RAND.txt`. This product is licensed under MIT OR Apache-2.0.
