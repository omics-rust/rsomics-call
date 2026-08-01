# rsomics-call

`rsomics-call` is the alignment-likelihood and lightweight small-variant
calling product in the rsomics portfolio.

The crate remains unpublished while its release performance gates are being
completed. Its implemented product surface contains three commands:

```text
rsomics-call pileup
rsomics-call call
rsomics-call run
```

`pileup` writes likelihood VCF/BCF from one or more coordinate-sorted
alignments, `call` consumes likelihood VCF/BCF, and `run` executes the same
typed stages without serializing an intermediate file. For example:

```sh
rsomics-call pileup -f reference.fa alignments.bam -Ou -o likelihoods.bcf
rsomics-call call likelihoods.bcf --model multiallelic -v -Ob -o calls.bcf
rsomics-call run -f reference.fa alignments.bam -v -Ob -o calls.bcf
```

All three commands use the shared `rsomics-help` command tree and
`rsomics-common` result and output contracts. Variant data may stream to
standard output, while named outputs are committed transactionally. `--json`
requires named variant output so the machine result envelope remains a
separate standard-output stream.

The current code establishes coordinate-merged SAM/BAM/CRAM input, multi-input
sample projection, bounded reference caching, and streaming SNP likelihoods
with per-input depth limits, deterministic deep-sample selection, and
allele-aligned site and sample quality evidence. An explicit reference-free
mode emits `N` reference alleles and rejects BAQ and indel configuration at the
builder boundary. Indexed region input merges
records across samples by coordinate, sorts selections by alignment-header
order, merges overlapping or adjacent intervals, and emits each selected site
once. Streaming inclusion targets require no alignment index, normalize
overlapping selections, preserve alignment order, and intersect with indexed
regions when both are present. Target files accept BED, one-based tabular, and
VCF coordinates with transparent plain, gzip, or BGZF input. The library
includes full and partial BAQ, insertion and
deletion likelihoods, sample-specific reference consensus, glocal realignment,
STR penalties, pooled or per-sample candidate support, explicit ambiguous-read
allele-depth policy, typed indel annotations, zero-copy, haploid, diploid, and
independently grouped multiallelic calls, a fused typed calling path, and strict
likelihood and called-record streams for plain VCF, BGZF VCF, raw BCF, and BGZF
BCF. The command boundary also binds alignment lists, sample projection,
indexed region files, streaming target files, flag filters, overlap and depth
policy, BAQ, indel policy, caller models, ploidy, grouped samples, gVCF blocks,
prior-frequency tags, and all four output encodings. Pileup records carry the
bcftools default bias metrics plus strand-aware
allele depths, quality means and sums, strand bias, mismatch, and soft-clip
annotations at sample and site scope. These paths are checked against bcftools
and HTSlib 1.24. Called records retain those annotations, replace the internal
`I16` evidence with model-specific `DP4`, `MQ`, and `PV4`, and preserve
allele-indexed fields when callers trim alternates. Typed ploidy resolution
supports constant, GRCh37, GRCh38, and checked custom definitions with
sex-aware or fixed sample assignments and allocation-free repeated queries. It
also groups consecutive reference calls into bcftools-compatible gVCF depth
blocks with typed `END` and `MIN_DP` output. Likelihood sample selection is
projected during record decoding, while call sample files preserve requested
order and bind explicit sex or fixed ploidy through the same definition.
Multiallelic workflows can select alleles independently for validated sample
groups, retain unused alternate dimensions at variant sites, and incorporate
validated integer panel allele counts from configurable prior-frequency INFO
tags. Prior counts remain aligned when output alleles are projected. The call
stream applies masked-reference, variant-only, and SNP/indel skip policy before
serialization with bcftools-compatible record-type semantics. Call sample
files currently accept exactly one sample name per row; the contradictory
bcftools 1.24 optional sex and numeric-ploidy behavior is not exposed as a
command-line promise. Call-stage targets, target complements, and a separate
unseen-allele switch are likewise absent until their documented and installed
bcftools 1.24 behaviors can be reconciled. The complete implemented commands
remain available while these isolated options and the product performance
gate are unresolved.

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

gVCF reference blocking follows bcftools 1.24 `gvcf.c` under the same retained
bcftools MIT notice.

The revised MAQ error model follows HTSlib 1.24 `errmod.c`, Copyright Broad
Institute and Genome Research Ltd., under the MIT/Expat license retained in
`LICENSES/HTSLIB-MIT.txt`. Its deterministic `drand48` sequence follows the
FreeBSD implementation carried by HTSlib; the notice is retained in
`LICENSES/FREEBSD-RAND.txt`. This product is licensed under MIT OR Apache-2.0.
