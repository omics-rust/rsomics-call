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
sample projection, bounded reference caching, and a streaming SNP likelihood
path with per-input depth limits, deterministic deep-sample selection, and
allele-aligned site and sample quality evidence. The current slice is checked
against bcftools and HTSlib 1.24 and includes zero-copy, haploid, diploid, and
independently grouped multiallelic calls, a fused typed calling path, and
strict likelihood and called-record streams for plain VCF, BGZF VCF, raw BCF,
and BGZF BCF. The consensus allele-frequency posterior, diploid and haploid
genotypes, allele trimming, site quality, and genotype quality match bcftools
1.24 on the implemented oracle matrix. It does not yet expose a command-line
binary.
Incomplete commands stay absent until indel and BAQ likelihoods, consensus
annotations, regions and targets, gVCF behavior, complete command orchestration,
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

The revised MAQ error model follows HTSlib 1.24 `errmod.c`, Copyright Broad
Institute and Genome Research Ltd., under the MIT/Expat license retained in
`LICENSES/HTSLIB-MIT.txt`. Its deterministic `drand48` sequence follows the
FreeBSD implementation carried by HTSlib; the notice is retained in
`LICENSES/FREEBSD-RAND.txt`. This product is licensed under MIT OR Apache-2.0.
