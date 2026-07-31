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
against bcftools and HTSlib 1.24 and includes the single-group zero-copy,
haploid, and diploid multiallelic caller core. It does not yet expose a
command-line binary.
Incomplete commands stay absent until SNP/indel and BAQ likelihoods, complete
consensus and multiallelic policies, VCF/BCF, oracle, and performance gates
pass.

The historical `rsomics-vcf-mpileup` and `rsomics-vcf-call` repositories are
implementation and fixture sources, not dependencies. Their single-sample
shells are not the product boundary.

The multiallelic allele-selection, genotype, and quality model follows
bcftools 1.24 `mcall.c`, Copyright Genome Research Ltd., under its MIT/Expat
license retained in `LICENSES/BCFTOOLS-MIT.txt`.

The revised MAQ error model follows HTSlib 1.24 `errmod.c`, Copyright Broad
Institute and Genome Research Ltd., under the MIT/Expat license retained in
`LICENSES/HTSLIB-MIT.txt`. Its deterministic `drand48` sequence follows the
FreeBSD implementation carried by HTSlib; the notice is retained in
`LICENSES/FREEBSD-RAND.txt`. This product is licensed under MIT OR Apache-2.0.
