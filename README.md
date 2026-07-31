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

The current code establishes the typed likelihood boundary, multi-input sample
projection, and SNP error model against bcftools and HTSlib 1.24. It does not
yet expose a command-line binary.
Incomplete commands stay absent until the complete multi-input, SNP/indel,
BAQ, consensus/multiallelic, VCF/BCF, oracle, and performance gates pass.

The historical `rsomics-vcf-mpileup` and `rsomics-vcf-call` repositories are
implementation and fixture sources, not dependencies. Their single-sample
shells are not the product boundary.

The revised MAQ error model follows HTSlib 1.24 `errmod.c`, Copyright Broad
Institute and Genome Research Ltd., under the MIT/Expat license retained in
`LICENSES/HTSLIB-MIT.txt`. This product is licensed under MIT OR Apache-2.0.
