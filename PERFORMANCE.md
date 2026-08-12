# Performance

The release gate compares the fused `rsomics-call run` path with the equivalent
`bcftools mpileup | bcftools call` pipeline. The benchmark is single-threaded,
uses raw BCF on both sides, and verifies every reported call before accepting a
timing round.

## Release benchmark

The 2026-08-02 release benchmark used `rsomics-call` revision
`85579cb94f9aacc91bc878447b9241ef7e20c58e` and bcftools/HTSlib 1.24 on an
Ubuntu 22.04 machine with two Intel Xeon Gold 6238R processors. Both commands
were pinned to CPU 20. One warm-up preceded five alternating timed rounds.

```text
rsomics-call run -f reference.fa -v -O u -o ours.bcf alignments.bam
bcftools mpileup -f reference.fa -a FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/QM,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/FS,INFO/NMBZ,INFO/NM,INFO/SCR -Ou alignments.bam |
  bcftools call -mv -Ou -o bcftools.bcf
```

| Tool | Median wall time | Mean wall time | Range | Peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-call run` | 25.64 s | 25.70 s | 25.52–26.09 s | 22,400 KiB |
| `bcftools mpileup \| bcftools call` | 26.05 s | 26.19 s | 25.98–26.82 s | 47,488 KiB |

| Round | `rsomics-call` | bcftools |
|---:|---:|---:|
| 1 | 26.09 s | 26.09 s |
| 2 | 25.52 s | 26.03 s |
| 3 | 25.67 s | 26.05 s |
| 4 | 25.64 s | 25.98 s |
| 5 | 25.59 s | 26.82 s |

On this fixture, `rsomics-call` had 1.6% lower median wall time and 52.8% lower
peak RSS. This is evidence for this workload, not a claim that the same ratio
holds across data sets or machines.

The deterministic fixture contains a 5 Mb reference and 500,000 pairs of
150-base simulated reads (about 30x coverage). It was generated with wgsim
1.24 using:

```text
-S 20260802 -N 500000 -1 150 -2 150 -d 350 -s 50 -e 0.002 -r 0.001 -R 0.15 -X 0.30
```

Reads were aligned with minimap2 2.31 and coordinate-sorted and indexed with
samtools 1.24. The input checksums were:

```text
13bd65f4568d0a30bc0ee218db62223cc26d9593f2b116530aa5e0b78b5f34dc  reference.fa
fb41379eeef996a50dc4aa58104d696b76c731b8bfdb5a8584fc348f9041f952  reference.fa.fai
33b6780ec3758a8ccde746935366dec441e89aaafb5b0253a19cfa1af350282c  alignments.bam
8f615ebabf218f338447cbb2dc6202f7478d81fdb72eae1219a58565a253981e  alignments.bam.bai
```

Each implementation emitted 5,024 calls in every round. The sorted
`CHROM/POS/REF/ALT/GT` projection was byte-identical, with SHA-256
`f5b6f9d78aba94c81d4c263acc969bafbf19284b278dc1d239915326d2ac43ee`.
The comparison intentionally excludes tool-specific headers and does not claim
byte identity for all annotations.

The tested `rsomics-call` binary was built in release mode with rustc 1.95.0
and had SHA-256
`1ee8fcb4b9a6d74cc1d45132d2f2195ec3d37acdd6e2b8a0eb58e328e7f08793`.
The bcftools binary had SHA-256
`50b8f25b1dd20ca1b0c9ffa76f0f2e4684515764ee4af1c190debd9ece490c5d`.

## Pileup record-state regression gate

Revision `8f29a887dc96` moves per-record CIGAR metrics into the retained state
provided by `rsomics-pileup` 0.9.0 and updates `rsomics-bamio` and
`rsomics-common` to their current foundation releases. A Linux `x86_64`
regression gate rebuilt both this revision and published release head
`b34cc226242ba` with rustc 1.91.0, pinned each command to CPU 20, and used the
same 5 Mb, 30x fixture and fused `run` command as the release benchmark.

Both binaries emitted 5,024 normalized calls with SHA-256
`f5b6f9d78aba94c81d4c263acc969bafbf19284b278dc1d239915326d2ac43ee`.
After one warm-up, five candidate rounds had a 34.081 s median and
34.583 ± 0.861 s mean. The published baseline had a 50.635 s median and
48.296 ± 10.573 s mean under high concurrent machine load. Three alternating
RSS rounds per binary were all 22,400 KiB. This gate establishes no regression
from the foundation migration; the lower-variance upstream comparison above
remains the release performance claim.

The candidate and baseline binary SHA-256 values were respectively
`e8acf12c12331a841f34452cbb2191156d01900aab60350395d9d851cd11c204`
and `233c1bfcf66bba6b6086f63751c950f45698d0ca6b0afdd9fcb62e450089f116`.
The Hyperfine JSON and RSS ledger have SHA-256
`d77775dd8480d44e28f6ef58e49f029c911ce7137975680418c8b7a50801b253`
and `eedc13e2f8da730ba7c1079ef1be6f6146c5b61f92ad74baf4572deb412da096`.

## Indexed reference consolidation regression gate

Revision `7d7bb20e64a07ac38ae9738691ac87dbf6e9234e` replaces the product-local
indexed FASTA cache with `rsomics-seqio` 0.6.0. It was compared with published
release 0.1.1, whose crate records VCS revision
`7d20d6b119dcaea60638cc7da793bca61a47fedc`. Both binaries were built with
rustc 1.91.0 and ran the full fused command on the same 5 Mb, 30x fixture.

Both emitted 5,024 calls. Their raw BCF output was byte-identical with SHA-256
`6c1e4e96ac3c1a7a8b8268364473585f22438c0ca552ece8dad9ec44b4d05c81`; the
normalized call projection also retained the release checksum
`f5b6f9d78aba94c81d4c263acc969bafbf19284b278dc1d239915326d2ac43ee`.

Hyperfine 1.20.0 ran one warm-up followed by two order-reversed batches of five
runs per binary on an Apple M2 Mac14,3 with 8 GiB RAM and macOS 26.6.1 build
25G76.

| Batch | Candidate mean | Published 0.1.1 mean |
|---|---:|---:|
| Candidate first | 14.698 ± 0.096 s | 14.778 ± 0.133 s |
| Baseline first | 14.666 ± 0.089 s | 14.702 ± 0.053 s |
| Combined | 14.682 s | 14.740 s |

Three alternating `/usr/bin/time -lp` rounds gave candidate mean and maximum
RSS of 21,561,344 and 23,363,584 bytes, versus 21,949,099 and 23,805,952 bytes
for 0.1.1. The timing and memory ranges overlap; this gate establishes no
regression and does not support a speed or memory improvement claim.

The candidate and baseline binary SHA-256 values were respectively
`3e57562efffbd942a456fadfc0894e0743ee4bdc8b664dc8c3ebfaf1682bc46a`
and `09a1a82cb58c47714e23f43beb2121a2652acd6c38fcf649492e4edd933a6da4`.
The forward Hyperfine JSON, reverse Hyperfine JSON, and RSS ledger have SHA-256
`4cae632c3f721a153e94299acbd1c87b5f7768c40fc84bd05fa1c789d40c9a53`,
`6bee43c9e4e14746cb7b3c579579e6ea2646a10325190fa872ffccc94c783e33`, and
`cd751f26926b22bda5e5b143412646de3f672176f906fe5cfd0ff51d04d437ed`.

## Reproduction

`benchmarks/call-vs-bcftools.sh` records the machine, tool versions, binary and
input checksums, per-round GNU time results, and normalized calls. It fails if
any round disagrees:

```sh
benchmarks/call-vs-bcftools.sh \
  target/release/rsomics-call \
  /path/to/bcftools \
  /path/to/reference.fa \
  /path/to/alignments.bam \
  /path/to/results \
  5 20
```
