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
