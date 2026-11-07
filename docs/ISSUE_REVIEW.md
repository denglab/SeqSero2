# SeqSero2 software issue review (2026-11-07)

This local patch was evaluated against baseline commit ff1dce923159fe6a7a0ea2ff24cf1ddb2000071a (v1.3.2). It is not a published release. The separate issue-by-issue report includes all 50 issues and comments (25 open, 25 closed), drafted replies and status recommendations; no comments or issue statuses were changed.

## Changes

- Literal subprocess arguments and Python filesystem operations, safe temporary paths for SPAdes/BLAST, fresh output directories and checked tool failures. This ports the reviewed SeqSero2S reliability changes, related to its issue #5.
- CLI input/type validation, importable helpers, accurate input provenance, dependency checks, SalmID table validation, single-end mapping and compressed input handling.
- Empty-primary-O-list guard for [#43](https://github.com/denglab/SeqSero2/issues/43) and [#54](https://github.com/denglab/SeqSero2/issues/54). These samples now finish with unresolved O evidence instead of crashing.
- Partial/full O9,46 wzy markers count as one contamination family ([#45](https://github.com/denglab/SeqSero2/issues/45)); separate O families still trigger the warning.
- Header-aware TSV summary helper ([#40](https://github.com/denglab/SeqSero2/issues/40)).
- SeqSero2S issue #6 marker fix for [#58](https://github.com/denglab/SeqSero2/issues/58): two matched context probes, maximum state scores, explicit unresolved variant and final-contig k-mer boundary correction. Existing 379 reference sets and original FASTA bytes are preserved; two sets are added. The negative context probe is a designed SNP counterpart, not an observed natural ST28 negative isolate.
- Modern setuptools build, no distutils runtime, independent version constant, documented dependency minimums and a proposed Conda recipe. The external Bioconda builds remain unchanged.

## Validation

Native macOS ARM64, Python 3.12.16, Biopython 1.88, BWA 0.7.19, samtools 1.25, BLAST 2.17.0, SPAdes 4.3.0, SalmID 0.1.23, SRA Toolkit 3.4.1 and seqtk 1.5 were used. The existing seqsero environment supplied dependencies. The final wheel was installed into a fresh workspace target with no dependency installation; the existing SeqSero2S installation was preserved.

- 38 automated regression tests pass: `python -m unittest discover -s tests -v`.
- Final scientific dataset: 123 real workflow executions, comprising 86 patched-source and 37 baseline executions. All 86 patched-source cases pass. The two baseline crashes are intentionally retained as reproductions.
- 14 build/install/CLI integration checks pass, alongside reference-integrity and executable-free assembly checks.
- All 18 reported ST28 raw-read k-mer cases and 15 matched assembly cases change from plain Paratyphi B to the positive variant. All 18 patched microassemblies are positive. The 18 libraries represent 16 EnteroBase uberstrain identifiers, not 18 guaranteed independent isolates.
- Six non-ST28 controls (ST110/ST43/ST127 positive contexts; ST5113/ST86/ST8931 negative-marker variants) retain the expected classification in both modes. ST127 lacks an independently obtained biochemical phenotype record; it is a genotype-context control.
- Complete ERR6806453 (1,520,076 pairs) and SRR30326616 (4,666,195 pairs) reproduce the baseline empty-list exception. Patched profiles are `-:r:1,2` and `-:d:-`, respectively, with unresolved O.
- Complete SRR9326726 runs as Havana in baseline and patched workflows. SRR17736741 runs as Thompson in both baseline and patched allele mode with modern SPAdes; the reporter's trimmed input/older-tool hang was not reproduced.
- Numeric/CASAVA header transforms of all ERR6806453 reads preserve sequence/quality and leave the compared reports unchanged. NCBI denies access to SRR2200041 as controlled-access data; this is a surrogate comparison, not the exact SeqSero2S #3 reproduction.
- Native osx-arm64 Conda dependency dry solving succeeds. No Conda package was built or published and no other Python version was executed locally.

## Remaining limits

The supplied issue #56 BLAST fragments select e,n,z15 in the current decision code; the original raw-input misclassification is unverified. The O2:g,m:-/Enteritidis compatibility mapping in #52 is reproduced but unchanged. New subspecies references (#60), private inputs and mixture interpretations require further data/curation. Keep those issues open. The external exact Biopython/bedtools recipe pins (#35/#59/#61) require a separate Bioconda update and rebuilt package.

Prepared source fixes should only be treated as released after merge/release and deployment verification. A passing workflow does not establish correct phenotype or universal serotype accuracy.
