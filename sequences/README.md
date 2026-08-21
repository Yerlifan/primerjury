# Input sequences

Put your `.fasta` / `.fastq` (optionally gzipped) files here, then run the
pipeline from the repository root.

```
sequences/
├── sample-A1/
│   └── barcode01.fastq.gz
└── sample-B2/
    └── barcode02.fastq.gz
```

Notes:

- Both flat files and one-directory-per-sample layouts are accepted.
- **Base-name collisions are refused, not resolved.** `barcode03.fastq` and
  `barcode03.fastq.gz` would map to the same output name; in the original study
  exactly this caused a silent double-count, so the pipeline now stops instead.
- Nothing in this directory is committed to git except this file.
