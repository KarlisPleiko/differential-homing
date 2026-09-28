# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - unreleased

First public release. **`app.R` contains the analysis code used for the differential homing
results in the accompanying manuscript** and is published unchanged, so those figures can be
reproduced from this tag.

### Added

- Shiny interface (`app.R`) running locally: peptide count tables in, differential abundance
  and homing ratios out, with no data leaving the machine.
- edgeR analysis of target versus control tissue: TMM normalisation, `filterByExpr`,
  quasi-likelihood F-test where replicates allow and an assumed-dispersion likelihood-ratio
  test where they do not.
- Abundance relative to the input library, CPM-normalised, with propagated uncertainty.
- Two volcano plots, two ranked-ratio bar plots, and export controls for size, DPI, font size
  and axis limits.
- Four diagnostics tables covering per-file parsing, per-column value ranges, count-matrix
  integrity, and column eligibility.
- Reference peptide list accepted as an Excel sheet or as pasted text.
- Documentation: `README.md`, `docs/METHODS.md`, `KNOWN_ISSUES.md`, `CITATION.cff`.
- `examples/`: six synthetic count files in phader's format with planted effects, a matching
  reference list in both accepted forms, and the four figures they produce.
- `renv.lock` recording the validated environment, and `scripts/capture_renv.R` for capturing
  a complete lockfile from the machine that produced the published figures.

### Changed

Launcher only. Neither change can affect an analysis result; both are documented in
[KNOWN_ISSUES.md](KNOWN_ISSUES.md).

- `run.R` serves the app on `127.0.0.1` instead of `0.0.0.0`, so it is no longer reachable from
  other machines on the network (H12).
- `run.R` installs missing packages without upgrading unrelated Bioconductor packages, and asks
  for confirmation first in an interactive session (H13).

### Known issues

See [KNOWN_ISSUES.md](KNOWN_ISSUES.md), starting with H1: significantly depleted peptides are
drawn as "Not significant" in both volcano plots. Nothing listed there is fixed in this release
by design; fixes are deferred to keep the manuscript's figures reproducible.

## Planned

### [2.0.0]

Coordinated with phader 2.0.0, since the two tools share a file format.

- Adopt an explicit interchange format with phader: TSV with a real header row and a
  `#`-prefixed provenance preamble, read with `read_tsv(comment = "#")`. This removes the
  hardcoded five-line skip (H8) and the parsing warning it produces (H11), and lets the reader
  refuse a file it does not understand.
- Label depletion as well as enrichment in the volcano plots (H1).
- Record the analysis path taken, the dispersion assumed, and the software versions in the
  results files (H2).
- Key results by peptide sequence, with the display name as an attribute (H3).
- Validate that target, control and input groups are disjoint (H5).
- Consume phader's rejection statistics so read-level and peptide-level QC appear in one
  report.
- Write outputs to a user-chosen directory rather than to a temporary one (H16).
