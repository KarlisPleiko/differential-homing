# diffhoming

**Differential homing analysis for phage display peptide counts.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R >= 4.2](https://img.shields.io/badge/R-%3E%3D4.2-blue.svg)](https://www.r-project.org/)

diffhoming takes peptide count tables from a biopanning experiment — one per sequenced sample —
and identifies peptides that are differentially abundant in a target tissue relative to a
control tissue, using edgeR's negative-binomial framework. It also reports each peptide's
abundance relative to the injected input library, which is the quantity most often plotted in
a homing experiment.

It is a local Shiny app: data are read from files you select on your own machine, analysed in
your own R session, and written to your own disk. Nothing is uploaded anywhere.

Companion tool: [phader](https://github.com/KarlisPleiko/phader), which produces the count
tables this app consumes.

## Reproducibility notice

**v1.0.0 contains the analysis code used for the differential homing results in the
accompanying manuscript.** `app.R` — every statistic and every plot — is published exactly as
it was run, and the issues found while preparing this release are documented in
[KNOWN_ISSUES.md](KNOWN_ISSUES.md) rather than fixed, so the published figures can be
reproduced from this tag.

`run.R`, the launcher, carries two changes that cannot affect any result: the app is served on
`127.0.0.1` rather than on every network interface, and installing dependencies no longer
upgrades unrelated Bioconductor packages. Both are described in KNOWN_ISSUES (H12, H13).

**Read [KNOWN_ISSUES.md](KNOWN_ISSUES.md) item H1 before interpreting any volcano plot from
this version.** Peptides that are significantly *depleted* in the target tissue are drawn grey
and labelled "Not significant", even when they are the most significant points in the figure.
The results CSV is correct; the figure understates what was found.

## Installation and launch

R 4.2 or newer. From the repository root:

```bash
Rscript run.R
```

`run.R` reports any missing packages, asks before installing them, and then opens the app at
<http://127.0.0.1:3838>. Dependencies are shiny, tidyverse, janitor, readr, ggrepel, ggplot2,
readxl, scales, zip, DT and shinyFeedback from CRAN, plus edgeR from Bioconductor.

To reproduce a specific software environment rather than installing current versions, see
[Dependency versions](#dependency-versions) below.

## Input

### Count tables

One file per sequenced sample, in the output format of
[phader](https://github.com/KarlisPleiko/phader): five header lines, then a column header, then
whitespace-padded rows of `peptide  count`. The app reads them with
`readr::read_table(skip = 5)`, so **the header must be exactly five lines** — the format is
positional, not self-describing.

Each phader file contributes one count column, named after the FASTQ file it came from. The app
joins all files on the peptide column, keeps only peptides present in your reference list, and
replaces missing combinations with zero.

Despite the "Upload CSV files" label, these are normally phader's `.txt` outputs.

### Sample naming is the experimental design

This is the part that surprises people. The app decides which samples are target, control and
input by **searching the count column names for the keywords you type**, and those column
names come from your original FASTQ filenames. Renaming a FASTQ therefore changes the
analysis.

With the default keywords `brain`, `liver` and `input`, a set of files named

```
S01_brain_rep1_R1_001_results.txt     -> target
S02_brain_rep2_R1_001_results.txt     -> target
S03_liver_rep1_R1_001_results.txt     -> control
S04_liver_rep2_R1_001_results.txt     -> control
S05_input_rep1_R1_001_results.txt     -> input
S06_input_rep2_R1_001_results.txt     -> input
```

is assigned as shown. Two rules follow:

1. **Each keyword must appear in exactly the samples that play that role.** A file named
   `brain_input_rep1` matches both `brain` and `input` and silently enters both roles; the app
   does not check that the three groups are disjoint (KNOWN_ISSUES H5).
2. **Check the assignment before trusting a result.** The Status tab lists the columns
   assigned to each role after every run. Read it.

### Reference peptide list

The peptides to analyse, with a display name for each. Either an Excel sheet with columns
`peptides` and `brand_name`, or pasted text with one peptide per line and the name after a tab,
comma, semicolon or two or more spaces. Peptides absent from this list are dropped, so it also
acts as the analysis whitelist.

Results are keyed by `brand_name`, so two peptides sharing a name will collide
(KNOWN_ISSUES H3).

## Options

| Option | Default | Effect |
|---|---|---|
| Target tissue keyword | `brain` | Substring identifying target samples among the count columns. |
| Control tissue keyword | `liver` | Substring identifying control samples. The edgeR test is target vs control. |
| Input keyword | `input` | Substring identifying the input library samples. **These never enter the edgeR test** — they are used only for the abundance-ratio plots. |
| log2FC cutoff | 1 | Vertical guide lines, and the threshold used to classify points as enriched. |
| log2 CPM cutoff | 11 | Horizontal guide line on the logCPM volcano, and the lower edge of its shaded region. |
| FDR cutoff | 0.05 | Horizontal guide line on the FDR volcano and the significance threshold for classification. |
| Top N for bar plots | 15 | How many peptides to show in each ratio plot. |
| Volcano axis limits | blank (auto) | Applied with `coord_cartesian`, so points outside the range are hidden but never dropped from the statistics. |
| Figure font size | 8 pt | Base size; point labels scale from it. |
| Export DPI | 300 | Applies to the on-screen preview size as well as the exported files. |
| Volcano / bar size (cm) | 6x5 / 7x5 | Export dimensions. |
| Strict counts mode | on | Refuses to run if any selected column name contains `cpm`, `log`, `ratio`, `rpm` or `tpm`. A guard against feeding already-normalised values into edgeR. Leave it on. |
| Clip negative values to 0 | off | Troubleshooting only. With it off, negative counts abort the run, which is the safer behaviour. |
| Round numeric columns | on | edgeR expects integer counts. |
| Use CPM for Top ratios | on | Computes the ratio plots from CPM-normalised values rather than raw counts, so samples of different depth are comparable. |

## Outputs

Everything is written to a timestamped directory under R's temporary folder and offered as a
single ZIP. **Download the ZIP before you close the app** — the temporary directory is removed
when the R session ends (KNOWN_ISSUES H16).

| File | Contents |
|---|---|
| `<target>_vs_<control>_results.csv` | edgeR table: `genes`, `logFC`, `logCPM`, `PValue`, `FDR`, one row per peptide that passed filtering |
| `organ_input_stats_<cpm\|counts>.csv` | per-peptide means, standard deviations, and target/input and control/input ratios with propagated uncertainty |
| `merged_data_filt_round.csv` | the joined count matrix that went into the analysis |
| `merged_data_filt_round_cpm.csv` | the same matrix, CPM-normalised |
| `<target>_vs_<control>_volcano_logCPM.pdf` / `.png` | fold change against abundance |
| `<target>_vs_<control>_volcano_FDR.pdf` / `.png` | fold change against significance |
| `<target>_over_input_top<N>_<base>.pdf` | target/input ratio for the top N peptides, with error bars |
| `<control>_over_input_top<N>_<base>.pdf` | the same for the control tissue |
| `diagnostics_files.csv` | per input file: detected peptide column, rows read, rows surviving the reference join |
| `diagnostics_columns.csv` | per selected column: minimum, NA count, Inf count, non-integer count |
| `diagnostics_count_matrix.csv` | counts of negative, NA, Inf and non-integer values in the matrix given to edgeR |
| `diagnostics_eligible_columns.csv` | every column, whether it was eligible as a count column, and why not |

Bar plots are exported as PDF only, while volcanoes are exported as both PDF and PNG
(KNOWN_ISSUES H15).

## Worked example

[`examples/`](examples/) contains six synthetic count files in phader's exact format — two
brain, two liver and two input replicates — with planted effects, plus a matching reference
list in both accepted forms. Nothing in it is experimental data.

Load `examples/counts/*.txt` as the count files and `examples/reference_peptides.xlsx` as the
reference, keep the default keywords, and press **Run analysis**. Expected results:

| Peptide | Planted | Recovered |
|---|---|---|
| Probe-A | brain-enriched, strong | logFC +4.32, FDR 5.1e-14 |
| Probe-B | brain-enriched, moderate | logFC +2.10, FDR 2.7e-11 |
| Probe-C | brain-depleted | logFC -4.22, FDR 5.1e-14 — **drawn as "Not significant"** |
| Probe-H | brain-depleted, strong | logFC -4.25, FDR 5.1e-14 — **drawn as "Not significant"** |
| Probe-D, Probe-E | no change | FDR 0.41, 0.36 |
| Probe-F | input only | removed by `filterByExpr` |
| Probe-G | very low everywhere | removed by `filterByExpr` |

`examples/expected_figures/` holds the four figures this produces. The FDR volcano is the
clearest illustration of H1: the most significant point in the plot is grey and unlabelled.

## Statistics

[docs/METHODS.md](docs/METHODS.md) describes the pipeline in the detail a Methods section
needs: how the count matrix is assembled and filtered, which edgeR path is taken with and
without replicates, the dispersion assumed when replicates are absent, and how the abundance
ratios and their error bars are defined.

## Dependency versions

edgeR's filtering and dispersion estimation have changed across releases, so the package
version is part of the method. Two lockfiles serve different purposes:

- **`renv.lock`** — the environment this release was validated in: R 4.5.3, edgeR 4.8.2,
  ggplot2 4.0.3, shiny 1.14.0 and the rest as recorded. It pins the twelve packages the app
  declares, not their full recursive dependency tree, and it was assembled from the validated
  environment rather than written by `renv::snapshot()`. Treat it as "known to work", not as a
  complete environment specification.
- **`renv.manuscript.lock`** — not present until you generate it. Run
  `Rscript scripts/capture_renv.R` on the machine that produced the published figures; that
  writes a complete lockfile, including every transitive dependency, describing the environment
  the manuscript's results actually came from. This is the file to cite.

Restore either with `renv::restore(lockfile = "<file>")`.

## Limitations

Sixteen items, with effects and workarounds, in [KNOWN_ISSUES.md](KNOWN_ISSUES.md). The ones
that change how you read a result: H1 (depleted peptides labelled not significant), H2 (an
assumed dispersion is used silently when a group has no replicates), H3 (results keyed by
display name, not peptide), H5 (role groups are not checked for overlap).

## How to cite

> Pleiko, K. diffhoming: differential homing analysis for phage display peptide counts
> (v1.0.0). Zenodo. DOI: *to be added on release*

Machine-readable metadata is in [CITATION.cff](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).

## Contact

Questions and bug reports:
[GitHub issues](https://github.com/KarlisPleiko/diffhoming/issues).
