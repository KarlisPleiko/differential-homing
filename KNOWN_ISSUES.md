# Known issues and limitations (v1.0.0)

v1.0.0 contains the analysis code used for the differential homing results in the accompanying
manuscript. `app.R` is therefore **deliberately frozen**: the items below are documented rather
than fixed, so that the published figures can be reproduced exactly from this tag. They are the
work list for the next release.

Each entry states the user-visible effect and, where one exists, the workaround. Items marked
*confirmed* were reproduced on the synthetic dataset in `examples/`.

---

## Read this one first

**H1. Significantly depleted peptides are labelled "Not significant".** *(confirmed)*

Both volcano plots classify points as

```
FDR <= fdr_cut AND logFC > lfc_cut  ->  "<target>-enriched"
everything else                     ->  "Not significant"
```

A peptide that is significantly **depleted** in the target tissue — passing the same FDR
threshold with a negative fold change — therefore falls into the second category and is drawn
grey, unlabelled, as "Not significant".

On the example dataset, Probe-C and Probe-H carry FDR = 5.1e-14 at logFC ~ -4.2. In the FDR
volcano they are the most significant points in the figure, and both are grey. See
`examples/expected_figures/volcano_FDR.png`.

**The statistics are not affected** — `<target>_vs_<control>_results.csv` reports the correct
logFC and FDR for every peptide. Only the figures misstate them.

**Workaround.** Read significance from the results CSV, not from the plot colours. If a figure
from this version is used in a publication, either state that only enrichment is highlighted,
or filter `FDR <= 0.05 & logFC < -1` from the CSV and annotate the depleted peptides manually.

---

## Statistical items

**H2. An assumed dispersion is used silently when a group lacks replicates.** *(confirmed)*

With fewer than two samples in either group, dispersion cannot be estimated, and the analysis
falls back to `glmFit` with a fixed dispersion of **0.1** and a likelihood-ratio test. Nothing
in any output file records that the fallback was taken or what value was assumed; the flag is
computed internally and discarded.

On the example data, reducing to one sample per group leaves fold changes almost unchanged
while FDRs weaken by five to six orders of magnitude (5.1e-14 to 2.0e-08 for Probe-A).

**Workaround.** If any group has one sample, say so explicitly in the Methods section, together
with the assumed dispersion of 0.1. See [docs/METHODS.md](docs/METHODS.md) section 2.

**H3. Results are keyed by display name, not by peptide sequence.**

The edgeR table is built with `brand_name` as the feature identifier. Two reference entries
sharing a display name collide irrecoverably in the results. The app warns about duplicated
peptide *sequences* in the reference list but not about duplicated *names*.

**Workaround.** Keep display names unique in the reference sheet. The merged count tables
retain the peptide sequences, so an affected run can be diagnosed from them.

**H4. `filterByExpr` is used with its RNA-seq defaults.**

Filtering requires roughly 10 counts in enough samples and 15 counts in total. These thresholds
were developed for transcript counts, where library composition differs substantially from a
peptide library dominated by a few high-abundance clones. They are used here unexamined.

**Workaround.** The status log reports how many peptides were kept out of how many. Check it:
if a large fraction is being discarded, the threshold rather than the biology may be
responsible.

**H5. Group assignment is not checked for overlap.** 

Target, control and input columns are each matched by independent substring search. A sample
named `brain_input` matches both a `brain` target keyword and an `input` keyword and silently
enters both roles. The app validates that each role matched at least one column, and that no
ineligible column was selected, but never that the three sets are disjoint.

**Workaround.** Choose keywords that cannot co-occur, and read the assignment printed on the
Status tab after every run.

**H6. The shaded region on the logCPM volcano is drawn to partly fixed coordinates.**

Its lower-left corner follows the user's cutoffs, but the upper-right corner is hardcoded. On
data whose axes fall outside those bounds the shading is misleading or invisible, independently
of the thresholds set.

**H7. Ratio error bars rest on undocumented assumptions.**

Uncertainty on the target/input and control/input ratios is propagated in quadrature, treating
tissue and input measurements as independent, with a pseudocount of one added to both
numerator and denominator. Where a group has a single replicate its contribution is zero, so
the interval understates the true uncertainty. These are descriptive error bars; no test is
applied to the ratios.

---

## Input handling

**H8. The input format is assumed rather than checked.** *(confirmed)*

Count files are read with `read_table(skip = 5)`. The five-line skip is hardcoded and matches
phader's output exactly. A file that is not phader output, or phader output whose header
length differs, parses into nonsense rather than raising an error — a data row silently
becomes the column header.

**Workaround.** Use count files from phader v1.0.0. After any run, check
`diagnostics_files.csv`: the detected peptide column should be `peptides` and the row counts
should look sensible.

**H9. `select(-starts_with("x"))` can remove a real sample.**

This line exists to drop the all-`NA` column that phader's trailing tab produces. After
`clean_names()`, any sample whose filename begins with `x` is removed along with it.

**Workaround.** Do not name samples with a leading `x`.

**H10. Diagnostics embed printed output as text.**

`diagnostics_files.csv` stores a captured print of each file's first rows in a single cell.
It is human-readable but not parseable, and it inflates the file.

**H11. Every phader count file produces one parsing warning.** *(confirmed)*

phader's output ends with a trailing blank line, which `read_table` reports as a row with one
column instead of three. The spurious row is then removed by the reference join, so it has
never affected a result — but it means a clean run still prints one warning per input file,
which trains users to ignore warnings.

---

## Changed in this release

These two were defects in the launcher, not in the analysis, so they are fixed in v1.0.0 rather
than documented. Neither can affect any result.

**H12. `run.R` used to bind every network interface.**

It called `runApp(host = "0.0.0.0")`, making the app reachable from any machine on the same
network — at odds with the app's own guarantee that no data leaves the computer. v1.0.0 binds
`127.0.0.1`. To deliberately expose it on a trusted network, change `host` in `run.R`.

**H13. `run.R` used to upgrade unrelated Bioconductor packages.**

It called `BiocManager::install(update = TRUE)`, which silently upgraded other Bioconductor
packages in the user's library. v1.0.0 passes `update = FALSE` and asks for confirmation before
installing anything in an interactive session.

---

## Partially addressed

**H14. Dependency versions are only partly pinned.**

edgeR's filtering and dispersion behaviour has changed across releases, so the package version
is part of the method. `renv.lock` records the twelve declared packages at the versions this
release was validated against, but not their full recursive dependency tree, and it was
assembled from the validated environment rather than produced by `renv::snapshot()`.

**Workaround.** Run `Rscript scripts/capture_renv.R` on the machine that produced your
published figures to write a complete `renv.manuscript.lock`, and cite that file.

---

## Output handling

**H15. Bar plots are exported as PDF only.**

The two volcanoes are written as both PDF and PNG; the two ratio bar plots are written as PDF
only. Inconsistent, and inconvenient if a raster version is needed.

**Workaround.** Export the PNG from the PDF, or screenshot the in-app preview, which is
rendered at the export DPI.

**H16. Outputs live in a temporary directory.**

Results are written to a timestamped folder under R's `tempdir()`, which R removes when the
session ends. The ZIP download is the only copy that persists.

**Workaround.** Download the ZIP before closing the app, and move it somewhere permanent
immediately.
