# Synthetic example

Nothing here is experimental data. These files exist so that a new user can confirm an
installation behaves correctly, and so that the expected input format is documented by example
rather than by prose.

## Contents

| Path | What it is |
|---|---|
| `counts/S0*_*_results.txt` | six peptide count tables, written by phader's own output routine so the layout is byte-faithful |
| `reference_peptides.xlsx` | the reference list as an Excel sheet (`peptides`, `brand_name`) |
| `reference_peptides.txt` | the same list as pasteable text, for the "Paste peptide + brand pairs" option |
| `expected_figures/*.png` | the four figures this dataset produces at default settings |

The six samples are named so that the default keywords assign them correctly — two `brain`
(target), two `liver` (control), two `input`. Each file holds 11 peptides, of which 8 are in
the reference list and 3 are not, so the reference join is exercised too.

## Planted effects

| Peptide | Name | Design |
|---|---|---|
| CSPWNRGVC | Probe-A | brain-enriched, strong |
| CKEFTYHAC | Probe-B | brain-enriched, moderate |
| CQQWNRGVC | Probe-C | brain-depleted (liver-enriched) |
| CPPWNRGVC | Probe-D | no change |
| CGGGGGGGC | Probe-E | no change, high abundance |
| CAAKTSRVC | Probe-F | present in the input library only |
| CWWDTYPQC | Probe-G | very low in every sample |
| CMMRSKLTC | Probe-H | brain-depleted, strong |

## Running it

Select all six files under `counts/` as the count files, `reference_peptides.xlsx` as the
reference sheet, leave every option at its default, and press **Run analysis**.

| Peptide | Expected outcome |
|---|---|
| Probe-A | logFC +4.32, FDR 5.1e-14 |
| Probe-B | logFC +2.10, FDR 2.7e-11 |
| Probe-C | logFC -4.22, FDR 5.1e-14 |
| Probe-H | logFC -4.25, FDR 5.1e-14 |
| Probe-D | logFC +0.07, FDR 0.41 — not significant |
| Probe-E | logFC +0.06, FDR 0.36 — not significant |
| Probe-F | removed by `filterByExpr` |
| Probe-G | removed by `filterByExpr` |

Six of eight peptides are kept. Probe-F is removed despite being abundant in the input library,
because the edgeR test uses only the target and control columns — a useful reminder that input
samples do not participate in the differential test.

In the ratio plots, Probe-A ranks highest for brain/input (about 4.2) and Probe-H and Probe-C
highest for liver/input (about 4.5 and 3.7). Probe-D and Probe-E sit at approximately 1.0.

## What the figures demonstrate

`expected_figures/volcano_FDR.png` is the clearest illustration of
[KNOWN_ISSUES.md](../KNOWN_ISSUES.md) H1. Probe-C and Probe-H are the most significant points
in the plot, and both are drawn grey and unlabelled as "Not significant", because the
classifier only recognises enrichment. Compare against the table above, and against the
results CSV, which reports them correctly.

Each count file also produces one `readr` parsing warning on load (H11) — expected, caused by
the trailing blank line in phader's output, and harmless.
