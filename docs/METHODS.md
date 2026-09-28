# Statistical methods

What diffhoming v1.0.0 computes, in the order it computes it, at the level of detail a Methods
section needs. Numeric values quoted here are the defaults in force unless the interface says
otherwise; the parameters a user can change are marked. Behaviour that is imperfect is listed
in [KNOWN_ISSUES.md](../KNOWN_ISSUES.md) and cross-referenced here.

## 1. Assembling the count matrix

Each input file is a peptide count table for one sequenced sample. Files are read
positionally (five header lines skipped), the peptide column is detected by name, and each
table is reduced to the peptides present in the user's reference list by an inner join.

The per-sample tables are then combined by successive full outer joins on the peptide
sequence, so the matrix spans the union of peptides observed in any sample, and peptide-sample
combinations that were never observed are set to zero. Display names are attached from the
reference list by a second inner join, which is also what keys the results (H3).

Counts are rounded to integers before analysis (user option, on by default), since edgeR's
negative-binomial model is defined on counts.

Samples are assigned to the target, control and input groups by substring matching on the
count column names, which derive from the original sequencing filenames. Only numeric columns
are eligible. With *strict counts mode* on (default), the run aborts if any selected column
name suggests already-normalised values (`cpm`, `log`, `ratio`, `rpm`, `tpm`), which prevents
transformed values reaching a count-based model. Group membership is not checked for overlap
(H5).

Before modelling, the matrix is checked for negative, missing, infinite and non-integer
values. Negative counts abort the run unless the user explicitly enables clipping; missing and
infinite values always abort. All four counts are reported in the diagnostics output.

## 2. Differential abundance: target versus control

The test compares the **target and control tissues only. Input samples take no part in it**
— they are used solely for the abundance ratios of section 3. A peptide abundant in the input
but absent from both tissues is therefore removed by filtering rather than reported as
depleted.

Counts for the control and target columns are assembled into a `DGEList`, in that order, with
a two-level group factor whose reference level is the control tissue. The contrast tested is
the second coefficient of `~ group`, so a **positive log2 fold change means higher abundance
in the target tissue**.

Filtering proceeds in two stages:

1. Peptides with zero counts across all target and control samples are removed, and library
   sizes are recomputed.
2. `edgeR::filterByExpr` is applied with its defaults — a peptide must reach roughly 10 counts
   in a sufficient number of samples (`min.count = 10`) and 15 counts in total
   (`min.total.count = 15`), evaluated with respect to the group structure. These thresholds
   were developed for RNA-seq libraries and are used here unexamined (H4). The numbers of
   peptides before and after filtering are reported in the app's status log.

Library sizes are normalised by the trimmed mean of M-values (TMM, `calcNormFactors`), which
corrects for differences in sequencing depth and in the degree to which a few dominant
peptides consume the library.

Dispersion estimation then takes one of two paths, determined automatically by the number of
samples per group:

- **With at least two samples in every group**, dispersion is estimated from the data
  (`estimateDisp`) and the test is edgeR's quasi-likelihood F-test (`glmQLFit` followed by
  `glmQLFTest`). This is the recommended path and the one to use where the design allows.
- **With fewer than two samples in a group**, dispersion cannot be estimated. The analysis
  falls back to `glmFit` with a **fixed assumed dispersion of 0.1** (equivalent to a
  biological coefficient of variation of about 32%) and a likelihood-ratio test (`glmLRT`).
  The value 0.1 is a convention for well-controlled experiments; it is an assumption, not a
  measurement, and results from this path should be treated as exploratory. The app does not
  record which path was taken or what dispersion was assumed (H2), so a design without
  replicates must be reported as such in the Methods section by hand.

P-values are adjusted for multiple testing across all peptides that survived filtering by the
Benjamini-Hochberg procedure, reported as `FDR`.

## 3. Abundance relative to the input library

This is a descriptive summary, computed independently of the model above, and it is what the
bar plots show.

Counts for the target, control and input columns are converted to counts per million after TMM
normalisation (`cpm(normalized.lib.sizes = TRUE)`); raw counts can be used instead as a user
option, in which case no depth correction is applied. For each peptide the mean and standard
deviation are taken across the replicates of each group, with the standard deviation defined
as zero where a group has one sample.

Ratios are formed with a pseudocount of one on both sides:

```
target_over_input  = (target_mean  + 1) / (input_mean + 1)
control_over_input = (control_mean + 1) / (input_mean + 1)
```

The pseudocount keeps the ratio finite for peptides absent from the input, at the cost of
shrinking ratios for low-abundance peptides toward 1.

Uncertainty is propagated as the quadrature sum of the relative standard deviations of
numerator and denominator:

```
sd(ratio) = ratio * sqrt( (sd_tissue / (mean_tissue + 1))^2 + (sd_input / (mean_input + 1))^2 )
```

This treats tissue and input measurements as independent and assumes the relative errors are
small; it is a first-order approximation, and where a group has one replicate the
corresponding term is zero, so the interval understates the true uncertainty (H7). Non-finite
results are set to zero. These are descriptive error bars, not confidence intervals, and no
test is applied to the ratios.

## 4. Figures

Both volcano plots show log2 fold change on the x axis, against log2 CPM in one and
-log10(FDR) in the other. Guide lines are drawn at the user's fold-change, CPM and FDR
cutoffs. Points are classified as

```
FDR <= fdr_cut AND logFC > lfc_cut   ->  "<target>-enriched"
everything else                      ->  "Not significant"
```

and only the first group is coloured and labelled. **Peptides that are significantly depleted
in the target tissue fall into the second group**, so they appear grey and unlabelled
regardless of how small their FDR is (H1). The underlying statistics in the results CSV are
unaffected. Axis limits, when set, are applied with `coord_cartesian`, which hides points
outside the range without removing them from the statistics. The shaded region on the logCPM
volcano is drawn to partly fixed coordinates and may not correspond to the plotted range (H6).

The bar plots show the top N peptides by target/input or control/input ratio with the error
bars of section 3, lower whiskers clipped at zero.

## 5. What is not modelled

- No covariates, batch terms or paired designs: the model is a single two-group comparison.
- No comparison across more than two tissues in one run; additional tissues require separate
  runs, and the multiple-testing correction is then per run, not across all of them.
- The input samples are not modelled statistically, only summarised.
- Peptide length, composition and sequence similarity are ignored; peptides are treated as
  independent features, so families of related sequences are not collapsed.
