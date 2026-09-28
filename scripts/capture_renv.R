# capture_renv.R — record the exact package versions on THIS machine.
#
# Run this from the repository root on the machine that produced your published figures.
# It writes renv.manuscript.lock, a complete lockfile including the full recursive
# dependency tree, which the shipped renv.lock does not contain (see README.md).
#
#   Rscript scripts/capture_renv.R
#
# Then commit renv.manuscript.lock and cite it in your Methods section. Anyone can restore
# that exact environment with:
#
#   renv::restore(lockfile = "renv.manuscript.lock")

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

out <- "renv.manuscript.lock"

renv::snapshot(
  project  = ".",       # scans app.R and run.R for the packages they load
  lockfile = out,
  prompt   = FALSE,
  force    = TRUE
)

if (file.exists(out)) {
  message("Wrote ", out, " — commit this file alongside the release.")
} else {
  stop("Snapshot did not produce ", out,
       ". Check that renv installed correctly and that you ran this from the repo root.")
}
