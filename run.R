# run.R — install dependencies (CRAN + Bioconductor) and launch the Shiny app
#
# Two deliberate differences from the version used for the manuscript analysis. Neither can
# affect any analysis result; both are documented in KNOWN_ISSUES.md (H12, H13).
#
#   1. The app is served on 127.0.0.1 (this machine only) instead of 0.0.0.0. The original
#      bound every network interface, which made the app reachable from other machines on
#      the same network — at odds with the app's own guarantee that no data leaves the
#      computer. If you deliberately want to reach it from another machine on a trusted
#      network, change `host` below and be aware of what that exposes.
#   2. Bioconductor installation no longer passes `update = TRUE`, so running this script
#      installs what is missing without silently upgrading unrelated Bioconductor packages
#      in your library. Installation is also confirmed with you first in an interactive
#      session.
#
# Package versions matter here: edgeR's filtering and dispersion behaviour has changed
# across releases. See README.md on renv.lock before relying on results across machines.

options(bitmapType = "cairo")
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "shiny","tidyverse","janitor","readr","ggrepel","ggplot2",
  "readxl","scales","zip","DT","shinyFeedback"
)
bioc_pkgs <- c("edgeR")

installed <- rownames(installed.packages())
to_install_cran <- setdiff(cran_pkgs, installed)
to_install_bioc <- setdiff(bioc_pkgs, installed)

if (length(to_install_cran) || length(to_install_bioc)) {
  message("Missing packages: ",
          paste(c(to_install_cran, to_install_bioc), collapse = ", "))

  proceed <- TRUE
  if (interactive()) {
    ans <- readline("Install these now? Existing packages will not be upgraded. [y/N] ")
    proceed <- tolower(trimws(ans)) %in% c("y", "yes")
  }
  if (!proceed) {
    stop("Dependencies not installed. Install them yourself, or re-run and answer 'y'.",
         call. = FALSE)
  }

  if (length(to_install_cran)) {
    message("Installing from CRAN: ", paste(to_install_cran, collapse = ", "))
    install.packages(to_install_cran, dependencies = TRUE)
  }
  if (length(to_install_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", dependencies = TRUE)
    }
    message("Installing from Bioconductor: ", paste(to_install_bioc, collapse = ", "))
    BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
  }
} else {
  message("All dependencies present.")
}

missing_after <- setdiff(c(cran_pkgs, bioc_pkgs), rownames(installed.packages()))
if (length(missing_after)) {
  stop(paste("These packages failed to install:", paste(missing_after, collapse = ", ")))
}

message("Launching Shiny app on http://127.0.0.1:3838 ...")
shiny::runApp(appDir = ".", host = "127.0.0.1", port = 3838, launch.browser = TRUE)
