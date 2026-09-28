# app.R — Differential Homing Analysis (Shiny Local App) with DIAGNOSTICS + CPM bars
# ----------------------------------------------------------------------------------
# Runs fully on the user's machine. No data leaves the computer.
# Dependencies: shiny, tidyverse, janitor, readr, ggrepel, ggplot2,
#               readxl, edgeR, scales, zip, DT, shinyFeedback
# Optional: renv for reproducibility

options(bitmapType = "cairo")  # better UTF-8 glyph rendering in PNGs

suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(janitor)
  library(readr)
  library(ggrepel)
  library(ggplot2)
  library(readxl)
  library(edgeR)
  library(scales)
  library(zip)
  library(DT)
  if (requireNamespace("shinyFeedback", quietly = TRUE)) {
    library(shinyFeedback)
  }
})

# ---------- Helpers ----------
msg_ok   <- function(...) { if (requireNamespace("shinyFeedback", quietly=TRUE)) shinyFeedback::showToast("success", paste(...), .options = list(timeOut = 2500)) }
msg_warn <- function(...) { if (requireNamespace("shinyFeedback", quietly = TRUE)) shinyFeedback::showToast("warning", paste(...), .options = list(timeOut = 3500)) }
msg_err  <- function(...) { if (requireNamespace("shinyFeedback", quietly = TRUE)) shinyFeedback::showToast("error",   paste(...), .options = list(timeOut = 5000)) }

.safe_require <- function(pkg) requireNamespace(pkg, quietly = TRUE)

read_peptide_csv <- function(path, skip = 5) {
  suppressMessages(readr::read_table(path, skip = skip, show_col_types = FALSE)) %>% clean_names()
}

get_peptide_col <- function(df) {
  peptide_col <- grep("pept", names(df), value = TRUE, ignore.case = TRUE)
  if (length(peptide_col) == 0) peptide_col <- names(df)[1]
  peptide_col[1]
}

parse_pasted_ref <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
  lines <- unlist(strsplit(txt, "\n", fixed = TRUE))
  rows <- lapply(lines, function(ln) {
    ln <- trimws(ln); if (!nzchar(ln)) return(NULL)
    parts <- strsplit(ln, "\\t|,|;|\\s{2,}")[[1]]
    parts <- trimws(parts[nzchar(parts)])
    if (!length(parts)) return(NULL)
    if (length(parts) == 1) tibble(peptides = parts[1], brand_name = parts[1])
    else tibble(peptides = parts[1], brand_name = paste(parts[-1], collapse = " "))
  })
  df <- bind_rows(rows)
  if (is.null(df) || !nrow(df)) return(NULL)
  df %>% mutate(across(everything(), as.character)) %>%
    filter(!is.na(peptides), nzchar(trimws(peptides))) %>%
    distinct(peptides, .keep_all = TRUE)
}

process_and_filter_file <- function(file_path, ref_peptide_df) {
  df <- read_peptide_csv(file_path)
  pep_col <- get_peptide_col(df)
  df2 <- df %>%
    rename(peptides = !!sym(pep_col)) %>%
    select(peptides, everything()) %>%
    select(-starts_with("x"), -matches("^x\\d+$", ignore.case = TRUE))
  filtered <- inner_join(df2, ref_peptide_df, by = "peptides")
  list(
    basename = basename(file_path),
    peptide_col = pep_col,
    n_rows_in = nrow(df),
    n_rows_after_ref = nrow(filtered),
    head_preview = utils::head(df, 5),
    classes = sapply(df, class),
    filtered = filtered
  )
}

merge_filtered_list <- function(filtered_list) {
  stopifnot(length(filtered_list) > 0)
  dfs <- lapply(filtered_list, function(x) x$filtered)
  merged <- Reduce(function(x, y) full_join(x, y, by = "peptides"), dfs)
  merged %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))
}

min_by_col <- function(df, cols) {
  if (!length(cols)) return(numeric())
  sapply(cols, function(cc) suppressWarnings(min(suppressWarnings(as.numeric(df[[cc]])), na.rm = TRUE)))
}
na_count_by_col <- function(df, cols) {
  if (!length(cols)) return(integer())
  sapply(cols, function(cc) sum(is.na(suppressWarnings(as.numeric(df[[cc]])))))
}
inf_count_by_col <- function(df, cols) {
  if (!length(cols)) return(integer())
  sapply(cols, function(cc) sum(is.infinite(suppressWarnings(as.numeric(df[[cc]])))))
}
non_integer_count_by_col <- function(df, cols) {
  if (!length(cols)) return(integer())
  sapply(cols, function(cc) {
    v <- suppressWarnings(as.numeric(df[[cc]])); v <- v[!is.na(v)]
    sum(abs(v - round(v)) > 1e-8, na.rm = TRUE)
  })
}

eligible_count_cols <- function(df, exclude = c("peptides","brand_name")) {
  all <- setdiff(names(df), exclude)
  all[sapply(df[all], is.numeric)]
}
get_cols <- function(keyword, pool) {
  if (is.null(keyword) || keyword == "") return(character())
  grep(keyword, pool, value = TRUE, ignore.case = TRUE)
}

run_edger <- function(counts_mat, genes_vec, target_cols, control_cols, target_keyword, control_keyword,
                      assumed_dispersion = 0.1, min_count = NULL) {
  group <- factor(
    c(rep(control_keyword, length(control_cols)),
      rep(target_keyword,  length(target_cols))),
    levels = c(control_keyword, target_keyword)
  )
  total_rows_before <- nrow(counts_mat)
  y <- DGEList(counts = counts_mat, genes = data.frame(genes = genes_vec, stringsAsFactors = FALSE), group = group)
  
  keep_nz <- rowSums(y$counts) > 0
  y <- y[keep_nz, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) stop("After removing all-zero rows, no peptides remain for analysis.")
  
  keep <- if (is.null(min_count)) filterByExpr(y, group = group) else filterByExpr(y, group = group, min.count = min_count)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) < 2) stop("Too few peptides after filtering. Relax filtering or check inputs.")
  
  y <- calcNormFactors(y)
  design <- model.matrix(~ group)
  
  if (any(table(group) < 2)) {
    fit <- glmFit(y, design, dispersion = assumed_dispersion)
    tst <- glmLRT(fit, coef = 2)
  } else {
    y  <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    tst <- glmQLFTest(fit, coef = 2)
  }
  tab <- topTags(tst, n = Inf)$table
  if (!"genes" %in% names(tab)) tab$genes <- rownames(tab)
  
  list(
    table = tab,
    kept_after_filter = nrow(y),
    total_rows_before = total_rows_before,
    used_fallback = any(table(group) < 2),
    n_control = sum(group == control_keyword),
    n_target  = sum(group == target_keyword)
  )
}

# ---------- Plot builders with base_size ----------
.label_size_from_base <- function(base_size) {
  max(1.4, base_size * 0.28)
}

make_volcano_logcpm <- function(tab, target_keyword, control_keyword,
                                lfc_cut = 1, cpm_cut = 11, base_size = 8,
                                xlim = NULL, ylim = NULL) {
  lab_sz <- .label_size_from_base(base_size)
  tab2 <- tab %>% mutate(
    Differential_homing = case_when(
      FDR <= 0.05 & logFC > lfc_cut ~ paste0(target_keyword, "-enriched"),
      TRUE ~ "Not significant"
    ),
    label = ifelse(Differential_homing != "Not significant", genes, NA)
  )
  fill_colors <- setNames(c("#D73027", "grey80"),
                          c(paste0(target_keyword, "-enriched"), "Not significant"))
  
  ggplot(tab2, aes(x = logFC, y = logCPM)) +
    annotate("rect", xmin = lfc_cut, xmax = 10, ymin = cpm_cut, ymax = 23, alpha = 0.12, fill = "lightgreen") +
    geom_point(aes(fill = Differential_homing), shape = 21, size = 2.5, stroke = 0.5, color = "black") +
    geom_text_repel(
      data = tab2 %>% filter(Differential_homing != "Not significant"),
      aes(label = genes), size = lab_sz, max.overlaps = Inf, force = 8,
      box.padding = 0.5, point.padding = 0.5, segment.color = "black", segment.size = 0.3
    ) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = cpm_cut, linetype = "dashed", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = fill_colors) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", size = base_size),
      axis.text  = element_text(color = "black", size = base_size - 1),
      legend.position = "none"
    ) +
    labs(x = paste0("log2 Fold Change (", target_keyword, " vs ", control_keyword, ")"),
         y = "log2 CPM")
}

make_volcano_fdr <- function(tab, target_keyword, control_keyword,
                             lfc_cut = 1, fdr_cut = 0.05, base_size = 8,
                             xlim = NULL, ylim = NULL) {
  lab_sz <- .label_size_from_base(base_size)
  tab2 <- tab %>% mutate(
    Differential_homing = case_when(
      FDR <= fdr_cut & logFC > lfc_cut ~ paste0(target_keyword, "-enriched"),
      TRUE ~ "Not significant"
    ),
    label = ifelse(Differential_homing != "Not significant", genes, NA)
  )
  fill_colors <- setNames(c("#D73027","grey80"),
                          c(paste0(target_keyword, "-enriched"), "Not significant"))
  
  ggplot(tab2, aes(x = logFC, y = -log10(FDR))) +
    annotate("rect", xmin = lfc_cut, xmax = 10, ymin = -log10(fdr_cut), ymax = Inf, alpha = 0.12, fill = "lightgreen") +
    geom_point(aes(fill = Differential_homing), shape = 21, size = 2.5, stroke = 0.5, color = "black") +
    geom_text_repel(
      data = tab2 %>% filter(Differential_homing != "Not significant"),
      aes(label = genes), size = lab_sz, max.overlaps = Inf, force = 8,
      box.padding = 0.5, point.padding = 0.5, segment.color = "black", segment.size = 0.3
    ) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = fill_colors) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", size = base_size),
      axis.text  = element_text(color = "black", size = base_size - 1),
      legend.position = "none"
    ) +
    labs(x = paste0("log2 Fold Change (", target_keyword, " vs ", control_keyword, ")"),
         y = expression(-log[10]("FDR")))
}

compute_stats <- function(merged_df_round, target_cols, control_cols, input_cols) {
  safe_sd <- function(x) { x <- x[!is.na(x)]; if (length(x) <= 1) 0 else stats::sd(x) }
  target_df  <- merged_df_round %>% select(peptides, all_of(target_cols))
  control_df <- merged_df_round %>% select(peptides, all_of(control_cols))
  input_df   <- merged_df_round %>% select(peptides, all_of(input_cols))
  tibble(
    peptides     = merged_df_round$peptides,
    brand_name   = merged_df_round$brand_name,
    target_mean  = rowMeans(target_df[,-1, drop = FALSE],  na.rm = TRUE),
    target_sd    = apply(target_df[,-1, drop = FALSE],  1, safe_sd),
    control_mean = rowMeans(control_df[,-1, drop = FALSE], na.rm = TRUE),
    control_sd   = apply(control_df[,-1, drop = FALSE], 1, safe_sd),
    input_mean   = rowMeans(input_df[,-1, drop = FALSE],   na.rm = TRUE),
    input_sd     = apply(input_df[,-1, drop = FALSE],   1, safe_sd)
  ) %>%
    mutate(
      target_over_input  = (target_mean  + 1) / (input_mean + 1),
      control_over_input = (control_mean + 1) / (input_mean + 1),
      target_ratio_sd  = target_over_input  * sqrt((target_sd  / (target_mean  + 1))^2 + (input_sd / (input_mean + 1))^2),
      control_ratio_sd = control_over_input * sqrt((control_sd / (control_mean + 1))^2 + (input_sd / (input_mean + 1))^2)
    ) %>%
    mutate(across(c(target_ratio_sd, control_ratio_sd,
                    target_mean, control_mean, input_mean,
                    target_sd, control_sd, input_sd),
                  ~ dplyr::if_else(is.finite(.x), .x, 0)))
}

plot_topN <- function(df, var, var_sd, fill_color, label, topN = 15, base_size = 8) {
  df %>%
    arrange(desc({{ var }})) %>%
    slice_head(n = topN) %>%
    mutate(
      brand_name = factor(brand_name, levels = brand_name),
      y  = {{ var }},
      sd = tidyr::replace_na({{ var_sd }}, 0)
    ) %>%
    ggplot(aes(x = reorder(brand_name, -y), y = y)) +
    geom_col(fill = fill_color, width = 0.7) +
    geom_errorbar(aes(ymin = pmax(y - sd, 0), ymax = y + sd),
                  width = 0.25, linewidth = 0.3) +
    labs(x = NULL, y = label) +
    theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", size = base_size),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = base_size - 1),
      axis.text.y = element_text(color = "black", size = base_size - 1),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      plot.margin = margin(3, 3, 3, 3, "mm")
    )
}

# ---------------- UI ----------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .note {color:#555;} .ok{color:#2a7;} .warn{color:#b50;} .err{color:#b22;}
    .small {font-size: 12px;}
  "))),
  titlePanel("Differential Homing Analysis (edgeR, local)"),
  sidebarLayout(
    sidebarPanel(width = 4,
                 h4("Inputs"),
                 fileInput("csv_files", "Upload CSV files", multiple = TRUE, accept = c(".csv", ".tsv", ".txt")),
                 radioButtons("ref_source", "Reference list source",
                              choices = c("Upload Excel (.xlsx/.xls)" = "upload",
                                          "Paste peptide + brand pairs" = "paste"),
                              selected = "upload"),
                 conditionalPanel(
                   condition = "input.ref_source == 'upload'",
                   fileInput("ref_file", "Upload reference sheet (Excel)", accept = c(".xlsx", ".xls"))
                 ),
                 conditionalPanel(
                   condition = "input.ref_source == 'paste'",
                   helpText("Paste one item per line. Accepts tab/comma/semicolon/≥2 spaces between peptide and brand.",
                            "If only the peptide is given, it will be used as the brand name too."),
                   tags$details(tags$summary("Show example"),
                                tags$pre("CAGTGGT, BrainBinder\nRGDXXXX\tIntegrin-Probe\nVTVPPTN  MyBrand\nGFLGGGG")),
                   textAreaInput("ref_paste", "Pasted peptide + brand list",
                                 placeholder = "PEPTIDE, BrandName\nPEP2\tBrand2\nPEP3  Brand Three\nPEP4", rows = 8)
                 ),
                 textInput("target_keyword",  "Target tissue keyword",  value = "brain"),
                 textInput("control_keyword", "Control tissue keyword", value = "liver"),
                 textInput("input_keyword",   "Input keyword",          value = "input"),
                 hr(),
                 numericInput("lfc_cut", "log2FC cutoff (vertical lines)", value = 1, step = 0.1),
                 numericInput("cpm_cut", "log2 CPM cutoff (horizontal line)", value = 11, step = 0.5),
                 numericInput("fdr_cut", "FDR cutoff", value = 0.05, step = 0.01, min = 0, max = 1),
                 numericInput("topN", "Top N for bar plots", value = 15, min = 5, max = 50, step = 1),
                 
                 hr(),
                 h4("Volcano axis limits (optional)"),
                 helpText("Leave blank to auto-scale. Limits are applied with coord_cartesian (no points dropped)."),
                 fluidRow(
                   column(6, numericInput("volc_xmin", "x min (logFC)", value = NA, step = 0.5)),
                   column(6, numericInput("volc_xmax", "x max (logFC)", value = NA, step = 0.5))
                 ),
                 fluidRow(
                   column(6, numericInput("volc_ymin", "y min", value = NA, step = 0.5)),
                   column(6, numericInput("volc_ymax", "y max", value = NA, step = 0.5))
                 ),
                 
                 hr(),
                 h4("Figure size, DPI, fonts"),
                 sliderInput("font_base", "Figure font size (pt)", min = 6, max = 18, value = 8, step = 1),
                 numericInput("dpi", "Export DPI", value = 300, min = 72, max = 1200, step = 1),
                 tags$strong("Volcano size (cm)"),
                 fluidRow(
                   column(6, numericInput("volcano_w_cm", "Width",  value = 6, min = 3, max = 40, step = 0.1)),
                   column(6, numericInput("volcano_h_cm", "Height", value = 5, min = 3, max = 40, step = 0.1))
                 ),
                 tags$strong("Bar plot size (cm)"),
                 fluidRow(
                   column(6, numericInput("bar_w_cm", "Width",  value = 7, min = 3, max = 40, step = 0.1)),
                   column(6, numericInput("bar_h_cm", "Height", value = 5, min = 3, max = 40, step = 0.1))
                 ),
                 hr(),
                 checkboxInput("strict_counts", "Strict counts mode (block CPM/log/ratios)", value = TRUE),
                 checkboxInput("clip_neg", "Clip negative values to 0 (for troubleshooting only)", value = FALSE),
                 checkboxInput("round_counts", "Round numeric columns before edgeR", value = TRUE),
                 checkboxInput("bar_use_cpm", "Use CPM for Top ratios (recommended)", value = TRUE),
                 actionButton("run_btn", "Run analysis", class = "btn-primary"),
                 br(), br(),
                 downloadButton("download_zip", "Download all outputs (ZIP)")
    ),
    mainPanel(width = 8,
              tabsetPanel(id = "tabs",
                          tabPanel("Status", verbatimTextOutput("status"), hr(), uiOutput("detected_cols")),
                          tabPanel("Diagnostics",
                                   h5("Per-file summary"), DTOutput("diag_files"), br(),
                                   h5("Per-column summary (selected target/control/input)"), DTOutput("diag_cols"), br(),
                                   h5("Count matrix issues (what would go into edgeR)"), DTOutput("diag_counts"), br(),
                                   h5("Eligible count columns"), DTOutput("diag_eligible"),
                                   p(class="small note", "Tip: If negatives/NA appear, fix the source files or adjust toggles on the left.")
                          ),
                          tabPanel("Volcano (logCPM)", uiOutput("plot_volc_cpm_ui")),
                          tabPanel("Volcano (FDR)",    uiOutput("plot_volc_fdr_ui")),
                          tabPanel("Top ratios",
                                   uiOutput("plot_target_ui"),
                                   br(),
                                   uiOutput("plot_control_ui")
                          ),
                          tabPanel("Tables",
                                   h5("edgeR results"), DTOutput("tbl_edger"), br(),
                                   h5("Organ/input stats"), DTOutput("tbl_stats"))
              )
    )
  )
)

# ---------------- Server ----------------
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1024 * 1024^2)  # 1 GB
  
  cm_to_px <- function(cm, dpi) as.integer(round((as.numeric(cm) / 2.54) * as.numeric(dpi)))
  safe_num <- function(x, default) {
    x2 <- suppressWarnings(as.numeric(x))
    if (length(x2) != 1 || !is.finite(x2) || is.na(x2) || x2 <= 0) default else x2
  }
  safe_lim <- function(x) {
    x2 <- suppressWarnings(as.numeric(x))
    if (length(x2) != 1 || is.na(x2) || !is.finite(x2)) NULL else x2
  }
  
  volc_w_px <- reactive(cm_to_px(safe_num(input$volcano_w_cm, 6), safe_num(input$dpi, 300)))
  volc_h_px <- reactive(cm_to_px(safe_num(input$volcano_h_cm, 5), safe_num(input$dpi, 300)))
  bar_w_px  <- reactive(cm_to_px(safe_num(input$bar_w_cm, 7),  safe_num(input$dpi, 300)))
  bar_h_px  <- reactive(cm_to_px(safe_num(input$bar_h_cm, 5),  safe_num(input$dpi, 300)))
  
  volc_xlim <- reactive({
    a <- safe_lim(input$volc_xmin); b <- safe_lim(input$volc_xmax)
    if (is.null(a) && is.null(b)) return(NULL)
    c(a, b)
  })
  volc_ylim <- reactive({
    a <- safe_lim(input$volc_ymin); b <- safe_lim(input$volc_ymax)
    if (is.null(a) && is.null(b)) return(NULL)
    c(a, b)
  })
  
  rv <- reactiveValues(
    merged_df = NULL, merged_round = NULL, edger_tab = NULL, organ_stats = NULL,
    target_cols = character(), control_cols = character(), input_cols = character(),
    eligible_pool = character(), files_written = character(), outdir = NULL,
    log_msgs = character(), diag_files_df = NULL, diag_cols_df = NULL,
    diag_counts_df = NULL, diag_eligible_df = NULL
  )
  
  log_line <- function(...) {
    rv$log_msgs <- c(rv$log_msgs, paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(...)))
  }
  
  observeEvent(input$run_btn, {
    req(input$csv_files)
    if (identical(input$ref_source, "upload")) req(input$ref_file) else req(input$ref_paste)
    
    rv$log_msgs <- character()
    rv$files_written <- character()
    rv$diag_files_df <- rv$diag_cols_df <- rv$diag_counts_df <- rv$diag_eligible_df <- NULL
    rv$edger_tab <- rv$organ_stats <- NULL
    
    withProgress(message = "Running analysis", value = 0, {
      
      incProgress(0.05, detail = "Preparing reference list")
      
      ref_raw <- if (identical(input$ref_source, "upload")) {
        tryCatch({
          suppressMessages(readxl::read_excel(input$ref_file$datapath)) %>% clean_names()
        }, error = function(e) {
          log_line("Error reading reference sheet: ", e$message)
          showNotification("Failed to read reference sheet", type = "error")
          NULL
        })
      } else parse_pasted_ref(input$ref_paste)
      validate(need(!is.null(ref_raw) && nrow(ref_raw) > 0, "Reference sheet/list could not be read."))
      
      if (!all(c("peptides","brand_name") %in% names(ref_raw))) {
        nm <- tolower(names(ref_raw))
        idx_p <- which(nm %in% c("peptides","peptide","pep","seq","sequence"))
        idx_b <- which(nm %in% c("brand_name","brand","name","label"))
        if (length(idx_p) >= 1) names(ref_raw)[idx_p[1]] <- "peptides"
        if (length(idx_b) >= 1) names(ref_raw)[idx_b[1]] <- "brand_name"
      }
      validate(need(all(c("peptides","brand_name") %in% names(ref_raw)),
                    "Reference must contain 'peptides' and 'brand_name' columns."))
      
      if (anyDuplicated(ref_raw$peptides)) {
        dups <- ref_raw$peptides[duplicated(ref_raw$peptides)] %>% unique() %>% head(5)
        log_line("WARNING: Duplicate peptide IDs in reference (up to 5): ", paste(dups, collapse=", "))
        msg_warn("Reference contains duplicated peptides. Keeping first occurrence.")
        ref_raw <- ref_raw %>% distinct(peptides, .keep_all = TRUE)
      }
      validate(need(!any(is.na(ref_raw$peptides)), "Reference has NA in 'peptides'."))
      
      ref_pep_only <- ref_raw %>% select(peptides)
      ref_full     <- ref_raw %>% select(peptides, brand_name)
      
      incProgress(0.15, detail = "Filtering uploaded CSV files")
      files <- input$csv_files$datapath
      names(files) <- basename(input$csv_files$name)
      
      diag_files <- list()
      filtered_list <- lapply(seq_along(files), function(i) {
        p <- files[i]; nm <- names(files)[i]
        log_line("Processing file: ", nm)
        tryCatch({
          fx <- process_and_filter_file(p, ref_pep_only)
          diag_files[[nm]] <<- tibble(
            file = nm,
            peptide_col_detected = fx$peptide_col,
            rows_in = fx$n_rows_in,
            rows_after_ref_join = fx$n_rows_after_ref,
            n_cols = ncol(fx$filtered),
            first_5_rows_preview = paste(capture.output(print(fx$head_preview)), collapse = " | "),
            classes = paste(paste(names(fx$classes), unlist(fx$classes), sep=":"), collapse = ", ")
          )
          fx
        }, error = function(e) {
          log_line("Error processing ", nm, ": ", e$message); NULL
        })
      })
      filtered_list <- Filter(Negate(is.null), filtered_list)
      validate(need(length(filtered_list) > 0, "No valid CSV files after processing."))
      
      rv$diag_files_df <- bind_rows(diag_files)
      
      incProgress(0.35, detail = "Merging filtered data")
      merged <- merge_filtered_list(filtered_list)
      merged_filt <- inner_join(merged, ref_full, by = "peptides")
      
      if (nrow(merged_filt) < nrow(merged)) {
        log_line("NOTE: Rows reduced after adding 'brand_name' via reference (inner_join).")
      }
      if (anyDuplicated(merged_filt$peptides)) {
        dups <- merged_filt$peptides[duplicated(merged_filt$peptides)] %>% unique() %>% head(5)
        log_line("WARNING: Duplicate peptides after merge (up to 5): ", paste(dups, collapse=", "))
      }
      
      merged_round <- if (isTRUE(input$round_counts)) merged_filt %>% mutate(across(where(is.numeric), round)) else merged_filt
      rv$merged_df <- merged_filt
      rv$merged_round <- merged_round
      
      eligible_pool <- eligible_count_cols(merged_round, exclude = c("peptides","brand_name"))
      rv$eligible_pool <- eligible_pool
      
      all_cols <- names(merged_round)
      classes  <- sapply(merged_round, function(x) class(x)[1])
      is_num   <- sapply(merged_round, is.numeric)
      is_ref   <- all_cols %in% c("peptides","brand_name")
      rv$diag_eligible_df <- tibble(
        column   = all_cols,
        class    = classes,
        numeric  = is_num,
        reserved = is_ref,
        eligible = all_cols %in% eligible_pool,
        reason   = dplyr::case_when(
          is_ref ~ "reserved/ref column",
          !is_num ~ "non-numeric",
          TRUE ~ ""
        )
      )
      
      target_cols  <- get_cols(input$target_keyword, eligible_pool)
      control_cols <- get_cols(input$control_keyword, eligible_pool)
      input_cols   <- get_cols(input$input_keyword,  eligible_pool)
      
      rv$target_cols  <- target_cols
      rv$control_cols <- control_cols
      rv$input_cols   <- input_cols
      
      validate(need(length(target_cols)  > 0, paste0("No columns found for target keyword within eligible pool: ", input$target_keyword)))
      validate(need(length(control_cols) > 0, paste0("No columns found for control keyword within eligible pool: ", input$control_keyword)))
      validate(need(length(input_cols)   > 0, paste0("No columns found for input keyword within eligible pool: ",   input$input_keyword)))
      
      all_sel <- c(target_cols, control_cols, input_cols)
      rv$diag_cols_df <- tibble(
        column = all_sel,
        min_value = as.numeric(min_by_col(merged_round, all_sel)[all_sel]),
        n_na      = as.integer(na_count_by_col(merged_round, all_sel)[all_sel]),
        n_inf     = as.integer(inf_count_by_col(merged_round, all_sel)[all_sel]),
        n_non_integer = as.integer(non_integer_count_by_col(merged_round, all_sel)[all_sel]),
        class     = sapply(all_sel, function(cc) class(merged_round[[cc]])[1])
      )
      
      non_eligible_hits <- setdiff(all_sel, eligible_pool)
      validate(need(!length(non_eligible_hits),
                    paste0("Non-eligible columns selected: ", paste(non_eligible_hits, collapse=", "))))
      if (isTRUE(input$strict_counts)) {
        bad_named <- grep("cpm|log|ratio|rpm|tpm", all_sel, ignore.case = TRUE, value = TRUE)
        validate(need(!length(bad_named),
                      paste0("Strict counts mode: remove transformed columns (CPM/log/ratio): ",
                             paste(bad_named, collapse=", "))))
      }
      
      incProgress(0.50, detail = "Building count matrix & validating")
      count_cols <- c(control_cols, target_cols)
      counts_mat <- as.matrix(merged_round[, count_cols, drop = FALSE])
      mode(counts_mat) <- "numeric"
      
      neg_n   <- sum(counts_mat < 0, na.rm = TRUE)
      na_n    <- sum(is.na(counts_mat))
      inf_n   <- sum(is.infinite(counts_mat))
      nonint_n<- sum(abs(counts_mat - round(counts_mat)) > 1e-8, na.rm = TRUE)
      
      if (neg_n > 0) log_line("Found ", neg_n, " negative values in count columns.")
      if (na_n  > 0) log_line("Found ", na_n,  " NA values in count columns.")
      if (inf_n > 0) log_line("Found ", inf_n, " Inf values in count columns.")
      if (nonint_n > 0 && !isTRUE(input$round_counts)) log_line("Notice: ", nonint_n, " non-integers in counts.")
      
      rv$diag_counts_df <- tibble(issue = c("negatives","NA","Inf","non_integer"),
                                  count = c(neg_n, na_n, inf_n, nonint_n))
      
      if (neg_n > 0 && isTRUE(input$clip_neg)) counts_mat[counts_mat < 0] <- 0
      if (neg_n > 0 && !isTRUE(input$clip_neg)) validate(need(FALSE, "Negative counts detected. Enable clip or fix data."))
      validate(need(na_n  == 0, "NA values detected in counts matrix."))
      validate(need(inf_n == 0, "Infinite values detected in counts matrix."))
      
      incProgress(0.62, detail = "Running edgeR")
      res_edger <- run_edger(
        counts_mat   = counts_mat,
        genes_vec    = rv$merged_round$brand_name,
        target_cols  = target_cols,
        control_cols = control_cols,
        target_keyword  = input$target_keyword,
        control_keyword = input$control_keyword
      )
      tab <- res_edger$table
      rv$edger_tab <- tab
      
      log_line(sprintf("edgeR: kept %d/%d peptides after filtering.",
                       res_edger$kept_after_filter, res_edger$total_rows_before))
      
      stats_cols <- unique(c(target_cols, control_cols, input_cols))
      stats_df   <- rv$merged_round[, c("peptides","brand_name", stats_cols), drop = FALSE]
      
      if (isTRUE(input$bar_use_cpm)) {
        incProgress(0.70, detail = "Computing CPM for Top ratios")
        dge_stats <- DGEList(counts = as.matrix(stats_df[, stats_cols, drop = FALSE]))
        dge_stats <- calcNormFactors(dge_stats)
        cpm_mat   <- cpm(dge_stats, normalized.lib.sizes = TRUE)
        merged_cpm <- stats_df
        merged_cpm[, colnames(cpm_mat)] <- cpm_mat
        stats_base <- "cpm"
        log_line("Top ratios: using CPM.")
      } else {
        merged_cpm <- stats_df
        stats_base <- "counts"
        log_line("Top ratios: using raw counts.")
      }
      
      outdir <- file.path(tempdir(), paste0("diff_homing_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      
      write_csv(rv$diag_files_df,    file.path(outdir, "diagnostics_files.csv"))
      write_csv(rv$diag_cols_df,     file.path(outdir, "diagnostics_columns.csv"))
      write_csv(rv$diag_counts_df,   file.path(outdir, "diagnostics_count_matrix.csv"))
      write_csv(rv$diag_eligible_df, file.path(outdir, "diagnostics_eligible_columns.csv"))
      
      res_file <- file.path(outdir, paste0(input$target_keyword, "_vs_", input$control_keyword, "_results.csv"))
      readr::write_csv(tab, res_file)
      rv$files_written <- c(rv$files_written, res_file)
      
      incProgress(0.74, detail = "Computing ratios and stats")
      organ_stats <- compute_stats(merged_cpm, target_cols, control_cols, input_cols)
      rv$organ_stats <- organ_stats
      
      stats_file <- file.path(outdir, paste0("organ_input_stats_", stats_base, ".csv"))
      readr::write_csv(organ_stats, stats_file)
      rv$files_written <- c(rv$files_written, stats_file)
      
      incProgress(0.85, detail = "Rendering plots")
      base_pt <- safe_num(input$font_base, 8)
      
      p1 <- make_volcano_logcpm(tab, input$target_keyword, input$control_keyword,
                                lfc_cut = input$lfc_cut, cpm_cut = input$cpm_cut, base_size = base_pt,
                                xlim = volc_xlim(), ylim = volc_ylim())
      p2 <- make_volcano_fdr(tab, input$target_keyword, input$control_keyword,
                             lfc_cut = input$lfc_cut, fdr_cut = input$fdr_cut, base_size = base_pt,
                             xlim = volc_xlim(), ylim = volc_ylim())
      
      dpi_exp <- safe_num(input$dpi, 300)
      volc_w  <- safe_num(input$volcano_w_cm, 6)
      volc_h  <- safe_num(input$volcano_h_cm, 5)
      bar_w   <- safe_num(input$bar_w_cm, 7)
      bar_h   <- safe_num(input$bar_h_cm, 5)
      
      volc1_pdf <- file.path(outdir, paste0(input$target_keyword, "_vs_", input$control_keyword, "_volcano_logCPM.pdf"))
      volc1_png <- sub(".pdf$", ".png", volc1_pdf)
      ggsave(volc1_pdf, p1, units = "cm", width = volc_w, height = volc_h, dpi = dpi_exp)
      ggsave(volc1_png, p1, units = "cm", width = volc_w, height = volc_h, dpi = dpi_exp)
      rv$files_written <- c(rv$files_written, volc1_pdf, volc1_png)
      
      volc2_pdf <- file.path(outdir, paste0(input$target_keyword, "_vs_", input$control_keyword, "_volcano_FDR.pdf"))
      volc2_png <- sub(".pdf$", ".png", volc2_pdf)
      ggsave(volc2_pdf, p2, units = "cm", width = volc_w, height = volc_h, dpi = dpi_exp)
      ggsave(volc2_png, p2, units = "cm", width = volc_w, height = volc_h, dpi = dpi_exp)
      rv$files_written <- c(rv$files_written, volc2_pdf, volc2_png)
      
      pt <- plot_topN(organ_stats, target_over_input,  target_ratio_sd,  "#D73027",
                      paste(input$target_keyword,  "/ Input ratio (", stats_base, ")", sep=""),
                      topN = input$topN, base_size = base_pt)
      pc <- plot_topN(organ_stats, control_over_input, control_ratio_sd, "#4575B4",
                      paste(input$control_keyword, "/ Input ratio (", stats_base, ")", sep=""),
                      topN = input$topN, base_size = base_pt)
      
      pt_pdf <- file.path(outdir, paste0(input$target_keyword,  "_over_input_top", input$topN, "_", stats_base, ".pdf"))
      pc_pdf <- file.path(outdir, paste0(input$control_keyword, "_over_input_top", input$topN, "_", stats_base, ".pdf"))
      ggsave(pt_pdf, pt, width = bar_w, height = bar_h, units = "cm", dpi = dpi_exp)
      ggsave(pc_pdf, pc, width = bar_w, height = bar_h, units = "cm", dpi = dpi_exp)
      rv$files_written <- c(rv$files_written, pt_pdf, pc_pdf)
      
      # --- Existing: merged rounded counts table ---
      merged_file <- file.path(outdir, "merged_data_filt_round.csv")
      readr::write_csv(merged_round, merged_file)
      rv$files_written <- c(rv$files_written, merged_file)
      
      # --- NEW: merged CPM table (same shape as merged_data_filt_round) ---
      incProgress(0.90, detail = "Writing merged CPM table")
      cpm_cols <- eligible_count_cols(merged_round, exclude = c("peptides","brand_name"))
      dge_all  <- DGEList(counts = as.matrix(merged_round[, cpm_cols, drop = FALSE]))
      dge_all  <- calcNormFactors(dge_all)
      cpm_all  <- cpm(dge_all, normalized.lib.sizes = TRUE)
      
      merged_round_cpm <- merged_round
      merged_round_cpm[, cpm_cols] <- cpm_all
      
      merged_cpm_file <- file.path(outdir, "merged_data_filt_round_cpm.csv")
      readr::write_csv(merged_round_cpm, merged_cpm_file)
      rv$files_written <- c(rv$files_written, merged_cpm_file)
      
      rv$outdir <- outdir
      log_line("Done. Outputs saved to ", outdir)
      incProgress(1, detail = "Complete")
    })
    
    updateTabsetPanel(session, "tabs", selected = "Diagnostics")
  })
  
  output$status <- renderText({
    if (!length(rv$log_msgs)) return("Idle. Upload files and press 'Run analysis'.")
    paste(rv$log_msgs, collapse = "\n")
  })
  
  output$detected_cols <- renderUI({
    req(rv$merged_round)
    tagList(
      h5("Detected columns"),
      p(HTML(paste0("<b>Target</b>: ", paste(rv$target_cols, collapse = ", ")))),
      p(HTML(paste0("<b>Control</b>: ", paste(rv$control_cols, collapse = ", ")))),
      p(HTML(paste0("<b>Input</b>: ", paste(rv$input_cols, collapse = ", "))))
    )
  })
  
  output$plot_volc_cpm_ui <- renderUI({
    div(style = "overflow:hidden;",
        plotOutput("plot_volc_cpm",
                   width  = paste0(volc_w_px(), "px"),
                   height = paste0(volc_h_px(), "px")))
  })
  output$plot_volc_fdr_ui <- renderUI({
    div(style = "overflow:hidden;",
        plotOutput("plot_volc_fdr",
                   width  = paste0(volc_w_px(), "px"),
                   height = paste0(volc_h_px(), "px")))
  })
  output$plot_target_ui <- renderUI({
    div(style = "overflow:hidden;",
        plotOutput("plot_target",
                   width  = paste0(bar_w_px(), "px"),
                   height = paste0(bar_h_px(), "px")))
  })
  output$plot_control_ui <- renderUI({
    div(style = "overflow:hidden;",
        plotOutput("plot_control",
                   width  = paste0(bar_w_px(), "px"),
                   height = paste0(bar_h_px(), "px")))
  })
  
  output$plot_volc_cpm <- renderPlot({
    req(rv$edger_tab)
    make_volcano_logcpm(rv$edger_tab, input$target_keyword, input$control_keyword,
                        lfc_cut = input$lfc_cut, cpm_cut = input$cpm_cut,
                        base_size = safe_num(input$font_base, 8),
                        xlim = volc_xlim(), ylim = volc_ylim())
  })
  
  output$plot_volc_fdr <- renderPlot({
    req(rv$edger_tab)
    make_volcano_fdr(rv$edger_tab, input$target_keyword, input$control_keyword,
                     lfc_cut = input$lfc_cut, fdr_cut = input$fdr_cut,
                     base_size = safe_num(input$font_base, 8),
                     xlim = volc_xlim(), ylim = volc_ylim())
  })
  
  output$plot_target <- renderPlot({
    req(rv$organ_stats)
    plot_topN(rv$organ_stats, target_over_input,  target_ratio_sd,  "#D73027",
              paste(input$target_keyword,  "/ Input ratio"),
              topN = input$topN, base_size = safe_num(input$font_base, 8))
  })
  
  output$plot_control <- renderPlot({
    req(rv$organ_stats)
    plot_topN(rv$organ_stats, control_over_input, control_ratio_sd, "#4575B4",
              paste(input$control_keyword, "/ Input ratio"),
              topN = input$topN, base_size = safe_num(input$font_base, 8))
  })
  
  output$diag_files     <- renderDT({ req(rv$diag_files_df);    datatable(rv$diag_files_df, options = list(pageLength = 10, scrollX = TRUE)) })
  output$diag_cols      <- renderDT({ req(rv$diag_cols_df);     datatable(rv$diag_cols_df,  options = list(pageLength = 15, scrollX = TRUE)) })
  output$diag_counts    <- renderDT({ req(rv$diag_counts_df);   datatable(rv$diag_counts_df, options = list(pageLength = 10, scrollX = TRUE)) })
  output$diag_eligible  <- renderDT({ req(rv$diag_eligible_df); datatable(rv$diag_eligible_df, options = list(pageLength = 20, scrollX = TRUE)) })
  output$tbl_edger      <- renderDT({ req(rv$edger_tab);        datatable(rv$edger_tab, options = list(pageLength = 15, scrollX = TRUE)) })
  output$tbl_stats      <- renderDT({ req(rv$organ_stats);      datatable(rv$organ_stats, options = list(pageLength = 15, scrollX = TRUE)) })
  
  output$download_zip <- downloadHandler(
    filename = function() paste0("diff_homing_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      req(rv$files_written)
      owd <- setwd(rv$outdir); on.exit(setwd(owd), add = TRUE)
      zip::zip(zipfile = file, files = list.files(rv$outdir, recursive = TRUE, all.files = FALSE, full.names = FALSE))
    }
  )
}

shinyApp(ui = ui, server = server)
