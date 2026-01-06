# changer_app — SERVER

library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(DT)
library(ggplot2)
library(grid)              # arrows
library(shinycssloaders)   # (7) spinners
library("shinyalert")

# NEW: SQLite
library(DBI)
library(RSQLite)

# --- Rutas (app corre dentro de Dashboard/) -----------------------------------
DATA_DIR    <- "."
GENE_DIR    <- file.path(DATA_DIR, "by_gene_MANE_SELECT")
SUMMARY_TSV <- file.path(DATA_DIR, "gene_summary.tsv")
GENE_NAMES  <- file.path(DATA_DIR, "gene_names.txt")
# NEW: base consolidada
GENE_DB     <- file.path(DATA_DIR, "by_gene.sqlite")

# --- Utilidades ---------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_col <- function(df, nm, default = NA) {
  if (!nm %in% names(df)) default else df[[nm]]
}

# Flags y AF (sin imputar)
compute_flags <- function(df) {
  REVEL      <- suppressWarnings(as.numeric(safe_col(df, "REVEL")))
  CADD_PHRED <- suppressWarnings(as.numeric(safe_col(df, "CADD_PHRED")))
  AM_CLASS   <- tolower(as.character(safe_col(df, "am_class")))
  EVE_CLASS  <- tolower(as.character(safe_col(df, "EVE_CLASS")))
  LoF        <- as.character(safe_col(df, "LoF"))
  AC         <- suppressWarnings(as.numeric(safe_col(df, "AC")))
  AN         <- suppressWarnings(as.numeric(safe_col(df, "AN")))
  gAC        <- suppressWarnings(as.numeric(safe_col(df, "gnomADj_AC_joint")))
  
  CHANGER_AF <- ifelse(!is.na(AC) & !is.na(AN) & AN > 0, AC / AN, NA_real_)
  not_in_gnomad <- is.na(gAC) | gAC == 0
  
  miss_damaging <- (!is.na(REVEL)      & REVEL >= 0.75) |
    (!is.na(CADD_PHRED) & CADD_PHRED >= 20) |
    (!is.na(AM_CLASS)   & AM_CLASS %in% c("pathogenic","likely_pathogenic")) |
    (!is.na(EVE_CLASS)  & EVE_CLASS == "pathogenic") |
    (!is.na(LoF)        & LoF == "HC")
  
  common_changer <- (!is.na(AC) & AC >= 9)
  
  tibble(
    CHANGER_AF       = CHANGER_AF,
    FLAG_NOT_GNOMAD  = not_in_gnomad,
    FLAG_DAMAGING    = miss_damaging,
    FLAG_COMMON      = common_changer
  )
}

clean_consequence <- function(x) {
  if (is.null(x)) return(x)
  y <- as.character(x)
  y <- gsub("&intron_variant", "", y, fixed = TRUE)
  y <- gsub("_tract_variant", "", y, fixed = TRUE)
  y <- gsub("_variant", "", y, fixed = TRUE)
  y <- gsub("__+", "_", y)
  y <- gsub("^_", "", y)
  y <- gsub("_$", "", y)
  y
}

read_gene_summary <- function() {
  if (!file.exists(SUMMARY_TSV)) return(NULL)
  suppressMessages(read_tsv(SUMMARY_TSV, show_col_types = FALSE))
}

read_gene_names <- function() {
  if (!file.exists(GENE_NAMES)) return(NULL)
  df <- tryCatch(suppressMessages(readr::read_delim(GENE_NAMES, delim = "\t", show_col_types = FALSE)),
                 error = function(...) NULL)
  if (is.null(df)) df <- tryCatch(suppressMessages(readr::read_delim(GENE_NAMES, delim = ",", show_col_types = FALSE)),
                                  error = function(...) NULL)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  nms <- tolower(names(df))
  col_sym <- which(nms %in% c("symbol","gene","genesymbol","symbol_id"))
  col_nm  <- which(nms %in% c("name","gene_name","fullname","description"))
  if (length(col_sym) == 0) col_sym <- 1
  if (length(col_nm)  == 0) col_nm  <- 2
  out <- df[, c(col_sym[1], col_nm[1])]
  names(out) <- c("SYMBOL","FULL_NAME")
  out
}

get_full_gene_name <- function(sym, tbl) {
  if (is.null(tbl) || is.null(sym)) return(sym)
  hit <- tbl$FULL_NAME[match(sym, tbl$SYMBOL)]
  ifelse(is.na(hit) | !nzchar(hit), sym, hit)
}

# NEW: preferir SQLite para la lista de genes
list_available_genes <- function() {
  if (file.exists(GENE_DB)) {
    con <- dbConnect(RSQLite::SQLite(), GENE_DB)
    on.exit(dbDisconnect(con), add = TRUE)
    genes <- dbGetQuery(con, "SELECT DISTINCT SYMBOL FROM variants ORDER BY SYMBOL")[[1]]
    return(genes %||% character(0))
  }
  gs <- read_gene_summary()
  if (!is.null(gs)) {
    col <- if ("gene" %in% names(gs)) "gene" else if ("SYMBOL" %in% names(gs)) "SYMBOL" else NULL
    if (!is.null(col)) return(sort(unique(gs[[col]])))
  }
  if (dir.exists(GENE_DIR)) {
    fs <- list.files(GENE_DIR, pattern = "\\.tsv$", full.names = FALSE)
    return(unique(sub("\\.tsv$", "", fs)))
  }
  character(0)
}

# NEW: preferir SQLite para leer variantes
read_gene_tsv <- function(symbol) {
  if (file.exists(GENE_DB)) {
    con <- dbConnect(RSQLite::SQLite(), GENE_DB)
    on.exit(dbDisconnect(con), add = TRUE)
    return(dbGetQuery(con, "SELECT * FROM variants WHERE SYMBOL = ?", params = list(symbol)))
  }
  f <- file.path(GENE_DIR, paste0(symbol, ".tsv"))
  if (!file.exists(f)) return(NULL)
  suppressMessages(read_tsv(f, show_col_types = FALSE))
}

# --- HGVSc helpers (PB) -------------------------------------------------------
parse_hgvsc <- function(x) {
  x <- as.character(x)
  m_intr <- regexec("c\\.([0-9]+)([+-])([0-9]+)", x);  r_intr <- regmatches(x, m_intr)
  m_utr5 <- regexec("c\\.-([0-9]+)", x);               r_utr5 <- regmatches(x, m_utr5)
  m_utr3 <- regexec("c\\.\\*([0-9]+)", x);             r_utr3 <- regmatches(x, m_utr3)
  m_exon <- regexec("c\\.([0-9]+)(?![+-]|\\*)", x, perl = TRUE); r_exon <- regmatches(x, m_exon)
  n <- length(x)
  cpos <- rep(NA_real_, n); sign <- rep(NA_character_, n); offset <- rep(NA_real_, n)
  is_int <- rep(FALSE, n); utr <- rep(NA_character_, n); is_ex <- rep(FALSE, n)
  for (i in seq_len(n)) {
    if (length(r_intr[[i]]) == 4) { cpos[i] <- as.numeric(r_intr[[i]][2]); sign[i] <- r_intr[[i]][3]; offset[i] <- as.numeric(r_intr[[i]][4]); is_int[i] <- TRUE; next }
    if (length(r_utr5[[i]]) == 2) { utr[i] <- "5p"; next }
    if (length(r_utr3[[i]]) == 2) { utr[i] <- "3p"; next }
    if (length(r_exon[[i]]) >= 2) { cpos[i] <- as.numeric(r_exon[[i]][2]); is_ex[i] <- TRUE; next }
  }
  data.frame(cpos=cpos, sign=sign, offset=offset, is_intronic=is_int, utr=utr, is_exonic=is_ex, stringsAsFactors=FALSE)
}

build_transcript_layout <- function(gs_row) {
  stopifnot(nrow(gs_row) == 1)
  pick <- function(nms) { nm <- intersect(nms, names(gs_row)); if (length(nm) == 0) return(NULL); gs_row[[nm[1]]] }
  n_exons   <- suppressWarnings(as.numeric(pick(c("num_coding_exons","NUM_CODING_EXONS"))))
  cds_len   <- suppressWarnings(as.numeric(pick(c("cds_length","CDS_LENGTH"))))
  blk_sizes <- as.character(pick(c("blockSizes","BLOCKSIZES")))
  strand    <- as.character(pick(c("strand","STRAND")))
  validate(need(!is.na(cds_len)  && cds_len  > 0, "cds_length missing/invalid"))
  validate(need(!is.na(blk_sizes) && nzchar(blk_sizes), "blockSizes missing/invalid"))
  validate(need(!is.na(strand) && strand %in% c("+","-"), "strand missing/invalid (+/- expected)"))
  exon_sizes_raw <- strsplit(blk_sizes, ",", fixed = TRUE)[[1]]
  exon_sizes <- suppressWarnings(as.numeric(trimws(exon_sizes_raw)))
  exon_sizes <- exon_sizes[!is.na(exon_sizes) & exon_sizes > 0]
  validate(need(length(exon_sizes) > 0, "blockSizes has no valid exon sizes"))
  n_ex <- length(exon_sizes)
  intron_size <- ((cds_len + 3) * 0.2) / n_ex
  half_i      <- intron_size * 1.0
  starts <- numeric(n_ex); ends <- numeric(n_ex); x <- half_i
  for (j in seq_len(n_ex)) { starts[j] <- x; ends[j] <- x + exon_sizes[j]; x <- ends[j] + intron_size }
  total_len <- ends[n_ex] + half_i
  exons_df <- data.frame(exon_index_genomic = seq_len(n_ex), xmin = starts, xmax = ends, xmid = (starts + ends)/2, stringsAsFactors = FALSE)
  tr2gen <- if (strand == "+") seq_len(n_ex) else rev(seq_len(n_ex))
  labels_exon <- if (strand == "+") seq_len(n_ex) else rev(seq_len(n_ex))
  exon_sizes_tr <- if (strand == "+") exon_sizes else rev(exon_sizes)
  list(exons=exons_df, exon_sizes_tr=exon_sizes_tr, tr2gen=tr2gen, intron_size=intron_size, half_i=half_i, total_len=total_len, strand=strand, labels_exon=labels_exon)
}

map_variants_to_x <- function(hgvsc, layout, consequence = NULL) {
  cons_lc <- if (is.null(consequence)) rep(NA_character_, length(hgvsc)) else tolower(as.character(consequence))
  p <- parse_hgvsc(hgvsc)
  ex_sizes_tr <- layout$exon_sizes_tr
  n_ex        <- length(ex_sizes_tr)
  cum_end_tr  <- cumsum(ex_sizes_tr)
  cum_prev_tr <- c(0, cum_end_tr[-n_ex])
  exon_start_x <- function(j_gen) layout$exons$xmin[j_gen]
  exon_end_x   <- function(j_gen) layout$exons$xmax[j_gen]
  exon_of_cpos_tr <- function(cp) { if (is.na(cp)) return(NA_integer_); which(cum_end_tr >= cp)[1] }
  x_var <- rep(NA_real_, nrow(p))
  intron_id <- rep(NA_character_, nrow(p))
  for (i in seq_len(nrow(p))) {
    if (!is.na(cons_lc[i]) && grepl("upstream_gene_variant", cons_lc[i], fixed = TRUE)) {
      left_mid  <- layout$half_i * 0.5
      right_mid <- layout$total_len - layout$half_i * 0.5
      x_var[i] <- if (layout$strand == "+") left_mid else right_mid
      intron_id[i] <- "ext_first"; next
    }
    if (!is.na(cons_lc[i]) && grepl("downstream_gene_variant", cons_lc[i], fixed = TRUE)) {
      left_mid  <- layout$half_i * 0.5
      right_mid <- layout$total_len - layout$half_i * 0.5
      x_var[i] <- if (layout$strand == "+") right_mid else left_mid
      intron_id[i] <- "ext_last"; next
    }
    if (!is.na(p$utr[i])) {
      if (layout$strand == "+") x_var[i] <- if (p$utr[i] == "5p") layout$half_i else layout$total_len - layout$half_i
      else                      x_var[i] <- if (p$utr[i] == "5p") layout$total_len - layout$half_i else layout$half_i
      intron_id[i] <- paste0("utr_", p$utr[i]); next
    }
    cp <- p$cpos[i]; j_tr <- exon_of_cpos_tr(cp); if (is.na(j_tr)) next
    j_gen <- layout$tr2gen[j_tr]; if (is.na(j_gen)) next
    if (isTRUE(p$is_intronic[i])) {
      if (p$sign[i] == "+") x_var[i] <- if (layout$strand == "+") exon_end_x(j_gen) + layout$intron_size/2 else exon_start_x(j_gen) - layout$intron_size/2
      else                  x_var[i] <- if (layout$strand == "+") exon_start_x(j_gen) - layout$intron_size/2 else exon_end_x(j_gen) + layout$intron_size/2
      intron_id[i] <- paste0("int_tr_", j_tr)
    } else if (isTRUE(p$is_exonic[i])) {
      delta <- cp - cum_prev_tr[j_tr]
      x_var[i] <- if (layout$strand == "+") exon_start_x(j_gen) + delta else exon_end_x(j_gen) - delta
    }
  }
  if (any(!is.na(intron_id))) {
    df <- data.frame(idx = seq_along(x_var), intron_id = intron_id, x = x_var)
    df <- df %>% dplyr::group_by(intron_id) %>% dplyr::mutate(n = dplyr::n(), rn = dplyr::row_number(),
                                                              offset = ifelse(n > 1, (rn - (n + 1)/2) * (layout$intron_size * 0.2 / pmax(n - 1, 1)), 0),
                                                              x = x + ifelse(is.na(intron_id), 0, offset)) %>% dplyr::ungroup()
    x_var[df$idx] <- df$x
  }
  x_var
}

# --- Server -------------------------------------------------------------------
shinyServer(function(input, output, session) {
  shinyalert(
    title = "Welcome/Bienvenido! ",
    text = "CHANGER permite explorar variantes genéticas en la población chilena y acceder a nuestros datos de referencia para investigación e interpretación de variantes. Toda la data proviene de individuos deidentificados y debidamente consentidos. 
    Si estás interesado en contribuir con datos de exoma o genoma para ayudar a expandir este valioso recurso para comunidades subrepresentadas. Escríbenos a changer@udd.cl!",
    size = "m", closeOnEsc = TRUE, closeOnClickOutside = FALSE,
    html = FALSE, type = "success", showConfirmButton = TRUE, showCancelButton = FALSE,
    confirmButtonText = "Comenzar", confirmButtonCol = "#AEDEF4", timer = 0, imageUrl = "", animation = TRUE
  )
  
  GS     <- read_gene_summary()
  GN_TBL <- read_gene_names()
  
  # Autocomplete
  genes <- reactiveVal(character(0))
  observe({
    choices <- list_available_genes()
    genes(choices)
    updateSelectizeInput(session, "gene_search", choices = choices, selected = character(0), server = TRUE)
  })
  
  # Navegación
  observeEvent(input$go_search, { req(input$gene_search, nzchar(input$gene_search)); updateTabsetPanel(session, "tabs", selected = "results") })
  observeEvent(input$new_search, { updateTabsetPanel(session, "tabs", selected = "home") })
  
  # Estado
  current_gene <- reactive({ req(input$gene_search, nzchar(input$gene_search)); input$gene_search })
  
  variants_raw <- reactive({
    g  <- current_gene()
    df <- read_gene_tsv(g)
    validate(need(!is.null(df), paste0("No TSV/DB rows found for gene: ", g)))
    df
  })
  
  variants <- reactive({
    df <- variants_raw()
    flags <- compute_flags(df)
    df %>% mutate(
      CHANGER_AF       = flags$CHANGER_AF,
      FLAG_NOT_GNOMAD  = flags$FLAG_NOT_GNOMAD,
      FLAG_DAMAGING    = flags$FLAG_DAMAGING,
      FLAG_COMMON      = flags$FLAG_COMMON
    )
  })
  
  # PA - Consequence
  output$consequence_picker <- renderUI({
    df <- variants_raw()
    cons <- sort(na.omit(unique(as.character(safe_col(df, "Consequence")))))
    cons_labels <- clean_consequence(cons)
    checkboxGroupInput("consequence", "Consequence (annotation):",
                       choices = setNames(cons, cons_labels), selected = NULL, inline = FALSE)
  })
  
  # Filtros
  filtered_variants <- reactive({
    df <- variants()
    if (!is.null(input$consequence) && length(input$consequence) > 0) {
      df <- df %>% filter(as.character(Consequence) %in% input$consequence)
    }
    fset <- input$filters %||% character(0)
    if ("f_lofhc"    %in% fset) df <- df %>% filter(LoF == "HC")
    if ("f_dmg"      %in% fset) df <- df %>% filter(FLAG_DAMAGING)
    if ("f_nognomad" %in% fset) df <- df %>% filter(FLAG_NOT_GNOMAD)
    if ("f_common"   %in% fset) df <- df %>% filter(FLAG_COMMON)
    df
  })
  
  # PC — resumen
  output$gene_summary_card <- renderDT({
    gs <- GS; g <- current_gene()
    if (is.null(gs)) return(datatable(data.frame(Message = "gene_summary.tsv not found"), options = list(dom="t"), rownames = FALSE))
    col_gene <- if ("gene" %in% names(gs)) "gene" else if ("SYMBOL" %in% names(gs)) "SYMBOL" else NA_character_
    if (is.na(col_gene)) return(datatable(data.frame(Message = "Column 'gene' not found in gene_summary.tsv"), options = list(dom="t"), rownames = FALSE))
    row <- gs %>% filter(.data[[col_gene]] == g)
    if (nrow(row) == 0) return(datatable(data.frame(Message = paste("No summary found for", g)), options = list(dom="t"), rownames = FALSE))
    keep <- intersect(
      c("gene","Transcript","cds_length","num_coding_exons",
        "n_variants","n_synonymous","n_missense","n_missense_damaging",
        "n_loftee_hc","n_not_in_gnomAD","n_changer_common"),
      names(row)
    )
    if (length(keep) == 0) {
      df0 <- as.data.frame(t(row[1, , drop = FALSE]))
      return(datatable(df0, options = list(dom="t", paging=FALSE, ordering=FALSE, info=FALSE),
                       class = "display stripe hover compact", escape = FALSE))
    }
    labels <- c(
      gene = "<b>Gene</b>",
      Transcript = "<b>MANE transcript:</b>",
      cds_length = "<b>CDS length (bp):</b>",
      num_coding_exons = "<b>Coding exons:</b>",
      n_variants = "<b>Total variants:</b>",
      n_synonymous = "<b>Synonymous:</b>",
      n_missense = "<b>Missense:</b>",
      n_missense_damaging = "<b>Missense damaging:</b>",
      n_loftee_hc = "<b>Loss of function (HC):</b>",
      n_not_in_gnomAD = "<b>Not in gnomAD:</b>",
      n_changer_common = "<b>Common in Changer (≥1%):</b>"
    )
    df <- data.frame(Metric = unname(labels[keep]), Value  = as.vector(t(row[1, keep, drop = FALSE])), stringsAsFactors = FALSE)
    sym_col <- if ("gene" %in% names(row)) "gene" else if ("SYMBOL" %in% names(row)) "SYMBOL" else NA_character_
    if (!is.na(sym_col)) {
      idx_gene <- which(keep == sym_col)[1]
      if (length(idx_gene) == 1 && !is.na(idx_gene)) {
        sym <- as.character(row[[sym_col]][1])
        if (!is.na(sym) && nzchar(sym)) {
          url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sym)
          df$Value[idx_gene] <- sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', url, sym)
        }
      }
    }
    datatable(df, escape=FALSE, rownames=FALSE, class="display stripe hover compact",
              options=list(dom="t", paging=FALSE, ordering=FALSE, info=FALSE, scrollX=TRUE,
                           columnDefs=list(list(width="55%", targets=0))),
              callback = DT::JS("$(table.table().header()).hide();"))
  })
  
  # PD — tabla
  output$variants_table <- renderDT({
    df <- filtered_variants()
    if ("HGVSc" %in% names(df)) {
      v <- as.character(df$HGVSc); parts <- strsplit(v, ":", fixed = TRUE)
      last  <- sapply(parts, function(x) if (length(x) == 0) NA_character_ else tail(x, 1))
      df$HGVSc <- ifelse(is.na(last) | !nzchar(last), df$HGVSc, last)
    }
    if ("HGVSp" %in% names(df)) {
      v <- as.character(df$HGVSp); parts <- strsplit(v, ":", fixed = TRUE)
      last  <- sapply(parts, function(x) if (length(x) == 0) NA_character_ else tail(x, 1))
      last  <- gsub("%3D", "=", last, fixed = TRUE)
      df$HGVSp <- ifelse(is.na(last) | !nzchar(last), df$HGVSp, paste0(last))
    }
    if ("CHANGER_AF" %in% names(df)) {
      df$CHANGER_AF <- ifelse(is.na(df$CHANGER_AF), NA_character_,
                              formatC(signif(as.numeric(df$CHANGER_AF)), format="e", digits=2))
    }
    if ("gnomADj_AF_joint" %in% names(df)) {
      df$gnomADj_AF_joint <- ifelse(is.na(df$gnomADj_AF_joint), NA_character_,
                                    formatC(signif(as.numeric(df$gnomADj_AF_joint)), format="e", digits=2))
    }
    if ("Consequence" %in% names(df)) df$Consequence <- clean_consequence(df$Consequence)
    chrom_clean <- if ("#CHROM" %in% names(df)) gsub("^chr", "", as.character(df[["#CHROM"]]), ignore.case = TRUE) else NA_character_
    pos <- safe_col(df, "POS"); ref <- safe_col(df, "REF"); alt <- safe_col(df, "ALT")
    variant_id <- ifelse(!is.na(chrom_clean) & !is.na(pos) & !is.na(ref) & !is.na(alt),
                         paste(chrom_clean, pos, ref, alt, sep = "-"), NA_character_)
    flag_not_g <- safe_col(df, "FLAG_NOT_GNOMAD")
    in_gnomad  <- !(flag_not_g %in% c(TRUE, "TRUE", 1)); in_gnomad[is.na(in_gnomad)] <- FALSE
    df$Variant <- ifelse(!is.na(variant_id) & in_gnomad,
                         paste0('<a href="https://gnomad.broadinstitute.org/variant/', variant_id,
                                '?dataset=gnomad_r4" target="_blank" rel="noopener">', variant_id, '</a>'),
                         variant_id)
    preferred_cols <- c("Variant","AC","AN","AC_Hom","CHANGER_AF","gnomADj_AF_joint",
                        "Consequence","HGVSc","HGVSp","EXON","Protein_position","LoF")
    show_cols <- intersect(preferred_cols, names(df)); if (length(show_cols) == 0) show_cols <- names(df)
    show_cols <- setdiff(show_cols, grep("^FLAG_", names(df), value = TRUE))
    show_cols <- setdiff(show_cols, c("SYMBOL","Gene","gnomADj_AC_joint","gnomADj_AN_joint","gnomADj_nhomalt_joint"))
    df_sub <- df[, show_cols, drop = FALSE]
    to_badge <- function(txt, bg, title) {
      paste0('<span title="', title,
             '" style="display:inline-block;min-width:1.2em;padding:0 6px;margin-right:4px;border-radius:12px;',
             'background:', bg, ';color:#fff;font-weight:600;line-height:1.2;">', txt, '</span>')
    }
    flag_ng <- ifelse(safe_col(df, "FLAG_NOT_GNOMAD") %in% TRUE, to_badge("Novel", "#6c757d", "Not in gnomAD"), "")
    flag_d  <- ifelse(safe_col(df, "FLAG_DAMAGING")  %in% TRUE, to_badge("Damaging", "#d62728", "Damaging"), "")
    flag_c  <- ifelse(safe_col(df, "FLAG_COMMON")    %in% TRUE, to_badge("Common", "#2ca02c", "Common (AF ≥ 1%)"), "")
    df_sub$Flags <- paste0(flag_ng, flag_d, flag_c); df_sub$Flags[nchar(df_sub$Flags) == 0] <- ""
    df_sub <- dplyr::relocate(df_sub, Flags, .after = "Variant")
    display_names <- c(
      "Variant"="Variant","Flags"="Flags","AC"="Allele count","AN"="Allele number","AC_Hom"="Number homozygotes",
      "CHANGER_AF"="Allele frequency (AF)","gnomADj_AF_joint"="gnomAD AF","Consequence"="VEP Consequence",
      "HGVSc"="HGVSc","HGVSp"="HGVSp","EXON"="Exon","Protein_position"="Protein position","LoF"="Loss of Function"
    )
    labels <- names(df_sub)
    labels[labels %in% names(display_names)] <- unname(display_names[labels[labels %in% names(display_names)]])
    datatable(df_sub, rownames=FALSE, filter="none", colnames=labels, escape=FALSE,
              class="display stripe hover nowrap",
              options=list(paging=FALSE, scrollY="40vh", scrollX=TRUE, scrollCollapse=TRUE, dom="t",
                           deferRender=TRUE,
                           columnDefs=list(list(className="dt-center", width="120px", targets=1),
                                           list(orderable=FALSE, targets=1)))
    )
  })
  
  # PB — título y plot
  output$pb_title <- renderUI({
    g <- current_gene()
    full_name <- get_full_gene_name(g, GN_TBL)
    h3(paste0(g, ": ", full_name), class = "card-title", style = "text-align:center;")
  })
  
  output$pb_lollipop <- renderPlot({
    g  <- current_gene()
    gs <- read_gene_summary()
    validate(need(!is.null(gs), "gene_summary.tsv not found"))
    col_gene <- if ("gene" %in% names(gs)) "gene" else if ("SYMBOL" %in% names(gs)) "SYMBOL" else NA_character_
    validate(need(!is.na(col_gene), "Column 'gene' not found in gene_summary.tsv"))
    row <- gs %>% dplyr::filter(.data[[col_gene]] == g)
    validate(need(nrow(row) == 1, paste("No summary found for", g)))
    L <- build_transcript_layout(row)
    dfv <- filtered_variants()
    validate(need(!is.null(dfv) && nrow(dfv) > 0, "No variants (after filters) for this gene"))
    validate(need("HGVSc" %in% names(dfv), "HGVSc column is required"))
    cons_raw <- tolower(as.character(safe_col(dfv, "Consequence")))
    lof_col  <- as.character(safe_col(dfv, "LoF"))
    var_class <- ifelse(!is.na(lof_col) & lof_col == "HC", "LoF",
                        ifelse(!is.na(cons_raw) & grepl("missense", cons_raw), "missense",
                               ifelse(!is.na(cons_raw) & grepl("synonymous", cons_raw), "synonymous", "other")))
    dfv$VarClass <- factor(var_class, levels = c("LoF","missense","synonymous","other"))
    dfv$vx <- map_variants_to_x(dfv$HGVSc, L, consequence = safe_col(dfv, "Consequence"))
    dfv <- dfv[!is.na(dfv$vx), , drop = FALSE]
    pparse <- parse_hgvsc(dfv$HGVSc)
    cons_lc <- tolower(as.character(safe_col(dfv, "Consequence")))
    intr_flag <- (pparse$is_intronic %in% TRUE) | (!is.na(pparse$utr)) |
      grepl("upstream_gene_variant", cons_lc, fixed = TRUE) |
      grepl("downstream_gene_variant", cons_lc, fixed = TRUE)
    intr_flag[is.na(intr_flag)] <- FALSE
    df_ex  <- dfv[!intr_flag, , drop = FALSE]
    df_int <- dfv[ intr_flag, , drop = FALSE]
    exon_h <- 0.25; y0 <- 0; pt_sz <- 3.6
    pal <- c("LoF"="#D62728","missense"="#FF7F0E","synonymous"="#1F77B4","other"="grey80")
    p <- ggplot() +
      geom_segment(aes(x = 0, xend = L$total_len, y = y0, yend = y0), linewidth = 0.4) +
      geom_rect(data = L$exons, aes(xmin = xmin, xmax = xmax, ymin = -exon_h, ymax = exon_h),
                fill = "grey70", color = NA) +
      { if (nrow(df_ex) > 0) geom_segment(data = df_ex, aes(x = vx, xend = vx, y = exon_h, yend = exon_h + 0.7, color = VarClass)) else NULL } +
      { if (nrow(df_ex) > 0) geom_point(  data = df_ex, aes(x = vx, y = exon_h + 0.7, color = VarClass), size = pt_sz) else NULL } +
      { if (nrow(df_int) > 0) geom_segment(data = df_int, aes(x = vx, xend = vx, y = y0, yend = y0 + 0.7, color = VarClass), linewidth = 0.9) else NULL } +
      { if (nrow(df_int) > 0) geom_point(  data = df_int, aes(x = vx, y = y0 + 0.7, color = VarClass), size = pt_sz) else NULL } +
      geom_text(data = L$exons, aes(x = xmid, y = -exon_h - 0.4, label = L$labels_exon), size = 6) +
      coord_cartesian(xlim = c(0, L$total_len), ylim = c(-1.2, 2)) +
      theme_minimal(base_size = 11) +
      theme(axis.title = element_blank(), axis.text = element_blank(),
            panel.grid = element_blank(), legend.position = "top",
            legend.text  = element_text(size = rel(1.75)),
            legend.title = element_text(size = rel(1.75)),
            plot.margin = margin(5, 10, 5, 10)) +
      scale_color_manual(values = pal, drop = FALSE, name = "Variants observed:")
    ah <- arrow(type = "closed", length = unit(12, "pt"))
    arrow_len <- 0.35 * L$intron_size
    if (L$strand == "+") p <- p + annotate("segment", x = L$total_len - arrow_len, xend = L$total_len, y = y0, yend = y0, arrow = ah, linewidth = 0.6)
    else                 p <- p + annotate("segment", x = arrow_len, xend = 0, y = y0, yend = y0, arrow = ah, linewidth = 0.6)
    p
  })
  
  # Download (DB o TSV)
  output$dl_gene_tsv <- downloadHandler(
    filename = function() paste0(current_gene(), ".tsv"),
    content = function(file) {
      df <- read_gene_tsv(current_gene())
      if (is.null(df)) stop("No TSV found for current gene.")
      readr::write_tsv(df, file)
    }
  )
})
