# CHANGER Dashboard (R Shiny)
Explore aggregated coding genetic variants from 902 Chilean exomes

---

## Contents
- [`ui.R`](ui.R) — UI layout.
- [`server.R`](server.R) — Server logic (loading, flags, filters, plot, variants table).
- [`build_gene_db.R`](build_gene_db.R) — Builds a single SQLite DB (`by_gene.sqlite`) from per-gene TSVs.
- [`gene_summary.tsv`](gene_summary.tsv) — Gene-level summary (MANE transcript metadata and counts).
- [`gene_names.txt`](gene_names.txt) — Optional map from `SYMBOL` → full gene name.
- [`by_gene_MANE_SELECT/`](by_gene_MANE_SELECT/) — with three examples per-gene TSVs for local testing.

---

## Requirements
R (≥ 4.1) and packages:
`shiny`, `shinythemes`, `shinyWidgets`, `DT`, `shinycssloaders`, `shinyalert`,
`readr`, `dplyr`, `stringr`, `ggplot2`, `grid`, `DBI`, `RSQLite`.

Install:
```r
pkgs <- c("shiny","shinythemes","shinyWidgets","DT","shinycssloaders","shinyalert",
          "readr","dplyr","stringr","ggplot2","grid","DBI","RSQLite")
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

## Quick start (with example TSVs)

Ensure gene_summary.tsv and the folder by_gene_MANE_SELECT/ are in the project root.
(Recommended) Build a single SQLite database:
Rscript build_gene_db.R
The app prefers SQLite when present; otherwise it falls back to TSVs + gene_summary.tsv.

Run the app:
shiny::runApp(".")
In the Search box, type a gene symbol matching one of the example TSV filenames (without .tsv).

**Data expectations (minimal)**
Gene summary: SYMBOL/gene, Transcript, strand, num_coding_exons, cds_length, blockSizes.
Per-gene TSVs: standard annotations (VEP, LOFTEE, REVEL, CADD, AlphaMissense, EVE, gnomAD) and AC, AN, #CHROM, POS, REF, ALT.

**Deploy (optional)**
Deploy to shinyapps.io:
rsconnect::deployApp(".")

If using SQLite, you can omit by_gene_MANE_SELECT/ from the bundle and deploy only the app files, by_gene.sqlite, and gene_summary.tsv.
