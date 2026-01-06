# build_gene_db.R  (run from Dashboard/)
library(DBI)
library(RSQLite)
library(readr)

gene_dir <- "by_gene_MANE_SELECT"
db_path  <- "by_gene.sqlite"

con <- dbConnect(SQLite(), db_path)
on.exit(dbDisconnect(con), add = TRUE)

files <- list.files(gene_dir, pattern = "\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)

first <- TRUE
for (f in files) {
  sym <- sub("\\.tsv$", "", basename(f))
  df  <- readr::read_tsv(f, show_col_types = FALSE)
  if (!"SYMBOL" %in% names(df)) df$SYMBOL <- sym
  if (first) {
    dbWriteTable(con, "variants", df, overwrite = TRUE)
    dbExecute(con, "CREATE INDEX idx_variants_symbol ON variants(SYMBOL)")
    first <- FALSE
  } else {
    dbWriteTable(con, "variants", df, append = TRUE)
  }
}
message("SQLite built at: ", normalizePath(db_path))

