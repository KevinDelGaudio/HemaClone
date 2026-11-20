# Render one HemaClone report (Snakemake rule ch_report).
#
#   Rscript render_report.R RMD OUTPUT_HTML RESULTS_DIR TMP_DIR EXCLUDE_FILE
#
# The project folder is the working directory. EXCLUDE_FILE is "" when no
# samples are excluded. Pandoc is taken from the same environment as Rscript.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("usage: Rscript render_report.R RMD OUTPUT_HTML RESULTS_DIR TMP_DIR EXCLUDE_FILE")
}
rmd <- args[[1]]
output_html <- args[[2]]
results_dir <- args[[3]]
tmp_dir <- args[[4]]
exclude_file <- args[[5]]

if (!file.exists(rmd)) stop("Report source not found: ", rmd)
# R.home() is <environment>/lib/R; pandoc is in <environment>/bin
env_bin <- normalizePath(file.path(R.home(), "..", "..", "bin"))
if (!file.exists(file.path(env_bin, "pandoc"))) stop("pandoc not found in ", env_bin)
Sys.setenv(RSTUDIO_PANDOC = env_bin)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Knit from a copy in the temporary folder: the built-in source sits in the
# read-only image, and knitr writes next to the file it renders
work_rmd <- file.path(tmp_dir, basename(rmd))
invisible(file.copy(rmd, work_rmd, overwrite = TRUE))

# Absolute paths, resolved here: render() changes the working directory before
# it evaluates its arguments
project_dir <- getwd()
results_dir <- normalizePath(results_dir)
output_dir <- normalizePath(dirname(output_html))

# Workflow version: VERSION in the workflow home, two folders above this script
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
workflow_home <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."))
version_file <- file.path(workflow_home, "VERSION")
workflow_version <- if (file.exists(version_file)) trimws(readLines(version_file, warn = FALSE)[1]) else "unknown"

rmarkdown::render(
  input = work_rmd,
  output_file = basename(output_html),
  output_dir = output_dir,
  intermediates_dir = tmp_dir,
  knit_root_dir = project_dir,
  params = list(project_dir = project_dir,
                results_dir = results_dir,
                exclude_samples = exclude_file,
                workflow_version = workflow_version),
  envir = new.env(),
  quiet = TRUE
)
message("Report written: ", output_html)
