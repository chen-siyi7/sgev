# Run from any directory: Rscript /path/to/sgev/simulations/run.R full
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "smoke"
studies <- if (length(args) > 1L) args[-1L] else
  c("monitoring", "drift", "paired", "stopping", "near_margin", "adaptive_check")
if (!(mode %in% c("smoke", "full"))) stop("First argument must be smoke or full")
if (any(!studies %in% c("monitoring", "drift", "paired", "stopping", "near_margin", "adaptive_check"))) {
  stop("Unknown study name")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) stop("Run this program with Rscript")
# Some Rscript builds encode spaces in --file arguments as ~+~.
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
options(stringsAsFactors = FALSE)
options(sgev.simulation.mode = mode)
# Separate modes so a smoke run cannot replace the full-run output.
destination <- file.path(root, "output", mode)
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
setwd(destination)
for (directory in c("results", "figures")) {
  dir.create(directory, showWarnings = FALSE)
}
for (study in studies) {
  message("Running ", study, " (", mode, ")")
  env <- new.env(parent = globalenv())
  for (source in c("R/core.R", "R/api.R", "R/scalar.R",
                   "simulations/helpers.R", "simulations/geometry.R")) {
    sys.source(file.path(root, source), envir = env)
  }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  sys.source(file.path(root, "simulations", paste0(study, ".R")), envir = env)
}
needed <- c("sequential_comparisons.csv", "batching_comparisons.csv",
            "drift_comparisons.csv", "stopping_efficiency.csv", "near_margin.csv")
if (all(file.exists(file.path("results", needed)))) {
  sys.source(file.path(root, "simulations/plots.R"), envir = new.env(parent = globalenv()))
}
writeLines(capture.output(sessionInfo()), "session_info.txt")
message("Results written to ", destination)
