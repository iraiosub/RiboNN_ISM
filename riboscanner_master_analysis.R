#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (!capabilities("X11")) {
  options(bitmapType = "cairo")
}
options(error = function() {
  msg <- geterrmessage()
  cat("=== ERROR TRACEBACK ===\n", file = stderr())
  traceback(3, max.lines = 10)
  cat("======================\n", file = stderr())
})

analysis_script_version <- "2026-06-10-post-predict-join-fix"
message("RiboScanner analysis script version: ", analysis_script_version)

defaults <- list(
  mode = "full",
  riboscanner_dir = "/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/ag/riboscanner",
  fasta_gz = "master_table.context_m40_p40.fa.gz",
  test_fasta = "scn2a_mouse_centered.fa",
  master_table = "master_table.tsv.gz",
  prediction_output = "",
  analysis_dir = "riboscanner_analysis",
  riboscanner_exe = "RiboScanner",
  riboscanner_env = "",
  cxx_lib_dir = "",
  run_prediction = TRUE,
  min_sequence_length = 81,
  id_col = "",
  prediction_id_col = "",
  class_col = "",
  detected_col = "",
  uorf_col = "",
  gfp_col = "",
  score_col = "",
  uncertainty_col = ""
)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript riboscanner_master_analysis.R --mode test\n",
    "  Rscript riboscanner_master_analysis.R --mode full\n",
    "  Rscript riboscanner_master_analysis.R --mode full --no-predict --prediction-output PATH\n",
    "\n",
    "Main options:\n",
    "  --mode test|full\n",
    "  --riboscanner-dir PATH\n",
    "  --riboscanner-exe PATH_OR_COMMAND\n",
    "  --riboscanner-env PATH_TO_CONDA_ENV\n",
    "  --cxx-lib-dir PATH_CONTAINING_LIBSTDCXX\n",
    "  --analysis-dir NAME\n",
    "  --fasta PATH_OR_FILENAME\n",
    "  --test-fasta PATH_OR_FILENAME\n",
    "  --master-table PATH_OR_FILENAME\n",
    "  --prediction-output PATH\n",
    "  --no-predict\n",
    "  --min-sequence-length N\n",
    "\n",
    "Column overrides:\n",
    "  --id-col COL --prediction-id-col COL --class-col COL --detected-col COL\n",
    "  --uorf-col COL --gfp-col COL --score-col COL --uncertainty-col COL\n",
    sep = ""
  )
}

parse_args <- function(args) {
  cfg <- defaults
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      usage()
      quit(status = 0)
    }
    if (arg == "--no-predict") {
      cfg$run_prediction <- FALSE
      i <- i + 1
      next
    }
    if (!grepl("^--", arg)) {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for ", arg, call. = FALSE)
    }
    value <- args[[i + 1]]
    key <- sub("^--", "", arg)
    key <- gsub("-", "_", key)
    if (key == "fasta") key <- "fasta_gz"
    if (!key %in% names(cfg)) {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    cfg[[key]] <- value
    i <- i + 2
  }

  cfg$mode <- tolower(cfg$mode)
  if (!cfg$mode %in% c("test", "full")) {
    stop("--mode must be 'test' or 'full'", call. = FALSE)
  }
  if (cfg$mode == "test") {
    cfg$analysis_dir <- if (identical(cfg$analysis_dir, defaults$analysis_dir)) {
      "riboscanner_analysis_scn2a_test"
    } else {
      cfg$analysis_dir
    }
  }
  cfg$run_prediction <- tolower(as.character(cfg$run_prediction)) %in% c("true", "t", "1", "yes")
  cfg$min_sequence_length <- as.integer(cfg$min_sequence_length)
  cfg
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing R package '", pkg, "'. Install it before running this script.", call. = FALSE)
  }
}

for (pkg in c("dplyr", "ggplot2", "readr", "stringr", "tibble", "tidyr")) {
  need_pkg(pkg)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

is_blank <- function(x) {
  is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))
}

resolve_input_path <- function(x, root = cfg$riboscanner_dir) {
  if (is_blank(x)) return(NULL)
  path <- as.character(x)
  if (!grepl("^/|^[A-Za-z]:", path)) {
    path <- file.path(root, path)
  }
  normalizePath(path, mustWork = FALSE)
}

path_for_command <- function(path, root = cfg$riboscanner_dir) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = FALSE)
  prefix <- paste0(root, .Platform$file.sep)
  if (startsWith(path, prefix)) {
    return(file.path(".", substring(path, nchar(prefix) + 1)))
  }
  path
}

resolve_executable <- function(exe) {
  if (is_blank(exe)) return("")
  if (grepl("/", exe)) return(normalizePath(exe, mustWork = FALSE))
  found <- Sys.which(exe)
  if (!nzchar(found)) return(exe)
  normalizePath(found, mustWork = FALSE)
}

prepend_env_path <- function(name, value) {
  if (is_blank(value)) return(invisible(FALSE))
  current <- Sys.getenv(name)
  new_value <- value
  if (nzchar(current)) {
    new_value <- paste(value, current, sep = .Platform$path.sep)
  }
  do.call(Sys.setenv, stats::setNames(list(new_value), name))
  invisible(TRUE)
}

configure_riboscanner_env <- function(exe, env_prefix = "") {
  exe_path <- resolve_executable(exe)
  prefix <- env_prefix
  if (is_blank(prefix) && nzchar(exe_path)) {
    bin_dir <- dirname(exe_path)
    candidate <- dirname(bin_dir)
    if (basename(bin_dir) == "bin" && dir.exists(file.path(candidate, "lib"))) {
      prefix <- candidate
    }
  }

  if (!is_blank(prefix)) {
    prefix <- normalizePath(prefix, mustWork = FALSE)
    prepend_env_path("PATH", file.path(prefix, "bin"))
    prepend_env_path("LD_LIBRARY_PATH", file.path(prefix, "lib"))
    Sys.setenv(CONDA_PREFIX = prefix)
    message("Configured RiboScanner env: ", prefix)
    message("RiboScanner env lib first: ", file.path(prefix, "lib"))
  } else {
    message("RiboScanner env prefix was not detected; using current PATH/LD_LIBRARY_PATH.")
  }
  invisible(prefix)
}

torch_lib_paths <- function(prefix) {
  if (is_blank(prefix)) return(character())
  pattern <- file.path(prefix, "lib", "python*", "site-packages", "torch", "lib")
  Sys.glob(pattern)
}

libstdcxx_path <- function(lib_dir) {
  exact <- file.path(lib_dir, "libstdc++.so.6")
  if (file.exists(exact)) return(exact)
  matches <- Sys.glob(file.path(lib_dir, "libstdc++.so.6*"))
  matches <- matches[file.exists(matches)]
  if (length(matches) == 0) return(exact)
  sort(matches, decreasing = TRUE)[1]
}

libstdcxx_versions <- function(lib_dir) {
  lib <- libstdcxx_path(lib_dir)
  strings <- Sys.which("strings")
  if (!file.exists(lib) || !nzchar(strings)) return(character())
  suppressWarnings(system2(strings, lib, stdout = TRUE, stderr = FALSE))
}

libstdcxx_has_cxxabi <- function(lib_dir) {
  any(grepl("^CXXABI_1\\.3\\.11$", libstdcxx_versions(lib_dir)))
}

cxx_lib_candidates <- function(prefix, override = "") {
  candidates <- character()
  if (!is_blank(override)) {
    candidates <- c(candidates, override)
  }
  if (!is_blank(prefix)) {
    prefix <- normalizePath(prefix, mustWork = FALSE)
    candidates <- c(candidates, file.path(prefix, "lib"))
    base_prefix <- normalizePath(file.path(prefix, "..", ".."), mustWork = FALSE)
    candidates <- c(candidates, file.path(base_prefix, "lib"))
  }
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  if (nzchar(conda_prefix)) {
    candidates <- c(candidates, file.path(conda_prefix, "lib"))
  }
  unique(normalizePath(candidates, mustWork = FALSE))
}

choose_cxx_lib_dir <- function(prefix, override = "") {
  candidates <- cxx_lib_candidates(prefix, override)
  candidates <- candidates[vapply(candidates, function(candidate) {
    file.exists(libstdcxx_path(candidate))
  }, logical(1))]

  if (!is_blank(override) && !file.exists(libstdcxx_path(override))) {
    stop("--cxx-lib-dir does not contain libstdc++.so.6*: ", override, call. = FALSE)
  }

  for (candidate in candidates) {
    if (libstdcxx_has_cxxabi(candidate)) {
      message("Selected C++ runtime lib dir: ", candidate)
      return(candidate)
    }
  }

  if (length(candidates) > 0) {
    message("Checked libstdc++ candidates but none had CXXABI_1.3.11:")
    for (candidate in candidates) {
      versions <- libstdcxx_versions(candidate)
      cxxabi <- tail(versions[grepl("^CXXABI_", versions)], 8)
      message("  ", libstdcxx_path(candidate), " -> ", paste(cxxabi, collapse = ", "))
    }
    warning(
      "No candidate libstdc++.so.6 advertises CXXABI_1.3.11. ",
      "Force a newer runtime with --cxx-lib-dir, or run: conda install -p ",
      prefix, " -c conda-forge --override-channels 'libstdcxx-ng>=12' 'libgcc-ng>=12'",
      call. = FALSE
    )
    return(candidates[1])
  }

  if (!is_blank(prefix)) {
    warning(
      "No libstdc++.so.6 was found near the RiboScanner env. ",
      "Install one with: conda install -p ", prefix,
      " -c conda-forge --override-channels 'libstdcxx-ng>=12' 'libgcc-ng>=12'",
      call. = FALSE
    )
  }
  ""
}

riboscanner_child_env <- function(prefix, cxx_lib_dir = "") {
  if (is_blank(prefix)) return(character())
  prefix <- normalizePath(prefix, mustWork = FALSE)
  cxx_lib <- if (!is_blank(cxx_lib_dir)) libstdcxx_path(cxx_lib_dir) else ""
  ld_paths <- unique(c(
    cxx_lib_dir,
    file.path(prefix, "lib"),
    torch_lib_paths(prefix)
  ))
  ld_paths <- ld_paths[dir.exists(ld_paths)]

  c(
    paste0("PATH=", paste(c(file.path(prefix, "bin"), Sys.getenv("PATH")), collapse = .Platform$path.sep)),
    paste0("LD_LIBRARY_PATH=", paste(ld_paths, collapse = .Platform$path.sep)),
    if (file.exists(cxx_lib)) paste0("LD_PRELOAD=", cxx_lib) else NULL,
    paste0("CONDA_PREFIX=", prefix)
  )
}

check_libstdcxx <- function(lib_dir, prefix) {
  if (is_blank(lib_dir)) return(invisible(FALSE))
  lib <- libstdcxx_path(lib_dir)
  if (!file.exists(lib)) {
    warning(
      "No libstdc++.so.6* found at ", lib, ". ",
      "Install it with: conda install -p ", prefix,
      " -c conda-forge --override-channels 'libstdcxx-ng>=12' 'libgcc-ng>=12'",
      call. = FALSE
    )
    return(invisible(FALSE))
  }

  message("Selected libstdc++: ", normalizePath(lib, mustWork = FALSE))
  message("RiboScanner child LD_PRELOAD: ", normalizePath(lib, mustWork = FALSE))
  if (!libstdcxx_has_cxxabi(lib_dir)) {
    warning(
      "The selected libstdc++.so.6 does not advertise CXXABI_1.3.11. ",
      "Update it with: conda install -p ", prefix,
      " -c conda-forge --override-channels 'libstdcxx-ng>=12' 'libgcc-ng>=12'",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

preflight_torch_import <- function(prefix, env_vars, analysis_dir) {
  if (is_blank(prefix)) return(invisible(TRUE))
  python <- file.path(prefix, "bin", "python")
  if (!file.exists(python)) return(invisible(TRUE))

  stdout_path <- file.path(analysis_dir, "riboscanner_torch_preflight_stdout.log")
  stderr_path <- file.path(analysis_dir, "riboscanner_torch_preflight_stderr.log")
  script_path <- file.path(analysis_dir, "riboscanner_torch_preflight.py")
  writeLines(
    c(
      "import torch",
      "print('torch import OK', torch.__version__)"
    ),
    script_path
  )

  status <- system2(
    python,
    args = shQuote(script_path),
    env = env_vars,
    stdout = stdout_path,
    stderr = stderr_path
  )
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "RiboScanner torch preflight failed before prediction. See: ", stderr_path, "\n",
      "If stderr mentions CXXABI/libstdc++, pass --cxx-lib-dir PATH_TO_A_NEWER_LIB_DIR ",
      "or update the conda runtime in: ", prefix,
      call. = FALSE
    )
  }
  message("Torch preflight OK. Log: ", stdout_path)
  invisible(TRUE)
}

discover_one <- function(root, patterns, label) {
  candidates <- unique(unlist(lapply(patterns, function(pattern) {
    list.files(root, pattern = pattern, full.names = TRUE, recursive = FALSE)
  })))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    stop("Could not find ", label, " in ", root, call. = FALSE)
  }
  if (length(candidates) > 1) {
    candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
    message("Found multiple ", label, " files. Using latest: ", candidates[1])
  }
  normalizePath(candidates[1], mustWork = TRUE)
}

save_plot <- function(plot, filename, width, height, dpi = 160) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) {
    grDevices::png(filename, width = width, height = height, units = "in", res = dpi, type = "cairo")
    print(plot)
    grDevices::dev.off()
    message("Saved ", filename)
    return(invisible(filename))
  }

  pdf_path <- sub("\\.png$", ".pdf", filename)
  grDevices::pdf(pdf_path, width = width, height = height)
  print(plot)
  grDevices::dev.off()
  message("Cairo PNG is unavailable; saved PDF instead: ", pdf_path)
  invisible(pdf_path)
}

riboscanner_dir <- normalizePath(cfg$riboscanner_dir, mustWork = FALSE)
analysis_dir <- file.path(riboscanner_dir, cfg$analysis_dir)
fig_dir <- file.path(analysis_dir, "figures")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("RiboScanner directory: ", riboscanner_dir)
message("Analysis directory:    ", analysis_dir)
message("Mode:                  ", cfg$mode)
riboscanner_env_prefix <- configure_riboscanner_env(cfg$riboscanner_exe, cfg$riboscanner_env)
riboscanner_cxx_lib_dir <- choose_cxx_lib_dir(riboscanner_env_prefix, cfg$cxx_lib_dir)
if (!is_blank(riboscanner_cxx_lib_dir)) {
  prepend_env_path("LD_LIBRARY_PATH", riboscanner_cxx_lib_dir)
}
riboscanner_env_vars <- riboscanner_child_env(riboscanner_env_prefix, riboscanner_cxx_lib_dir)
check_libstdcxx(riboscanner_cxx_lib_dir, riboscanner_env_prefix)
if (length(riboscanner_env_vars)) {
  message("RiboScanner child LD_LIBRARY_PATH: ", sub("^LD_LIBRARY_PATH=", "", riboscanner_env_vars[grepl("^LD_LIBRARY_PATH=", riboscanner_env_vars)]))
}

selected_fasta <- if (cfg$mode == "test") cfg$test_fasta else cfg$fasta_gz
fasta_path <- resolve_input_path(selected_fasta, riboscanner_dir)
if (is.null(fasta_path)) {
  fasta_path <- discover_one(
    riboscanner_dir,
    patterns = c("\\.(fa|fasta|fna)\\.gz$", "\\.(fa|fasta|fna)$"),
    label = "FASTA/FASTA.GZ"
  )
}
master_path <- resolve_input_path(cfg$master_table, riboscanner_dir)
using_default_master_table <- identical(cfg$master_table, defaults$master_table)
if (is.null(master_path) || (using_default_master_table && !file.exists(master_path))) {
  master_path <- discover_one(
    riboscanner_dir,
    patterns = c(
      "^master[._]table$",
      "^master[._]table\\.(tsv|txt|csv)(\\.gz)?$",
      "master.*table.*\\.(tsv|txt|csv)(\\.gz)?$"
    ),
    label = "master table"
  )
}
if (!file.exists(master_path)) {
  stop(
    "Master table not found: ", master_path,
    ". Pass the correct file with --master-table PATH.",
    call. = FALSE
  )
}
message("FASTA:        ", fasta_path)
message("Master table: ", master_path)

wrap_seq <- function(seq, width = 80) {
  starts <- seq(1, nchar(seq), by = width)
  substring(seq, starts, pmin(starts + width - 1, nchar(seq)))
}

write_record <- function(out, header, seq, min_len) {
  seq_len <- nchar(seq)
  keep <- seq_len >= min_len
  if (keep) {
    writeLines(paste0(">", header), out)
    writeLines(wrap_seq(seq), out)
  }
  tibble(
    fasta_header = header,
    fasta_id = str_replace(header, "\\s.*$", ""),
    seq_length = seq_len,
    kept_for_prediction = keep
  )
}

filter_fasta_by_length <- function(input, output, min_len) {
  con <- if (str_detect(input, "\\.gz$")) gzfile(input, "rt") else file(input, "rt")
  out <- file(output, "wt")
  on.exit(close(con), add = TRUE)
  on.exit(close(out), add = TRUE)
  rows <- list()
  current_header <- NULL
  current_seq <- character()

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (substr(line, 1, 1) == ">") {
      if (!is.null(current_header)) {
        rows[[length(rows) + 1]] <- write_record(out, current_header, paste0(current_seq, collapse = ""), min_len)
      }
      current_header <- str_replace(line, "^>", "")
      current_seq <- character()
    } else if (!is.null(current_header)) {
      current_seq <- c(current_seq, str_replace_all(line, "\\s+", ""))
    }
  }

  if (!is.null(current_header)) {
    rows[[length(rows) + 1]] <- write_record(out, current_header, paste0(current_seq, collapse = ""), min_len)
  }
  bind_rows(rows)
}

input_stem <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(fasta_path)))
filtered_fasta <- file.path(analysis_dir, paste0(input_stem, ".min", cfg$min_sequence_length, ".fa"))
fasta_lengths <- filter_fasta_by_length(fasta_path, filtered_fasta, cfg$min_sequence_length)

length_summary <- fasta_lengths %>%
  summarise(
    total_sequences = n(),
    removed_lt_min_length = sum(!kept_for_prediction),
    kept_for_prediction = sum(kept_for_prediction),
    min_length = min(seq_length, na.rm = TRUE),
    median_length = median(seq_length, na.rm = TRUE),
    max_length = max(seq_length, na.rm = TRUE)
  )
print(length_summary)

if (sum(fasta_lengths$kept_for_prediction) == 0) {
  stop("No FASTA records are at least ", cfg$min_sequence_length, " bases.", call. = FALSE)
}
message("Filtered FASTA written to: ", filtered_fasta)

p_length <- ggplot(fasta_lengths, aes(seq_length, fill = kept_for_prediction)) +
  geom_histogram(bins = 60, color = "white", size = 0.2) +
  geom_vline(xintercept = cfg$min_sequence_length, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "#b2182b", "TRUE" = "#2166ac")) +
  labs(
    title = "Input sequence length QC",
    x = "Sequence length (nt)",
    y = "Number of sequences",
    fill = paste0("Length >= ", cfg$min_sequence_length)
  ) +
  theme_bw()
save_plot(p_length, file.path(fig_dir, "sequence_length_qc.png"), width = 8, height = 5)

candidate_prediction_files <- function(roots) {
  roots <- unique(roots[file.exists(roots)])
  files <- unique(unlist(lapply(roots, function(root) {
    list.files(root, pattern = "\\.(tsv|csv|txt)(\\.gz)?$", full.names = TRUE, recursive = FALSE)
  })))
  files <- files[file.exists(files)]
  files <- files[
    str_detect(tolower(basename(files)), "pred|prediction|riboscanner|score") &
      !str_detect(tolower(basename(files)), "log$")
  ]
  files[order(file.info(files)$mtime, decreasing = TRUE)]
}

prediction_output_path <- if (is_blank(cfg$prediction_output)) {
  file.path(analysis_dir, "riboscanner_predictions.tsv")
} else {
  resolve_input_path(cfg$prediction_output, riboscanner_dir)
}
candidate_roots <- unique(c(analysis_dir, riboscanner_dir, dirname(filtered_fasta), dirname(prediction_output_path)))
before_candidates <- candidate_prediction_files(candidate_roots)
stdout_path <- file.path(analysis_dir, "riboscanner_predict_stdout.log")
stderr_path <- file.path(analysis_dir, "riboscanner_predict_stderr.log")

if (isTRUE(cfg$run_prediction)) {
  preflight_torch_import(riboscanner_env_prefix, riboscanner_env_vars, analysis_dir)
  filtered_fasta_cmd <- path_for_command(filtered_fasta, riboscanner_dir)
  prediction_output_cmd <- path_for_command(prediction_output_path, riboscanner_dir)
  cmd_args <- c("predict", "--input", filtered_fasta_cmd, "--output", prediction_output_cmd)
  message("Running: ", cfg$riboscanner_exe, " ", paste(cmd_args, collapse = " "))
  old_wd <- getwd()
  setwd(riboscanner_dir)
  status <- system2(cfg$riboscanner_exe, args = cmd_args, env = riboscanner_env_vars, stdout = stdout_path, stderr = stderr_path)
  setwd(old_wd)
  if (!is.null(status) && !identical(status, 0L)) {
    stop("RiboScanner predict failed with status ", status, ". See: ", stderr_path, call. = FALSE)
  }
} else {
  message("Skipping RiboScanner predict because --no-predict was set.")
}

after_candidates <- candidate_prediction_files(candidate_roots)
new_candidates <- setdiff(after_candidates, before_candidates)
prediction_candidates <- unique(c(
  prediction_output_path,
  new_candidates,
  after_candidates
))
prediction_candidates <- prediction_candidates[file.exists(prediction_candidates)]
if (length(prediction_candidates) > 0) {
  message("Prediction candidate files checked:")
  for (candidate in prediction_candidates) message("  ", candidate)
}

read_delim_errors <- list()
file_read_description <- function(path) {
  if (!file.exists(path)) return("exists=FALSE")
  info <- file.info(path)
  paste0(
    "exists=TRUE, readable=", file.access(path, 4) == 0,
    ", size_bytes=", info$size
  )
}

readr_call <- function(fun, path) {
  fn <- getExportedValue("readr", fun)
  args <- list(file = path)
  supported_args <- names(formals(fn))
  if ("show_col_types" %in% supported_args) args$show_col_types <- FALSE
  if ("progress" %in% supported_args) args$progress <- FALSE
  suppressWarnings(do.call(fn, args))
}

read_delim_flexible <- function(path) {
  errors <- character()
  record_error <- function(name, message) {
    errors <<- c(errors, paste0(name, ": ", message))
  }
  record_result <- function(name, out) {
    if (!is.null(out) && ncol(out) <= 1) {
      record_error(name, paste0("parsed ", ncol(out), " column(s); expected a table with >1 columns"))
    }
  }
  on.exit({
    read_delim_errors[[path]] <<- errors
  }, add = TRUE)

  if (!file.exists(path)) {
    record_error("file", "does not exist")
    return(NULL)
  }

  lower <- tolower(path)
  out <- tryCatch({
    if (str_detect(lower, "\\.csv(\\.gz)?$")) {
      readr_call("read_csv", path)
    } else {
      readr_call("read_tsv", path)
    }
  }, error = function(e) {
    record_error(if (str_detect(lower, "\\.csv(\\.gz)?$")) "read_csv" else "read_tsv", conditionMessage(e))
    NULL
  })
  if (!is.null(out) && ncol(out) > 1) return(out)
  record_result(if (str_detect(lower, "\\.csv(\\.gz)?$")) "read_csv" else "read_tsv", out)

  out <- tryCatch(readr_call("read_csv", path), error = function(e) {
    record_error("read_csv", conditionMessage(e))
    NULL
  })
  if (!is.null(out) && ncol(out) > 1) return(out)
  record_result("read_csv", out)

  out <- tryCatch(readr_call("read_table", path), error = function(e) {
    record_error("read_table", conditionMessage(e))
    NULL
  })
  if (!is.null(out) && ncol(out) > 1) return(out)
  record_result("read_table", out)
  NULL
}

read_error_details <- function(path) {
  errors <- read_delim_errors[[path]]
  if (length(errors) == 0) return("")
  paste0(" Read attempts: ", paste(errors, collapse = " | "))
}

pick_col <- function(df, override, patterns, required = TRUE, label = "column") {
  if (!is_blank(override)) {
    if (override %in% names(df)) return(as.character(override))
    stop("Requested ", label, " not found: ", override, call. = FALSE)
  }
  lower_names <- tolower(names(df))
  for (pattern in patterns) {
    exact <- which(lower_names == tolower(pattern))
    if (length(exact) > 0) return(names(df)[exact[1]])
  }
  for (pattern in patterns) {
    hit <- which(str_detect(lower_names, regex(pattern, ignore_case = TRUE)))
    if (length(hit) > 0) return(names(df)[hit[1]])
  }
  if (required) {
    stop("Could not infer ", label, ". Available columns: ", paste(names(df), collapse = ", "), call. = FALSE)
  }
  NULL
}

score_patterns <- c("^score$", "riboscanner.*score", "prediction.*score", "predicted.*score", "probability", "^prob$", "^p_", "prediction", "pred")
uncertainty_patterns <- c("uncertainty", "std", "sd", "stderr", "se$", "variance", "var$", "entropy", "confidence", "ci")
id_patterns <- c("^id$", "orf_id", "sequence_id", "seq_id", "fasta_id", "header", "^name$", "construct", "oligo", "transcript", "tx_id")

master <- read_delim_flexible(master_path)
if (is.null(master)) {
  stop(
    "Could not read master table: ", master_path,
    " (", file_read_description(master_path), ").",
    read_error_details(master_path),
    call. = FALSE
  )
}

read_prediction_candidate <- function(path) {
  df <- read_delim_flexible(path)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  score_col <- pick_col(df, cfg$score_col, score_patterns, required = FALSE, label = "score column")
  if (is.null(score_col)) return(NULL)
  list(path = path, data = df, score_col = score_col)
}

prediction_choice <- NULL
for (candidate in prediction_candidates) {
  prediction_choice <- read_prediction_candidate(candidate)
  if (!is.null(prediction_choice)) break
}
if (is.null(prediction_choice)) {
  stop(
    "Could not find a readable prediction table with a score column. ",
    "Set --prediction-output to the RiboScanner output path. stdout: ",
    stdout_path,
    " stderr: ",
    stderr_path,
    call. = FALSE
  )
}

predictions <- prediction_choice$data
prediction_path <- prediction_choice$path
message("Prediction table selected: ", prediction_path)

message("STEP: picking ID/score columns")
message("  master cols (first 10): ", paste(head(names(master), 10), collapse = ", "))
message("  prediction cols: ", paste(names(predictions), collapse = ", "))
master_id_col <- pick_col(master, cfg$id_col, id_patterns, required = FALSE, label = "master ID column")
prediction_id_col <- pick_col(predictions, cfg$prediction_id_col, id_patterns, required = FALSE, label = "prediction ID column")
if (is.null(master_id_col)) {
  master_id_col <- names(master)[1]
  message("Falling back to first master-table column as ID: ", master_id_col)
}
if (is.null(prediction_id_col)) {
  prediction_id_col <- names(predictions)[1]
  message("Falling back to first prediction column as ID: ", prediction_id_col)
}
message("  master_id_col=", master_id_col, "  prediction_id_col=", prediction_id_col)

class_col <- pick_col(master, cfg$class_col, c("negative.*class", "^class$", "class", "label", "category", "group"), required = FALSE, label = "class column")
detected_col <- pick_col(master, cfg$detected_col, c("detected", "is_detected", "called", "call", "positive"), required = FALSE, label = "detected column")
uorf_col <- pick_col(master, cfg$uorf_col, c("uoorf", "uo_orf", "uorf", "orf_type", "orf.class", "has.*orf"), required = FALSE, label = "uORF/uoORF column")
gfp_col <- pick_col(master, cfg$gfp_col, c("gfp", "green", "fluorescence", "signal", "reporter"), required = FALSE, label = "GFP column")
score_col <- prediction_choice$score_col
uncertainty_col <- pick_col(predictions, cfg$uncertainty_col, uncertainty_patterns, required = FALSE, label = "uncertainty column")
message("  score_col=", score_col, "  uncertainty_col=", if (is.null(uncertainty_col)) "NULL" else uncertainty_col)

column_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x[[1]])
}

column_choices <- tibble(
  role = c("master_id", "prediction_id", "class", "detected", "uORF_or_uoORF", "GFP_signal", "RiboScanner_score", "RiboScanner_uncertainty"),
  column = vapply(
    list(master_id_col, prediction_id_col, class_col, detected_col, uorf_col, gfp_col, score_col, uncertainty_col),
    column_or_na,
    character(1)
  )
)
print(column_choices)

normalize_join_id <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^>", "")
  x <- str_trim(x)
  str_replace(x, "\\s.*$", "")
}

parse_boolish <- function(x) {
  y <- str_to_lower(str_trim(as.character(x)))
  dplyr::case_when(
    y %in% c("true", "t", "yes", "y", "1", "detected", "positive", "pos") ~ TRUE,
    y %in% c("false", "f", "no", "n", "0", "not_detected", "undetected", "negative", "neg") ~ FALSE,
    TRUE ~ NA
  )
}

derive_group_status <- function(df, class_col, detected_col) {
  group <- rep("all_sequences", nrow(df))
  is_negative <- rep(FALSE, nrow(df))
  if (!is.null(class_col)) {
    class_value <- str_to_lower(str_trim(as.character(df[[class_col]])))
    is_negative <- str_detect(class_value, "neg") | class_value %in% c("0", "false", "no")
    group <- ifelse(is_negative, "negative_class", as.character(df[[class_col]]))
  }
  if (!is.null(detected_col)) {
    detected <- parse_boolish(df[[detected_col]])
    group <- dplyr::case_when(
      is_negative ~ "negative_class",
      detected %in% TRUE ~ "detected",
      detected %in% FALSE ~ "not_detected",
      TRUE ~ group
    )
  }
  str_replace_all(group, "\\s+", "_")
}

derive_uorf_status <- function(df, uorf_col) {
  if (is.null(uorf_col)) return(rep("unknown_uORF_status", nrow(df)))
  value <- df[[uorf_col]]
  if (is.logical(value)) return(ifelse(value, "uORF", "no_uORF"))
  if (is.numeric(value)) return(ifelse(value > 0, "uORF", "no_uORF"))
  text <- str_to_lower(str_trim(as.character(value)))
  dplyr::case_when(
    text %in% c("false", "f", "no", "none", "0", "no_uorf", "no_orf") ~ "no_uORF",
    str_detect(text, "uo\\s*_?orf|uoorf") ~ "uoORF",
    str_detect(text, "uorf") ~ "uORF",
    text %in% c("true", "t", "yes", "1") ~ "uORF",
    is.na(text) | text == "" ~ "unknown_uORF_status",
    TRUE ~ as.character(value)
  )
}

message("STEP: mutating master2")
master_join_ids  <- normalize_join_id(master[[master_id_col]])
master_gfp_vals  <- if (!is.null(gfp_col)) readr::parse_number(as.character(master[[gfp_col]])) else NA_real_
master2 <- master %>%
  mutate(
    join_id      = master_join_ids,
    group_status = derive_group_status(., class_col, detected_col),
    uorf_status  = derive_uorf_status(., uorf_col),
    gfp_signal   = master_gfp_vals
  )

message("STEP: mutating predictions2")
pred_join_ids   <- normalize_join_id(predictions[[prediction_id_col]])
pred_scores     <- readr::parse_number(as.character(predictions[[score_col]]))
pred_uncertainty <- if (!is.null(uncertainty_col)) {
  readr::parse_number(as.character(predictions[[uncertainty_col]]))
} else {
  NA_real_
}
predictions2 <- predictions %>%
  mutate(
    join_id                  = pred_join_ids,
    riboscanner_score        = pred_scores,
    riboscanner_uncertainty  = pred_uncertainty
  )

message("STEP: joining master_kept")
lengths2 <- fasta_lengths %>% mutate(join_id = normalize_join_id(fasta_id))
master_kept <- master2 %>%
  left_join(lengths2, by = "join_id") %>%
  filter(is.na(kept_for_prediction) | kept_for_prediction)

message("STEP: inner_join master_kept x predictions2 (", nrow(master_kept), " x ", nrow(predictions2), " rows)")
analysis_df <- master_kept %>%
  inner_join(
    predictions2 %>% select(join_id, riboscanner_score, riboscanner_uncertainty, everything()),
    by = "join_id",
    suffix = c("_master", "_prediction")
  )

if (nrow(analysis_df) == 0 && nrow(master_kept) == nrow(predictions2)) {
  message("ID join produced zero rows; falling back to row-order join because kept master rows match prediction rows.")
  analysis_df <- bind_cols(master_kept, predictions2 %>% select(riboscanner_score, riboscanner_uncertainty))
} else if (nrow(analysis_df) > 0 && nrow(analysis_df) < nrow(master_kept)) {
  message(
    "ID join kept ", nrow(analysis_df), " of ", nrow(master_kept),
    " length-filtered master rows. Check column_choices.tsv if this is unexpected."
  )
}

if (nrow(analysis_df) == 0) {
  stop(
    "No rows joined between master table and predictions. ",
    "Length-filtered master rows: ", nrow(master_kept),
    "; prediction rows: ", nrow(predictions2),
    ". If RiboScanner output preserves FASTA order but not FASTA IDs, the row counts must match for fallback joining.",
    call. = FALSE
  )
}

group_counts <- analysis_df %>% count(group_status, uorf_status, name = "n") %>% arrange(group_status, uorf_status)
join_summary <- tibble(
  metric = c("master_rows", "prediction_rows", "analysis_rows_joined_after_length_filter", "fasta_sequences_removed_lt_min_length", "fasta_sequences_kept_for_prediction"),
  value = c(nrow(master), nrow(predictions), nrow(analysis_df), sum(!fasta_lengths$kept_for_prediction), sum(fasta_lengths$kept_for_prediction))
)
print(join_summary)
print(group_counts)

analysis_df <- analysis_df %>%
  group_by(group_status) %>%
  mutate(group_status_label = paste0(group_status, "\nN=", n())) %>%
  ungroup()

has_uncertainty <- any(!is.na(analysis_df$riboscanner_uncertainty))
has_gfp <- any(!is.na(analysis_df$gfp_signal))
message("STEP: plotting (has_uncertainty=", has_uncertainty, ", has_gfp=", has_gfp, ", nrow=", nrow(analysis_df), ")")

p_score <- ggplot(analysis_df, aes(x = group_status_label, y = riboscanner_score, fill = uorf_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, position = position_dodge(width = 0.8)) +
  geom_point(aes(color = uorf_status, group = uorf_status), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.55, size = 1.4) +
  labs(title = "RiboScanner score by negative/detected group and uORF status", x = "Group", y = "RiboScanner score", fill = "uORF status", color = "uORF status") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
message("STEP: save p_score")
save_plot(p_score, file.path(fig_dir, "riboscanner_score_by_group_uorf.png"), width = 11, height = 6)
message("STEP: p_score done")

score_summary <- analysis_df %>%
  group_by(group_status, uorf_status) %>%
  summarise(
    n = n(),
    mean_score = mean(riboscanner_score, na.rm = TRUE),
    sd_score = sd(riboscanner_score, na.rm = TRUE),
    mean_uncertainty = mean(riboscanner_uncertainty, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    se_score = sd_score / sqrt(n),
    se_score = ifelse(is.na(se_score), 0, se_score),
    mean_uncertainty = ifelse(is.na(mean_uncertainty), 0, mean_uncertainty)
  )

pd <- position_dodge(width = 0.55)
p_score_summary <- ggplot(score_summary, aes(group_status, mean_score, color = uorf_status)) +
  geom_linerange(aes(ymin = mean_score - mean_uncertainty, ymax = mean_score + mean_uncertainty), position = pd, size = 2.2, alpha = if (has_uncertainty) 0.35 else 0) +
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), position = pd, width = 0.18, size = 0.7) +
  geom_point(position = pd, size = 2.8) +
  geom_text(aes(label = paste0("n=", n)), position = pd, vjust = -1.0, size = 3) +
  labs(title = "Mean RiboScanner score with uncertainty and n", subtitle = "Thin error bars: standard error. Thick translucent bars: mean RiboScanner uncertainty when available.", x = "Group", y = "Mean RiboScanner score", color = "uORF status") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_score_summary, file.path(fig_dir, "riboscanner_score_summary_uncertainty.png"), width = 11, height = 6)

if (has_uncertainty) {
  p_unc <- ggplot(analysis_df, aes(x = group_status_label, y = riboscanner_uncertainty, fill = uorf_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, position = position_dodge(width = 0.8)) +
    geom_point(aes(color = uorf_status, group = uorf_status), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.55, size = 1.4) +
    labs(title = "RiboScanner uncertainty by group and uORF status", x = "Group", y = "RiboScanner uncertainty", fill = "uORF status", color = "uORF status") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_unc, file.path(fig_dir, "riboscanner_uncertainty_by_group_uorf.png"), width = 11, height = 6)
} else {
  message("No uncertainty column was found; skipping uncertainty distribution plot.")
}

if (has_gfp) {
  p_gfp <- ggplot(analysis_df, aes(x = group_status_label, y = gfp_signal, fill = uorf_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, position = position_dodge(width = 0.8)) +
    geom_point(aes(color = uorf_status, group = uorf_status), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.55, size = 1.4) +
    labs(title = "GFP signal by negative/detected group and uORF status", x = "Group", y = "GFP signal", fill = "uORF status", color = "uORF status") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_gfp, file.path(fig_dir, "gfp_signal_by_group_uorf.png"), width = 11, height = 6)

  qc_df <- analysis_df %>% filter(!is.na(gfp_signal), !is.na(riboscanner_score))
  if (nrow(qc_df) >= 3) {
    overall_spearman <- cor.test(qc_df$gfp_signal, qc_df$riboscanner_score, method = "spearman")
    overall_pearson <- cor.test(qc_df$gfp_signal, qc_df$riboscanner_score, method = "pearson")
    corr_summary <- tibble(
      comparison = c("overall_spearman", "overall_pearson"),
      n = nrow(qc_df),
      estimate = c(unname(overall_spearman$estimate), unname(overall_pearson$estimate)),
      p_value = c(overall_spearman$p.value, overall_pearson$p.value)
    )
    corr_by_group <- qc_df %>%
      group_by(group_status, uorf_status) %>%
      summarise(
        n = n(),
        spearman_rho = if (sum(!is.na(gfp_signal) & !is.na(riboscanner_score)) >= 3) cor(gfp_signal, riboscanner_score, method = "spearman", use = "complete.obs") else NA_real_,
        pearson_r = if (sum(!is.na(gfp_signal) & !is.na(riboscanner_score)) >= 3) cor(gfp_signal, riboscanner_score, method = "pearson", use = "complete.obs") else NA_real_
      ) %>%
      ungroup()
    subtitle <- paste0(
      "Spearman rho=", signif(unname(overall_spearman$estimate), 3),
      ", p=", signif(overall_spearman$p.value, 3),
      "; Pearson r=", signif(unname(overall_pearson$estimate), 3),
      ", p=", signif(overall_pearson$p.value, 3),
      "; n=", nrow(qc_df)
    )
    p_corr <- ggplot(qc_df, aes(gfp_signal, riboscanner_score, color = group_status, shape = uorf_status)) +
      geom_errorbar(aes(ymin = riboscanner_score - riboscanner_uncertainty, ymax = riboscanner_score + riboscanner_uncertainty), alpha = if (has_uncertainty) 0.18 else 0, width = 0) +
      geom_point(alpha = 0.75, size = 2) +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", size = 0.8) +
      labs(title = "QC: RiboScanner score versus GFP signal", subtitle = subtitle, x = "GFP signal", y = "RiboScanner score", color = "Group", shape = "uORF status") +
      theme_bw()
    save_plot(p_corr, file.path(fig_dir, "qc_riboscanner_score_vs_gfp.png"), width = 9, height = 7)
    readr::write_tsv(corr_summary, file.path(analysis_dir, "score_gfp_correlation_summary.tsv"))
    readr::write_tsv(corr_by_group, file.path(analysis_dir, "score_gfp_correlation_by_group.tsv"))
  } else {
    message("Fewer than three complete GFP/score rows; skipping score-versus-GFP QC.")
  }
} else {
  message("No GFP column was found; skipping GFP group plot and score-versus-GFP QC.")
}

readr::write_tsv(fasta_lengths, file.path(analysis_dir, "sequence_length_filter.tsv"))
readr::write_tsv(length_summary, file.path(analysis_dir, "sequence_length_summary.tsv"))
readr::write_tsv(column_choices, file.path(analysis_dir, "column_choices.tsv"))
readr::write_tsv(join_summary, file.path(analysis_dir, "join_summary.tsv"))
readr::write_tsv(group_counts, file.path(analysis_dir, "group_counts.tsv"))
readr::write_tsv(score_summary, file.path(analysis_dir, "score_summary_by_group_uorf.tsv"))
readr::write_tsv(analysis_df, file.path(analysis_dir, "riboscanner_master_predictions_joined.tsv"))

message("Analysis tables written to: ", analysis_dir)
message("Figures written to:         ", fig_dir)
