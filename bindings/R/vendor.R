# Build src/rust/vendor.tar.xz, the offline crate source tree that CRAN compiles
# against. Run from bindings/R, *after* bootstrap.R.
#
# This replaces rextendr::vendor_crates(), which vendors and compresses in one
# call and so leaves nowhere to drop the files cargo never reads. Those files
# dominate the tarball and we hit CRAN's 10 MB limit for the entire package.
#
# Bootstrapping first is what keeps the Python binding out. Until bootstrap.R
# gives src/rust/Cargo.toml its own [workspace] table, cargo resolves against
# the repository root and vendors all three workspace members -- pulling in
# pyo3, numpy, ndarray, rustix, linux-raw-sys and 33 other crates that the R
# package never compiles.

manifest <- "src/rust/Cargo.toml"
vendor_dir <- "src/rust/vendor"
tarball <- "src/rust/vendor.tar.xz"

if (!any(grepl("^\\s*\\[workspace\\]", readLines(manifest)))) {
  stop(
    manifest, " still resolves against the root workspace; ",
    "run `Rscript bootstrap.R` first",
    call. = FALSE
  )
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("vendor.R needs the 'jsonlite' package", call. = FALSE)
}

# Directories at a crate's root that `cargo build` never reads. Matched at the
# root only, so a crate's own src/tests/ module is never touched.
root_drop <- c("tests", "test_data", "benches", "fuzz", ".github", "ci")

# cargo never builds a *dependency's* examples, so examples/ is dropped
# everywhere. A crate bundling its own build system could reference them; none in
# the current dependency set does. If one is added, name it here with the reason.
keep_examples <- character(0)

# Directories to drop from specific crates, beyond the shared rules. Empty: no
# dependency bundles an upstream source tree with its own website or docs.
extra_drop <- list()

message("* vendoring ", manifest)
unlink(vendor_dir, recursive = TRUE)
status <- system2(
  "cargo",
  c("vendor", "--versioned-dirs", "--manifest-path", manifest, vendor_dir),
  stdout = FALSE
)
if (status != 0L) {
  stop("cargo vendor failed", call. = FALSE)
}

message("* pruning")
removed <- 0
for (crate in list.dirs(vendor_dir, recursive = FALSE)) {
  checksums <- file.path(crate, ".cargo-checksum.json")
  if (!file.exists(checksums)) next

  # `--versioned-dirs` names each directory <crate>-<version>.
  name <- sub("-[0-9][^-]*$", "", basename(crate))
  drop <- root_drop
  if (!name %in% keep_examples) drop <- c(drop, "examples")
  prefixes <- c(paste0(drop, "/"), extra_drop[[name]])

  # Cargo verifies a directory source against .cargo-checksum.json, so pruned
  # files have to be dropped from the manifest as well as from disk.
  json <- jsonlite::fromJSON(checksums, simplifyVector = FALSE)
  doomed <- Filter(
    function(f) any(startsWith(f, prefixes)),
    names(json$files)
  )
  if (!length(doomed)) next

  paths <- file.path(crate, doomed)
  removed <- removed + sum(file.size(paths), na.rm = TRUE)
  unlink(paths)
  json$files[doomed] <- NULL
  writeLines(jsonlite::toJSON(json, auto_unbox = TRUE), checksums)
}

# Leave no empty directories behind for R CMD build to warn about.
repeat {
  dirs <- list.dirs(vendor_dir, recursive = TRUE)
  empty <- dirs[vapply(dirs, function(d) {
    length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0L
  }, logical(1))]
  if (!length(empty)) break
  unlink(empty, recursive = TRUE)
}

message(sprintf("  pruned %.1f MB", removed / 1024^2))

message("* compressing ", tarball)
unlink(tarball)
Sys.setenv(XZ_OPT = "-9e")
status <- system2("tar", c("-cJf", tarball, "-C", "src/rust", "vendor"))
if (status != 0L) {
  stop("tar failed", call. = FALSE)
}
unlink(vendor_dir, recursive = TRUE)

message(sprintf("* %s: %.2f MB", tarball, file.size(tarball) / 1024^2))
