# See https://github.com/extendr/rextendr/issues/398#issuecomment-3402930769

# see https://github.com/r-lib/pkgbuild/pull/157
library(tomledit)

# get R package crate name
# FIXME make this programattic
cargo_path <- "src/rust/Cargo.toml"
cargo_toml <- read_toml(cargo_path)
crate_name <- get_item(cargo_toml, c("package", "name"))

# get the pkg metadata as a list from cargo metadata
pkg_metadata <- rextendr:::read_cargo_metadata()

# get the index of the extendr package
crate_idx <- which(pkg_metadata$packages$name == crate_name)

# pkg dependencies
pkg_deps <- pkg_metadata$packages$dependencies[[crate_idx]]

# identify the dependencies that are not coming from a registry
non_reg_deps <- subset(pkg_deps, is.na(source))

# get package deps as a list
# we will need to modify the workspace / path in the Cargo.toml
pkg_deps <- get_item(cargo_toml, "dependencies")

# for each non_reg_dep:
# create a new directory src/rust/{dep}
workspace_toml <- read_toml(
  file.path(pkg_metadata$workspace_root, "Cargo.toml")
)

# extract workspace specific metadata from the manifest path
workspace_fields <- get_item(workspace_toml, "workspace")

# workspace members that are dependencies
workspace_member_deps <- pkg_metadata$packages$name[which(
  pkg_metadata$packages$name %in% names(pkg_deps)
)]

# "." is the R package workspace info
workspace_fields$members <- c(".", workspace_member_deps)

# iterate through deps and copy them to the workspace
for (.dep in non_reg_deps$path) {
  fs::dir_copy(.dep, file.path("src/rust", basename(.dep)))
}

# insert the workspace settings into the R package's Cargo.toml
new_cargo_toml <- insert_items(
  cargo_toml,
  workspace = workspace_fields
)

# if the dependency is in the workspace members
# we need to set workspace = true and remove path if possible
workspace_dep_idx <- which(names(pkg_deps) %in% workspace_fields$members)

for (idx in workspace_dep_idx) {
  .dep <- pkg_deps[[idx]]
  .dep$path <- NULL
  .dep$workspace <- TRUE
  pkg_deps[[idx]] <- .dep
}

# insert the new deps
new_cargo_toml <- insert_items(
  new_cargo_toml,
  dependencies = pkg_deps
)

# write the toml back
write_toml(new_cargo_toml, cargo_path)
