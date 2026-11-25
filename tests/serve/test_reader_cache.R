# Remote-cache isolation and correction refresh; local file URLs, no external IO.
source("ptmlink/R/ptmlink.R")
run <- function() {
  root <- tempfile("ptmlink-cache-test-"); dir.create(root)
  prior <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
  on.exit({ unlink(root, recursive = TRUE); if (is.na(prior)) Sys.unsetenv("R_USER_CACHE_DIR") else Sys.setenv(R_USER_CACHE_DIR = prior) })
  Sys.setenv(R_USER_CACHE_DIR = file.path(root, "cache"))
  a <- file.path(root, "a"); b <- file.path(root, "b")
  publish <- function(path, value) {
    dir.create(path, showWarnings = FALSE)
    jsonlite::write_json(list(value = value), file.path(path, "build.json"), auto_unbox = TRUE)
    hash <- unname(tools::sha256sum(file.path(path, "build.json")))
    jsonlite::write_json(setNames(list(hash), "build.json"), file.path(path, "SHA256SUMS.json"), auto_unbox = TRUE)
  }
  use <- function(path) ptm_set_data_url(paste0("file://", path))
  publish(a, 1); publish(b, 2)
  use(a); stopifnot(.ptm_json("build.json")$value == 1); path_a <- .ptm_path("build.json")
  use(b); stopifnot(.ptm_json("build.json")$value == 2, .ptm_path("build.json") != path_a)
  publish(a, 3)
  use(a); stopifnot(.ptm_json("build.json")$value == 3)
  writeLines("corrupt cache", .ptm_path("build.json"))
  use(a); stopifnot(.ptm_json("build.json")$value == 3)
  unlink(.ptm_path("build.json")); writeLines("corrupt upstream", file.path(a, "build.json"))
  use(a); stopifnot(inherits(try(.ptm_json("build.json"), silent = TRUE), "try-error"))
  stopifnot(inherits(try(.ptm_path("unlisted.json"), silent = TRUE), "try-error"))
}
run()
stopifnot(system2("python3", "tests/serve/test_reader_cache.py") == 0L)
cat("R/Python source isolation, correction refresh and checksum refusal passed\n")
