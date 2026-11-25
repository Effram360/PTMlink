# inst/reference.json is held to the package: every export is described, every described argument is a
# formal, every rule has a comparator the explainer knows, and the citation reads
test_that("the reference describes exactly the exported functions with their real arguments", {
  ref <- ptm_reference()
  ns <- getNamespaceExports("ptmlink"); described <- vapply(ref$functions, function(f) f$name, character(1))
  expect_setequal(described, ns)
  for (f in ref$functions) {
    formals_ <- names(formals(get(f$name, envir = asNamespace("ptmlink")))); if (is.null(formals_)) formals_ <- character(0)
    args <- vapply(f$args, function(a) a$name, character(1))
    expect_equal(unname(args), formals_, info = f$name)
  }
  known <- c("always", "informational", "<", "<=", ">", ">=", "abs>=", "==", "is_true")
  for (r in ref$rules) { expect_true(r$comparator %in% known, info = r$id); expect_true(nzchar(r$statement), info = r$id) }
  expect_true(all(c("ranking", "direction", "perturbation", "default_view") %in% names(ref$release_contract)))
  js <- ptm_reference(as = "json"); expect_true(is.character(js) && grepl('"rules"', js))
  # the filter names on the rules are real settings
  for (r in ref$rules) if (!is.null(r$filter)) expect_true(r$filter %in% names(ptm_default_filters()), info = r$id)
})
test_that("the reference's filters block and tool conventions match the package", {
  ref <- ptm_reference(); fl <- ref$filters; defaults <- ptm_default_filters()
  expect_equal(vapply(fl, function(f) f$name, character(1)), names(defaults))
  for (f in fl) { d <- defaults[[f$name]]; expect_equal(f$default, if (is.na(d)) "NA" else if (is.logical(d)) toupper(as.character(d)) else as.character(d), info = f$name); expect_true(f$type %in% names(ptmlink:::.PTM_JSON_TYPES), info = f$name) }
  expect_equal(ref$tool_conventions$call_first, "ptm_status"); expect_equal(ref$tool_conventions$row_limit$name, "limit")
  tools <- Filter(function(f) isTRUE(f$tool), ref$functions)
  expect_true(all(vapply(tools, function(f) f$name, character(1)) %in% getNamespaceExports("ptmlink")))
  for (f in tools) for (a in f$args) if (!isFALSE(a$tool)) expect_true(a$name == "filters" || a$type %in% names(ptmlink:::.PTM_JSON_TYPES), info = paste(f$name, a$name))
})
test_that("ptm_cite reads the shipped citation file and the release's when one is set", {
  assign("release_dir", NULL, envir = ptmlink:::.ptm_env); withr::local_envvar(PTMLINK_RELEASE_DIR = "")   # no release directory left by an earlier test
  c1 <- ptm_cite(); expect_match(c1$title, "PTMlink"); expect_null(c1$authors); expect_match(c1$text, "^PTMlink"); expect_match(c1$text, "CC BY 4.0")   # the citation carries no author list
  expect_match(c1$license, "CC-BY-4.0")
  # the shipped file names the preprint: the line cites it by journal, year and DOI, then the resource's version and address
  expect_equal(c1$paper$journal, "bioRxiv"); expect_equal(c1$paper$status, "preprint"); expect_match(c1$paper$doi, "^10\\.[0-9]+/[0-9.]+$"); expect_true(is.na(c1$doi))
  expect_match(c1$text, paste0("in human tumors. bioRxiv, 2026. doi:", c1$paper$doi, ". PTMlink release v1.0.0, released "), fixed = TRUE)
  expect_match(c1$text, ", https://ptmlink.biomarkercenter.org/. Data CC BY 4.0, code MIT.", fixed = TRUE)
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d); ptm_set_release_dir(file.path(d, "release"), verify = FALSE)
  c2 <- ptm_cite(); expect_equal(c2$title, "PTMlink fixture"); expect_equal(c2$version, "0.0.0"); expect_equal(c2$release, "fixture"); expect_null(c2$authors); expect_null(c2$paper); expect_match(c2$text, "^PTMlink fixture \\(2026\\)\\. Version 0\\.0\\.0")   # a citation file without a paper: the resource alone
})
test_that("the package's own words name no fixed cohort count, no release and no never-released method", {
  expect_false(grepl("\\bten\\b", utils::packageDescription("ptmlink")$Description))
  cit <- readLines(system.file("CITATION", package = "ptmlink"))
  expect_false(any(grepl("Release v1\\.1", cit))); expect_true(any(grepl("ptm_release", cit)))
  ref <- ptm_reference()
  expect_false(grepl("ten CPTAC", ref$purpose))
  ins <- ref$tool_conventions$instructions
  expect_false(grepl("ten CPTAC|student_normal_by|q_bh|purity-adjusted", ins)); expect_true(grepl("student_normal_bh", ins))
  expect_false(any(vapply(ref$columns, function(c) c$name %in% c("q_bh", "q_bh_anchored"), logical(1))))
})

test_that("citation('ptmlink') gives the article's title, a key and no placeholder, in step with ptm_cite()", {
  assign("release_dir", NULL, envir = ptmlink:::.ptm_env); withr::local_envvar(PTMLINK_RELEASE_DIR = "")
  cit <- utils::citation("ptmlink"); b <- cit[[1]]
  title <- "PTMlink: a database resource for co-variation of glycosylation and phosphorylation in human tumors"
  expect_equal(b$title, title); expect_equal(ptm_cite()$title, title)
  expect_true(nzchar(attr(unclass(cit)[[1]], "key")))
  txt <- format(cit, style = "text"); bib <- paste(format(cit, style = "bibtex"), collapse = "\n")
  expect_false(any(grepl("<(?!https?://)[^>]+>", c(txt, bib), perl = TRUE))); expect_match(bib, "@Misc\\{ptmlink,")   # no placeholder; the address may sit in angle brackets
  expect_false(any(grepl("author", bib))); expect_match(bib, "url = \\{https://ptmlink"); expect_match(bib, paste0("doi = {", ptm_cite()$paper$doi, "}"), fixed = TRUE); expect_match(bib, "howpublished = {bioRxiv}", fixed = TRUE)
  expect_equal(attr(unclass(cit)[[1]], "textVersion"), ptm_cite()$text)
})
