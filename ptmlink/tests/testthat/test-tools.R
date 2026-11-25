# ptm_tools generates the tool descriptions from the reference: the tool functions, their arguments as
# JSON Schema, the view's settings, the row limit convention, three shapes and a JSON text.
test_that("ptm_tools describes the tool functions from the reference in three shapes", {
  ref <- ptm_reference(); want <- vapply(Filter(function(f) isTRUE(f$tool), ref$functions), function(f) f$name, character(1))
  t <- ptm_tools()
  expect_equal(vapply(t, function(x) x$name, character(1)), want); expect_match(attr(t, "instructions"), "symmetric")
  byname <- setNames(t, want)
  g <- byname$ptm_gene$inputSchema
  expect_equal(g$type, "object"); expect_equal(unclass(g$required), "gene"); expect_false(g$additionalProperties)
  expect_equal(g$properties$gene$anyOf[[1]]$type, "string"); expect_equal(g$properties$gene$anyOf[[2]]$type, "array"); expect_equal(g$properties$gene$anyOf[[2]]$items$type, "string"); expect_equal(unclass(g$properties$side$enum), c("glyco", "phospho"))
  expect_equal(names(g$properties$filters$properties), names(ptm_default_filters())); expect_equal(g$properties$filters$properties$max_i2$anyOf[[1]]$type, "number"); expect_equal(g$properties$filters$properties$max_i2$anyOf[[2]]$type, "null")   # null is "no limit"
  # defaults read as a tool caller writes them (true, false, null), never as R writes them
  all_text <- paste(unlist(t), collapse = " ")
  expect_false(grepl("default (TRUE|FALSE|NA)\\b|character\\(0\\)|list\\(\\)|ptm_default_filters\\(\\)|\\$", all_text)); expect_match(g$properties$filters$properties$include_cis$description, "(default false)", fixed = TRUE)
  expect_equal(byname$ptm_network$inputSchema$properties$top_n$anyOf[[2]]$type, "null")
  expect_equal(g$properties$limit$default, 25L); expect_match(byname$ptm_gene$description, "sorted on evidence")
  expect_null(byname$ptm_query$inputSchema$properties$limit$default)                       # ptm_query's own limit, described in the reference
  expect_equal(byname$ptm_query$inputSchema$properties$where$items$type, "array")
  expect_false("rows" %in% names(byname$ptm_explain$inputSchema$properties))              # a data.frame is not a tool argument
  expect_false("limit" %in% names(byname$ptm_status$inputSchema$properties)); expect_equal(unclass(byname$ptm_status$inputSchema$required), "gene")
  expect_equal(unclass(byname$ptm_list$inputSchema$required), "x"); expect_null(byname$ptm_cohorts$inputSchema$required)
  expect_false(any(c("ptm_set_data_dir", "ptm_tools", "ptm_plot_network", "ptm_igraph", "ptm_write_network", "ptm_default_filters") %in% want))
  a <- ptm_tools("anthropic"); expect_true("input_schema" %in% names(a[[1]])); expect_false("inputSchema" %in% names(a[[1]]))
  o <- ptm_tools("openai"); expect_equal(o[[1]]$type, "function"); expect_equal(o[[1]]$`function`$name, want[1]); expect_true("parameters" %in% names(o[[1]]$`function`))
  js <- ptm_tools(as = "json"); p <- jsonlite::fromJSON(js, simplifyVector = FALSE)
  expect_equal(length(p$tools), length(want)); expect_equal(p$tools[[which(want == "ptm_gene")]]$inputSchema$required, list("gene")); expect_match(p$instructions, "ptm_status")
  expect_match(js, '"required": \\["gene"\\]')                                            # a one-element required list stays an array
})
