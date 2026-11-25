# ptm_tools -- the package as a set of tool descriptions for an agent framework.
#
# Every description is generated from inst/reference.json: the functions marked `tool` there, their
# arguments typed in R terms and mapped to JSON Schema here, the filter settings from the reference's
# `filters` block, and the conventions (call ptm_status first, a row limit on every function that
# returns rows) from its `tool_conventions` block. The Python reader generates the same descriptions
# from the same file, and its MCP server (python/ptmlink_mcp.py) serves them; the parity test holds
# the two generators to the same names, arguments and required lists. Nothing is written by hand
# twice.

# R type words used in the reference -> JSON Schema fragments
.PTM_JSON_TYPES <- list(
  # one value or several: an agent may send "EGFR" or ["EGFR", "TFRC"], and both are read
  "character" = list(anyOf = list(list(type = "string"), list(type = "array", items = list(type = "string")))),
  "character(1)" = list(type = "string"),
  "integer(1)" = list(type = "integer"),
  "numeric(1)" = list(type = "number"),
  "logical(1)" = list(type = "boolean"),
  "list of list(column, operator, value)" = list(type = "array", items = list(type = "array", minItems = 2L, maxItems = 3L,
    items = list(anyOf = list(list(type = "string"), list(type = "number"), list(type = "boolean"), list(type = "null"),
      list(type = "array", items = list(anyOf = list(list(type = "string"), list(type = "number"),
                                                   list(type = "boolean"), list(type = "null")))))))))
# a default as a tool caller writes it (true, null), not as R does (TRUE, NA): the words come from the reference
.ptm_tool_default <- function(v, words) { s <- if (is.logical(v)) (if (isTRUE(v)) "TRUE" else "FALSE") else as.character(v); w <- words[[s]]; if (is.null(w)) s else w }
# an argument that also takes null (no limit, no cap) says so in its schema
.ptm_tool_type <- function(a) {
  s <- .PTM_JSON_TYPES[[a$type]]
  if (is.null(s)) stop("no JSON Schema mapping for the argument type ", a$type, " (", a$name, ")", call. = FALSE)
  if (isTRUE(a$nullable)) list(anyOf = list(s, list(type = "null"))) else s
}
.ptm_arg_schema <- function(a, filter_schema, words = list()) {
  s <- if (identical(a$name, "filters")) filter_schema else .ptm_tool_type(a)
  s$description <- if (!is.null(a$default)) paste0(a$meaning, " (default ", .ptm_tool_default(a$default, words), ")") else a$meaning
  if (!is.null(a$enum)) s$enum <- I(unlist(a$enum))
  s
}

#' The package as tool descriptions for an agent framework
#'
#' One description per function marked as a tool in the reference, generated from `inst/reference.json`:
#' a JSON Schema for the arguments (the view's settings spelled out), the purpose and the return shape,
#' and the conventions an agent should follow (call `ptm_status()` first; a `limit` on every function that
#' returns rows). The Python MCP server in the repository serves the same descriptions from the same file, and
#' it is that server which applies `limit`: the R functions themselves take no such argument and return every row.
#' @param format "mcp" (`name`, `description`, `inputSchema`: the Model Context Protocol tool shape), "anthropic"
#'   (`input_schema`), or "openai" (`type: function` with `parameters`)
#' @param as "list" or "json" (the text an agent framework loads)
#' @return a list of tool descriptions, or one JSON string; attribute `instructions` (or the `instructions` key of
#'   the JSON) carries the release contract an agent should read before calling anything
#' @examples
#' t <- ptm_tools()
#' vapply(t, function(x) x$name, character(1))
#' cat(ptm_tools("openai", as = "json"))
#' @export
ptm_tools <- function(format = c("mcp", "anthropic", "openai"), as = c("list", "json")) {
  format <- match.arg(format); as <- match.arg(as)
  ref <- .ptm_reference_raw(); conv <- ref$tool_conventions
  fl <- ref$filters; words <- if (is.null(conv$default_words)) list() else conv$default_words
  filter_schema <- list(type = "object", additionalProperties = FALSE,
                        properties = setNames(lapply(fl, function(f) { s <- .ptm_tool_type(f); s$description <- paste0(f$meaning, " (default ", .ptm_tool_default(f$default, words), ")"); s }), vapply(fl, function(f) f$name, character(1))))
  tools <- lapply(Filter(function(f) isTRUE(f$tool), ref$functions), function(f) {
    args <- Filter(function(a) !isFALSE(a$tool), f$args)
    props <- setNames(list(), character(0)); req <- character(0)   # named even when empty, so the JSON is {} and not []
    # an argument is required when the reference gives it no default at all; a default of null is a default
    for (a in args) { props[[a$name]] <- .ptm_arg_schema(a, filter_schema, words); if (!"default" %in% names(a)) req <- c(req, a$name) }
    if (isTRUE(f$rows) && !"limit" %in% names(props)) props[["limit"]] <- list(type = "integer", description = conv$row_limit$meaning, default = conv$row_limit$default)
    schema <- list(type = "object", properties = props, additionalProperties = FALSE)
    if (length(req)) schema$required <- I(req)
    # the note on the ranking and pair_id goes on the tools that return ranked pairs, not on every table
    desc <- paste0(f$purpose, " Returns ", f$returns$type, ".", if (isTRUE(f$ranked)) paste0(" ", conv$rows_note) else "")
    switch(format,
           mcp = list(name = f$name, description = desc, inputSchema = schema),
           anthropic = list(name = f$name, description = desc, input_schema = schema),
           openai = list(type = "function", `function` = list(name = f$name, description = desc, parameters = schema)))
  })
  if (as == "json") return(jsonlite::toJSON(list(instructions = conv$instructions, tools = tools), auto_unbox = TRUE, pretty = TRUE, null = "null"))
  attr(tools, "instructions") <- conv$instructions
  attr(tools, "format") <- format; class(tools) <- "ptm_tools"
  tools
}
