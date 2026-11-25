# The List page's connectivity excludes same-protein pairs (python/ptmlink.py::_connectivity).
# The check itself is Python; this wrapper is what puts it in the R suite, which discovers
# test_*.R only.
stopifnot(system2("python3", "tests/serve/test_list_connectivity.py") == 0L)
cat("connectivity same-protein exclusion passed\n")
