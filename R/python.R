.tcrdist_python <- new.env(parent = emptyenv())


get_tcrdist_module <- function() {
  if (!exists(
    "module",
    envir = .tcrdist_python,
    inherits = FALSE
  )) {
    python_path <- system.file(
      "python",
      package = "igsc"
    )

    if (!nzchar(python_path)) {
      stop(
        "Bundled Python code could not be found.",
        call. = FALSE
      )
    }

    module <- reticulate::import_from_path(
      module = "tcrdist_function",
      path = python_path,
      convert = TRUE
    )

    assign(
      "module",
      module,
      envir = .tcrdist_python
    )
  }

  get(
    "module",
    envir = .tcrdist_python,
    inherits = FALSE
  )
}
