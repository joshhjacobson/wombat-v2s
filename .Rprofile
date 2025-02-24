source("renv/activate.R")

# TODO: remove below in final version
if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
  if ("httpgd" %in% .packages(all.available = TRUE)) {
    options(vsc.plot = FALSE)
    options(device = function(...) {
      httpgd::hgd(silent = TRUE)
      .vsc.browser(httpgd::hgd_url(), viewer = "Beside")
    })
  }
}

options(
  languageserver.formatting_style = function(options) {
    style <- styler::tidyverse_style()
    style$token$fix_quotes <- NULL
    style
  },
  lintr.linters = lintr::linters_with_defaults(
    line_length_linter = NULL,
    quotes_linter = NULL
  )
)
