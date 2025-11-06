#' Path to the packaged peace_pre_match dataset
#'
#' Returns the absolute file path to the `peace_pre_match.rds` file
#' bundled in this package (under `inst/extdata/`). Use this if you need
#' the raw file path for reading or sharing.
#'
#' @return A length-1 character vector with a file path.
#' @export
peace_pre_match_path <- function() {
  system.file("extdata", "peace_pre_match.rds", package = "matchingguide", mustWork = TRUE)
}

#' Load the peace_pre_match dataset
#'
#' Convenience helper to read the packaged `peace_pre_match.rds` into memory.
#' @return A data frame (tibble) containing the peacekeeping example data.
#' @examples
#' d <- peace_pre_match()
#' str(d)
#' @export
peace_pre_match <- function() {
  path <- peace_pre_match_path()
  readRDS(path)
}
