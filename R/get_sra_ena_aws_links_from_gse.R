#' Get links to fastq files from GSE accession
#'
#' Uses ffq (https://github.com/pachterlab/ffq) and processes the json file
#' returned.
#'
#' @param gse_accession GEO GSE accession
#' @param save_file json file name
#' @param save_dir directory for save_file
#' @param python_to_PATH python dir on disk to add to PATH
#'
#' @returns list
#' @export
#'
#' @examples
#' \dontrun{
#' get_sra_ena_aws_links_from_gse(gse_accession = "GSE151302")
#' }
get_sra_ena_aws_links_from_gse <- function(gse_accession,
                                           save_file = paste0(gse_accession, ".json"),
                                           save_dir = tempdir(),
                                           python_to_PATH = "~/Library/Python/3.13/bin") {

  fullpath <- path.expand(python_to_PATH)
  if (!grepl(fullpath, Sys.getenv("PATH"), fixed = T)) {
    Sys.setenv(PATH = paste(fullpath, Sys.getenv("PATH"), sep=":"))
  }

  # have ffq installed and in PATH (https://github.com/pachterlab/ffq)
  gse_accession <- toupper(trimws(gse_accession))
  path <- file.path(save_dir, save_file)
  cmd <- paste0("ffq -o ", path, " ", gse_accession)

  system(cmd, intern = T)

  js <- tidyjson::read_json(path)[[2]][[1]]
  dflong <- js |>
    tidyjson::spread_all() |>
    as.data.frame() |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    tidyr::pivot_longer(-document.id)
  jspaths <- dflong |>
    dplyr::filter(grepl("experiment$", name)) |>
    dplyr::mutate(name = gsub("\\.experiment$", "", name)) |>
    dplyr::pull(name)
  parts <- purrr::map(strsplit(jspaths, "\\."), ~c(gse_accession, .x))

  paths <- purrr::map(parts, ~paste0('[["', .x, '"]]', collapse = ""))
  names(paths) <- purrr::map_chr(parts, ~paste(.x[seq(1, length(.x), 2)], collapse = "_"))
  links <- purrr::map(paths, function(x) {
    mirrors <- names(eval(parse(text = paste0("js", x)))$files)
    names(mirrors) <- mirrors
    purrr::map(mirrors, function(mir) {
      data <- eval(parse(text = paste0("js", x)))$files[[mir]]
      purrr::map_chr(data, `[[`, "url")
    })
  })

  links <- brathering::list_invert(links)

  df_out <- purrr::map_dfr(names(links), function(i) {
    parts3 <- purrr::list_flatten(purrr::map2(parts, links[[i]], function(x,y) purrr::map(y, ~c(x[-1], i, .x))))
    parts4 <- purrr::map(parts3, ~stats::setNames(.x[seq(2, length(.x), 2)], .x[seq(1, length(.x), 2)]))
    dplyr::bind_rows(parts4) |>
      dplyr::arrange(runs) |>
      tibble::remove_rownames() |>
      tidyr::pivot_longer(cols = !!rlang::sym(i), names_to = "mirror", values_to = "url")
  })
  if ("ftp" %in% names(links)) {
    df_out <- df_out |>
      dplyr::filter(mirror == "ftp") |>
      dplyr::mutate(mirror = "https") |>
      dplyr::mutate(url = gsub("^ftp", "https", url)) |>
      dplyr::bind_rows(df_out)
  }

  sep <- brathering::find_sep(c(dflong$name, dflong$value))
  attr <- dflong |>
    dplyr::filter(grepl("attributes\\.", name) | grepl("\\.title$", name)) |>
    dplyr::filter(grepl("\\.samples\\.", name)) |>
    dplyr::mutate(name = gsub("attributes\\.", "", name)) |>
    dplyr::mutate(name = gsub("\\.", sep, name)) |>
    dplyr::mutate(name = paste0(name, sep, value))

  attr_parts <-  strsplit(attr$name, sep)
  attr_parts4 <- purrr::map(attr_parts, ~stats::setNames(.x[seq(2, length(.x), 2)], .x[seq(1, length(.x), 2)]))
  # geo_samples         samples     experiments            runs
  sample_lvl_data <- dplyr::bind_rows(attr_parts4) |>
    # get sample level data
    dplyr::filter(is.na(experiments)) |>
    tidyr::pivot_longer(cols = -c(geo_samples, samples)) |>
    tidyr::drop_na(value) |>
    tidyr::pivot_wider()

  return(list(urls = df_out, sample_info = sample_lvl_data))
}
