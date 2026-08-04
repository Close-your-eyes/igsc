#' Determine Sequence Order from a Distance Matrix
#'
#' Computes an ordering of sequences using either a traveling salesperson
#' problem (TSP) solver or neighbor-joining tree clustering.
#'
#' @param distmat A distance matrix or an object coercible to a
#'   [TSP::TSP()] object. For `method = "treeline"`, it must be accepted by
#'   [DECIPHER::Treeline()].
#' @param method Character string specifying the ordering method. `"tsp"`
#'   uses a TSP solver, while `"treeline"` derives the order from clusters
#'   produced by a neighbor-joining tree.
#'
#' @return An integer vector containing the sequence indices in their
#'   calculated order.
#'
#' @export
#'
#' @examples
#' distances <- stats::dist(USArrests[1:10, ])
#' get_sequence_order(distances, method = "tsp")
#'
#' \dontrun{
#' get_sequence_order(as.matrix(distances), method = "treeline")
#' }
get_sequence_order <- function(distmat,
                               method = c("tsp", "seriate", "treeline"),
                               ...) {
  method <- rlang::arg_match(method)

  # install.packages("seriation")

  if (method == "tsp") {
    order <- distmat |>
      TSP::as.TSP() |>
      TSP::solve_TSP(control = list(rep = 10)) |>
      as.numeric()
  } else if (method == "treeline") {
    # do not prefer this one
    tree <- DECIPHER::Treeline(
      myDistMatrix = distmat,
      cutoff = 0.01,
      method = "NJ",
      showPlot = FALSE,
      type = "clusters",
      verbose = FALSE
    )
    order <- base::order(tree$cluster)
  } else {

    order <- seriation::get_order(seriation::seriate(stats::as.dist(distmat), ...))

  }

  return(order)
}
