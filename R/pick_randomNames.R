#' Pick n random names from the randomNames package
#'
#' The randomNames function from randomNames does not necessarily provide unique names (https://github.com/CenterForAssessment/randomNames/issues/55).
#' Hence this function has been written to fix that issue as elegant as possible.
#'
#' @param n numeric (integer) of how many unique names to return
#' @param names_to_avoid a vector of names to avoid
#' @param max_iter maximum number of while-loop iterations to avoid an infinite loop
#' @param min_name_nchar min length of name returned by randomNames
#' @param avoid_chars character vector of forbidden characters in names or NULL;
#' @param randomNames_args arguments to randomNames::randomNames
#' @param number_duplicates allow duplicate names first, then make them unique by
#' trailing numbers
#' @param number_duplicates_pad pad trailing number with brathering::pad_num_in_str
#' reject names with one or more of these characters
#' @return a vector of n unique random names
#' @export
#'
#' @examples
#' pick_randomNames(n = 1000, names_to_avoid = c("Chris", "Diana", "Leonie"), randomNames_args = list(which.names = "first"))
#' # grep("[[:digit:]]{1,2}$", names, value = T)
pick_randomNames <- function(n,
                             names_to_avoid = NULL,
                             max_iter = 500,
                             min_name_nchar = 3,
                             avoid_chars = c(" ", "-", "'", ",", ";"),
                             randomNames_args = list(which.names = "first"),
                             number_duplicates = F,
                             number_duplicates_pad = T) {

  if ("which.names" %in% names(randomNames_args)) {
    if (randomNames_args$which.names == "both") {
      if ("name.sep" %in% names(randomNames_args)) {
        sep <- randomNames_args$name.sep
      } else {
        sep <- formals(randomNames::randomNames)$name.sep
      }
      if (length(intersect(trimws(sep), avoid_chars)) > 0) {
        collision <- intersect(trimws(sep), avoid_chars)
        stop("name.sep and avoid_chars are colliding: ", collision, ". Please change name.sep or remove this character from avoid_chars.")
      }
    }
  }

  names_to_avoid <- names_to_avoid[which(!is.na(names_to_avoid))]
  names_to_avoid <- trimws(names_to_avoid)
  names <- trimws(unique(Gmisc::fastDoCall(randomNames::randomNames, args = c(randomNames_args, list(n = n)))))
  #names <- trimws(unique(randomNames::randomNames(n = n, ...)))
  names <- names[which(!names %in% names_to_avoid)]
  names <- names[sapply(names, nchar, simplify = T) >= min_name_nchar]
  if (!is.null(avoid_chars)) {
    names <- names[which(!grepl(paste(avoid_chars, collapse = "|"), names))]
  }

  iters <- 0
  while (length(names) < n && iters <= max_iter) {
    # if (number_duplicates) {
    #   temp <- trimws(Gmisc::fastDoCall(randomNames::randomNames, args = c(randomNames_args,
    #                                                                       list(n = n - length(names)))))
    #   temp <- temp[which(!temp %in% c(names_to_avoid))]
    # } else {
    #   temp <- trimws(unique(Gmisc::fastDoCall(randomNames::randomNames, args = c(randomNames_args,
    #                                                                              list(n = n - length(names))))))
    #   temp <- temp[which(!temp %in% c(names_to_avoid, names))]
    # }

    temp <- trimws(unique(Gmisc::fastDoCall(randomNames::randomNames, args = c(randomNames_args,
                                                                               list(n = n - length(names))))))
    if (number_duplicates) {
      temp <- temp[which(!temp %in% c(names_to_avoid))]
    } else {
      temp <- temp[which(!temp %in% c(names_to_avoid, names))]
    }

    temp <- temp[sapply(temp, nchar, simplify = T) >= min_name_nchar]
    if (!is.null(avoid_chars)) {
      temp <- temp[which(!grepl(paste(avoid_chars, collapse = "|"), temp))]
    }
    names <- c(names, temp)
    iters <- iters + 1
  }

  if (length(names) < n) {
    stop(paste0("max_iter reached before n unique random names could be picked. ", length(names), " unique names were reached. Either increase max_iter or reduce stringency of name selection (more genders, more ethnicities etc.)."))
  }

  if (number_duplicates) {
    names <- make.unique(names, sep = "")
    if (number_duplicates_pad) {
      names <- brathering::pad_num_in_str(names)
    }
  }
  return(names)
}

# pick_random_names <- function(
    #     n,
#     names_to_avoid = NULL,
#     max_iter = 500L,
#     min_name_nchar = 3L,
#     avoid_chars = c(" ", "-", "'", ",", ";"),
#     randomNames_args = list(which.names = "first"),
#     oversample = 1.5,
#     case_sensitive = TRUE,
#     number_duplicates = FALSE
# ) {
#   stopifnot(
#     length(n) == 1L,
#     is.numeric(n),
#     is.finite(n),
#     n >= 0,
#     n == as.integer(n),
#     length(max_iter) == 1L,
#     is.numeric(max_iter),
#     is.finite(max_iter),
#     max_iter >= 1,
#     length(min_name_nchar) == 1L,
#     is.numeric(min_name_nchar),
#     is.finite(min_name_nchar),
#     min_name_nchar >= 0,
#     length(oversample) == 1L,
#     is.numeric(oversample),
#     is.finite(oversample),
#     oversample >= 1,
#     length(case_sensitive) == 1L,
#     is.logical(case_sensitive),
#     !is.na(case_sensitive),
#     length(number_duplicates) == 1L,
#     is.logical(number_duplicates),
#     !is.na(number_duplicates)
#   )
#
#   n <- as.integer(n)
#   max_iter <- as.integer(max_iter)
#   min_name_nchar <- as.integer(min_name_nchar)
#
#   if (n == 0L) {
#     return(character())
#   }
#
#   if ("n" %in% names(randomNames_args)) {
#     stop("Do not supply `n` in `randomNames_args`.")
#   }
#
#   compare_names <- function(x) {
#     if (case_sensitive) x else tolower(x)
#   }
#
#   normalize <- function(x) {
#     if (is.null(x)) {
#       return(character())
#     }
#
#     x <- trimws(as.character(x))
#     x <- x[!is.na(x) & nzchar(x)]
#     compare_names(x)
#   }
#
#   avoided <- unique(normalize(names_to_avoid))
#   browser()
#   # Check whether the separator used for full names is forbidden.
#   if (
#     identical(randomNames_args$which.names, "both") &&
#     length(avoid_chars) > 0L
#   ) {
#     sep <- randomNames_args$name.sep
#
#     if (is.null(sep)) {
#       sep <- eval(formals(randomNames::randomNames)$name.sep)
#     }
#
#     checked_chars <- avoid_chars[
#       !is.na(avoid_chars) & nzchar(avoid_chars)
#     ]
#
#     collisions <- checked_chars[
#       vapply(
#         checked_chars,
#         function(ch) grepl(ch, sep, fixed = TRUE),
#         logical(1)
#       )
#     ]
#
#     if (length(collisions) > 0L) {
#       stop(
#         "`name.sep` contains forbidden character(s): ",
#         paste(shQuote(unique(collisions)), collapse = ", "),
#         ". Change `name.sep` or remove the character(s) from ",
#         "`avoid_chars`."
#       )
#     }
#   }
#
#   filter_candidates <- function(x) {
#     x <- trimws(as.character(x))
#
#     valid <- !is.na(x) &
#       nzchar(x) &
#       nchar(x, type = "chars") >= min_name_nchar
#
#     checked_chars <- avoid_chars[
#       !is.na(avoid_chars) & nzchar(avoid_chars)
#     ]
#
#     if (length(checked_chars) > 0L) {
#       contains_forbidden <- Reduce(
#         `|`,
#         lapply(
#           checked_chars,
#           function(ch) grepl(ch, x, fixed = TRUE)
#         ),
#         init = rep(FALSE, length(x))
#       )
#
#       valid <- valid & !contains_forbidden
#     }
#
#     x[valid]
#   }
#
#   # Add trailing numbers only after sampling is complete.
#   #
#   # Existing sampled names are reserved so, for example:
#   # c("James", "James", "James2") becomes
#   # c("James", "James3", "James2"), not
#   # c("James", "James2", "James2").
#   make_unique_numbered <- function(x, unavailable = character()) {
#     original_keys <- compare_names(x)
#     reserved <- unique(c(original_keys, unavailable))
#     used <- character()
#     result <- character(length(x))
#
#     for (i in seq_along(x)) {
#       base <- x[[i]]
#       base_key <- original_keys[[i]]
#
#       if (!(base_key %in% used)) {
#         candidate <- base
#         candidate_key <- base_key
#       } else {
#         suffix <- 2L
#
#         repeat {
#           candidate <- paste0(base, suffix)
#           candidate_key <- compare_names(candidate)
#
#           if (
#             !(candidate_key %in% reserved) &&
#             !(candidate_key %in% used)
#           ) {
#             break
#           }
#
#           suffix <- suffix + 1L
#         }
#       }
#
#       result[[i]] <- candidate
#       used <- c(used, candidate_key)
#     }
#
#     result
#   }
#
#   # Sample valid names until n is reached. Duplicates are retained.
#   sampled <- character()
#
#   for (iteration in seq_len(max_iter)) {
#     needed <- n - length(sampled)
#
#     if (needed <= 0L) {
#       break
#     }
#
#     batch_size <- max(
#       needed,
#       as.integer(ceiling(needed * oversample))
#     )
#
#     candidates <- do.call(
#       randomNames::randomNames,
#       c(randomNames_args, list(n = batch_size))
#     )
#
#     candidates <- filter_candidates(candidates)
#
#     # Exclude prohibited names, but retain duplicates in the sample.
#     candidate_keys <- compare_names(candidates)
#     candidates <- candidates[!(candidate_keys %in% avoided)]
#
#     if (length(candidates) > 0L) {
#       sampled <- c(sampled, candidates)
#     }
#   }
#
#   if (length(sampled) < n) {
#     stop(
#       sprintf(
#         paste0(
#           "Could only sample %d of %d valid names after %d iterations. ",
#           "Increase `max_iter` or relax the name-selection constraints."
#         ),
#         length(sampled),
#         n,
#         max_iter
#       )
#     )
#   }
#
#   # Remove any excess names introduced by oversampling.
#   sampled <- sampled[seq_len(n)]
#
#   # Apply uniqueness only after all n names have been sampled.
#   if (number_duplicates) {
#     sampled <- make_unique_numbered(
#       sampled,
#       unavailable = avoided
#     )
#   }
#
#   sampled
# }

