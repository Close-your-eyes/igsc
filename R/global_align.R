#' (Semi-) Global alignment with Needleman-Wunsch algorithm
#'
#' The frame of this function has been taken from the "needles"-function from CRANs NameNeedle package.
#' Modifications were necessary to allow semi-global alignment and to implement additional arguments
#' for performance optimization.
#'
#' @param seq1 character sequence 1 (string) to align with seq2
#' @param seq2 character sequence 2 (string) to align with seq1
#' @param semi_global logical whether to perform a semi-global alignment which
#' means that there is not penalty for terminal gaps
#' @param score_matrix_centering NULL or numeric; only fill the diagonal center of the score matrix
#' to increase computational speed; this is a heuristic and may omit to find the optimal
#' alignment. score_matrix_centering = 10 means that for one residue in seq1 the scores for
#' +/- 10 residues away from that position of seq2 are calculated. if NULL the whole
#' matrix is caculated
#' @param min_score_for_backtrace NULL or numeric; minimal score to initiate the backtracing;
#' intended to increase performance if 1000s of alignments are to be computed and a minimal
#' score per alignment is expected/wished for
#' @param match numeric; score for match
#' @param mismatch numeric; score for mismatch
#' @param gap numeric; score for gap
#'
#' @return list of score, aligned sequences, score matrix and direction matrix
#' @export
#'
#' @examples
global_align <- function(seq1,
                         seq2,
                         semi_global = T,
                         score_matrix_centering = NULL,
                         min_score_for_backtrace = NULL,
                         match = 1,
                         mismatch = -1,
                         gap = -1) {

  MATCH <- match
  MISMATCH <- mismatch
  GAP <- gap
  GAPCHAR <- "-"

  #if (grepl(",", seq1) || grepl(",", seq2)) {stop("comma (,) are not allowed in seq1 or seq2.")}

  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]

  ## Initialize
  # s1 = cols
  # s2 = rows
  scoreMatrix <- matrix(0, ncol=1+length(s1), nrow=1+length(s2))
  direcMatrix <- matrix("none", ncol=1+length(s1), nrow=1+length(s2))

  for (j in 1:length(s1)) {
    if (!semi_global) {
      scoreMatrix[1, j+1] <- GAP*j
    }
    direcMatrix[1, j+1] <- "left"
  }
  for (i in 1:length(s2)) {
    if (!semi_global) {
      scoreMatrix[i+1,1] <- GAP*i
    }
    direcMatrix[i+1,1] <- "up"
  }

  ## Fill
  for (i in 1:length(s2)) {

    if (is.null(score_matrix_centering)) {
      jj <- 1:length(s1)
    } else {
      jj <- ((i-score_matrix_centering):(i+score_matrix_centering))[intersect(which(((i-score_matrix_centering):(i+score_matrix_centering)) <= length(s1)), which(((i-score_matrix_centering):(i+score_matrix_centering)) > 0))]
    }

    for (j in jj) {

      ## Calculate (mis)match scores
      if (s1[j] == s2[i]) {
        diagonalScore <- scoreMatrix[i,j] + MATCH
      } else {
        diagonalScore <- scoreMatrix[i,j] + MISMATCH
      }

      ## Calculate gap scores
      if (semi_global) {
        if (j == length(s1)) {
          upScore   <- scoreMatrix[i,j+1] # 0 for gap
        } else {
          upScore   <- scoreMatrix[i,j+1] + GAP
        }
      }
      if (semi_global) {
        if (i == length(s2)) {
          leftScore <- scoreMatrix[i+1,j] # 0 for gap
        } else {
          leftScore <- scoreMatrix[i+1,j] + GAP
        }
      }

      ## Choose best score
      if (diagonalScore >= upScore) {
        if (diagonalScore >= leftScore) {
          scoreMatrix[i+1,j+1] <- diagonalScore
          direcMatrix[i+1,j+1] <- "diagonal";
        } else {
          scoreMatrix[i+1,j+1] <- leftScore
          direcMatrix[i+1,j+1] <- "left";
        }
      } else {
        if (upScore >= leftScore) {
          scoreMatrix[i+1,j+1] <- upScore
          direcMatrix[i+1,j+1] <- "up";
        } else {
          scoreMatrix[i+1,j+1] <- leftScore
          direcMatrix[i+1,j+1] <- "left";
        }
      }
    }
  }
  score <- scoreMatrix[i+1,j+1]

  if (!is.null(min_score_for_backtrace)) {
    if (score < min_score_for_backtrace) {
      return(list(score = score,
                  seq1_align = NULL,
                  seq2_align = NULL,
                  scoreMatrix = scoreMatrix,
                  direcMatrix = direcMatrix))
    }
  }

  ## backtrace
  J <- length(s1)+1
  I <- length(s2)+1
  direc <- direcMatrix[I, J]


  align1 <- align2 <- c()
  while(direc != "none") {
    if (direc == "diagonal") {
      align1 <- c(s1[J-1], align1)
      align2 <- c(s2[I-1], align2)
      I <- I-1
      J <- J-1
    } else if (direc == "left") {
      align1 <- c(s1[J-1], align1)
      align2 <- c(GAPCHAR, align2)
      J <- J-1
    } else if (direc == "up") {
      align1 <- c(GAPCHAR, align1)
      align2 <- c(s2[I-1], align2)
      I <- I-1
    } else {
      stop("Error in backtracing.")
    }
    direc <- direcMatrix[I, J]
  }

 ' align1 <- align2 <- rep(",", sum(length(s1), length(s2)))
  while(direc != "none") {
    if (direc == "diagonal") {
      align1[min(which(align1 == ","))] <- s1[J-1]
      align2[min(which(align2 == ","))] <- s2[I-1]
      I <- I-1
      J <- J-1
    } else if (direc == "left") {
      align1[min(which(align1 == ","))] <- s1[J-1]
      align2[min(which(align2 == ","))] <- GAPCHAR
      J <- J-1
    } else if (direc == "up") {
      align1[min(which(align1 == ","))] <- GAPCHAR
      align2[min(which(align2 == ","))] <- s2[I-1]
      I <- I-1
    } else {
      stop("Error in backtracing.")
    }
    direc <- direcMatrix[I, J]
  }
  align1 <- rev(align1[which(align1 != ",")])
  align2 <- rev(align1[which(align2 != ",")])'


  return(list(score = score,
              seq1_align = paste(align1, collapse=''),
              seq2_align = paste(align2, collapse=''),
              scoreMatrix = scoreMatrix,
              direcMatrix = direcMatrix))
}
