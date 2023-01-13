#' Title
#'
#' @param pattern
#' @param subject
#' @param semi_global
#' @param match
#' @param mismatch
#' @param gap
#'
#' @return
#' @export
#'
#' @examples
global_align <- function(pattern,
                         subject,
                         semi_global = T,
                         match = 1,
                         mismatch = -1,
                         gap = -1) {

  MATCH <- match
  MISMATCH <- mismatch
  GAP <- gap
  GAPCHAR <- "-"

  patt <- strsplit(pattern, "")[[1]]
  subj <- strsplit(subject, "")[[1]]

  ## Initialize
  # subj = rows
  # patt = cols
  scoreMatrix <- matrix(0, ncol=1+length(patt), nrow=1+length(subj))
  direcMatrix <- matrix("none", ncol=1+length(patt), nrow=1+length(subj))

  for (j in 1:length(patt)) {
    if (!semi_global) {
      scoreMatrix[1, j+1] <- GAP*j
    }
    direcMatrix[1, j+1] <- "left"
  }
  for (i in 1:length(subj)) {
    if (!semi_global) {
      scoreMatrix[i+1,1] <- GAP*i
    }
    direcMatrix[i+1,1] <- "up"
  }

  ## Fill
  for (i in 1:length(subj)) {
    for (j in 1:length(patt)) {

    #for (j in ((i-10):(i+10))[intersect(which(((i-10):(i+10)) <= length(patt)), which(((i-10):(i+10)) > 0))]) {

      ## Translating from 0-based arrays and vectors to 1-based
      I <- i + 1
      J <- j + 1
      ## Calculate (mis)match scores
      if (patt[J-1] == subj[I-1]) {
        diagonalScore <- scoreMatrix[I-1, J-1] + MATCH
      } else {
        diagonalScore <- scoreMatrix[I-1, J-1] + MISMATCH
      }
      ## Calculate gap scores

      if (semi_global) {
        if (j == length(patt)) {
          upScore   <- scoreMatrix[I-1, J] # 0 for gap
        } else {
          upScore   <- scoreMatrix[I-1, J] + GAP
        }
      }
      if (semi_global) {
        if (i == length(subj)) {
          leftScore <- scoreMatrix[I, J-1] # 0 for gap
        } else {
          leftScore <- scoreMatrix[I, J-1] + GAP
        }
      }

      ## Choose best score
      if (diagonalScore >= upScore) {
        if (diagonalScore >= leftScore) {
          scoreMatrix[I, J] <- diagonalScore
          direcMatrix[I, J] <- "diagonal";
        } else {
          scoreMatrix[I, J] <- leftScore
          direcMatrix[I, J] <- "left";
        }
      } else {
        if (upScore >= leftScore) {
          scoreMatrix[I, J] <- upScore
          direcMatrix[I, J] <- "up";
        } else {
          scoreMatrix[I, J] <- leftScore
          direcMatrix[I, J] <- "left";
        }
      }
    }
  }
  theScore <- scoreMatrix[I, J]

  ## backtrace
  J <- length(patt)+1
  I <- length(subj)+1
  direc <- direcMatrix[I, J]

  '
  align1 <- align2 <- c()
  while(direc != "none") {
    if (direc == "diagonal") {
      align1 <- c(patt[J-1], align1)
      align2 <- c(subj[I-1], align2)
      I <- I-1
      J <- J-1
    } else if (direc == "left") {
      align1 <- c(patt[J-1], align1)
      align2 <- c(GAPCHAR, align2)
      J <- J-1
    } else if (direc == "up") {
      align1 <- c(GAPCHAR, align1)
      align2 <- c(subj[I-1], align2)
      I <- I-1
    } else {
      stop("Error in backtracing.")
    }
    direc <- direcMatrix[I, J]
  }'

  align1 <- align2 <- rep("X", sum(length(patt), length(subj)))
  while(direc != "none") {
    if (direc == "diagonal") {
      align1[min(which(align1 == "X"))] <- patt[J-1]
      align2[min(which(align2 == "X"))] <- subj[I-1]
      I <- I-1
      J <- J-1
    } else if (direc == "left") {
      align1[min(which(align1 == "X"))] <- patt[J-1]
      align2[min(which(align2 == "X"))] <- GAPCHAR
      J <- J-1
    } else if (direc == "up") {
      align1[min(which(align1 == "X"))] <- GAPCHAR
      align2[min(which(align2 == "X"))] <- subj[I-1]
      I <- I-1
    } else {
      stop("Error in backtracing.")
    }
    direc <- direcMatrix[I, J]
  }
  align1 <- rev(align1[which(align1 != "X")])
  align2 <- rev(align1[which(align2 != "X")])


  return(list(score = theScore,
              pattern_align = paste(align1, collapse=''),
              subject_align = paste(align2, collapse=''),
              scoreMatrix = scoreMatrix,
              direcMatrix = direcMatrix))
}
