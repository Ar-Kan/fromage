#' Create ARMA candidate models
#'
#' Creates a candidate set containing every non-seasonal `ARMA(p, q)` model
#' within the requested order limits. A candidate set is part of the
#' scientific specification of a multimodel analysis, so consider replacing
#' this convenience grid with a smaller set justified for the process being
#' studied.
#'
#' @param max_p Largest autoregressive order.
#' @param max_q Largest moving-average order.
#' @param include_white_noise Whether to include `ARMA(0, 0)`.
#'
#' @return A data frame with integer columns `p` and `q`, one row per candidate
#'   ARMA order.
#' @export
#'
#' @examples
#' arma_candidates(max_p = 1, max_q = 1)
arma_candidates <- function(max_p = 2L, max_q = 2L, include_white_noise = FALSE) {
  max_p <- validate_whole_number(max_p, "max_p", minimum = 0L)
  max_q <- validate_whole_number(max_q, "max_q", minimum = 0L)

  candidates <- expand.grid(
    p = seq.int(0L, max_p),
    q = seq.int(0L, max_q),
    KEEP.OUT.ATTRS = FALSE
  )
  candidates$p <- as.integer(candidates$p)
  candidates$q <- as.integer(candidates$q)

  if (!isTRUE(include_white_noise)) {
    candidates <- candidates[candidates$p != 0L | candidates$q != 0L, , drop = FALSE]
  }
  if (nrow(candidates) == 0L) {
    stop("The requested limits produce an empty candidate set.", call. = FALSE)
  }

  rownames(candidates) <- NULL
  candidates
}

#' Validate an ARMA candidate table
#'
#' Checks that a candidate specification is a nonempty data frame with unique,
#' non-negative, whole-number `p` and `q` values. It also discards any extra
#' columns so downstream functions receive one consistent representation.
#' This helper is called by [fit_qmm()].
#'
#' @param candidates Candidate-model data frame supplied by the user or created
#'   by [arma_candidates()].
#'
#' @return A data frame containing validated integer columns `p` and `q`, with
#'   row names removed.
#' @keywords internal
validate_candidates <- function(candidates) {
  if (!is.data.frame(candidates) || !all(c("p", "q") %in% names(candidates))) {
    stop("candidates must be a data frame with columns p and q.", call. = FALSE)
  }
  if (nrow(candidates) == 0L) {
    stop("candidates must contain at least one model.", call. = FALSE)
  }

  candidates <- candidates[, c("p", "q"), drop = FALSE]
  valid_order <- function(x) {
    is.numeric(x) && !anyNA(x) && all(is.finite(x)) && all(x >= 0) && all(x == floor(x))
  }
  if (!valid_order(candidates$p) || !valid_order(candidates$q)) {
    stop("Candidate orders p and q must be non-negative whole numbers.", call. = FALSE)
  }

  candidates$p <- as.integer(candidates$p)
  candidates$q <- as.integer(candidates$q)
  duplicated_rows <- duplicated(candidates)
  if (any(duplicated_rows)) {
    duplicated_labels <- sprintf(
      "ARMA(%d,%d)",
      candidates$p[duplicated_rows],
      candidates$q[duplicated_rows]
    )
    stop(
      "candidates contains duplicate models: ",
      paste(unique(duplicated_labels), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  rownames(candidates) <- NULL
  candidates
}
