#' Rare_Screening: resampling-based screening with limma
#'
#' Runs your rare-event resampling procedure on real data and returns the
#' selected predictors along with iteration-wise selection counts.
#'
#' @param predictor_list Numeric matrix of size n x p (rows = subjects, columns = predictors).
#'   Column names are treated as predictor IDs; if missing, they will be generated.
#' @param Outcome Integer or logical vector of length n with values in {0,1}.
#' @param iteration Integer, number of bootstrap iterations (e.g., 100).
#' @param cut Integer, selection count threshold (e.g., 70) used to define final selection.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{final_selection}: character vector of selected predictor IDs (selection count >= cut)
#'   \item \code{counts}: integer vector of selection counts for all predictors (names = predictors)
#'   \item \code{sel_mat}: p x iteration matrix of 0/1 selections per iteration
#'   \item \code{selected}: integer indices of the predictors that met \code{cut}
#' }
#'
#' @examples
#' \dontrun{
#' # predictor_list: n x p matrix; Outcome: 0/1 vector of length n
#' res <- Rare_Screening(predictor_list, Outcome, iteration = 100, cut = 70)
#' head(res$final_selection)
#' }
#'
#' @importFrom stats model.matrix
#' @importFrom limma lmFit eBayes
#' @export
Rare_Screening <- function(predictor_list, Outcome, iteration, cut) {
  # ---- basic checks (no name changes beyond your request) ----
  if (!is.matrix(predictor_list)) stop("predictor_list must be a numeric matrix (n x p).")
  n <- nrow(predictor_list); p <- ncol(predictor_list)
  if (length(Outcome) != n) stop("Outcome length must equal nrow(predictor_list).")
  if (!all(Outcome %in% c(0, 1))) stop("Outcome must be binary (0/1).")
  if (iteration < 1L) stop("iteration must be >= 1.")
  if (cut < 0L) stop("cut must be >= 0.")

  # ensure predictor IDs
  if (is.null(colnames(predictor_list))) {
    colnames(predictor_list) <- sprintf("predictor_%d", seq_len(p))
  }
  predictors <- colnames(predictor_list)

  # assemble data frame with requested names
  Outcome <- as.integer(Outcome)
  data <- cbind(Outcome = Outcome, predictor_list)

  data_positive <- data[data[, "Outcome"] == 1L, , drop = FALSE]
  data_control  <- data[data[, "Outcome"] == 0L, , drop = FALSE]

  npositive <- nrow(data_positive)
  ncontrol  <- nrow(data_control)
  if (npositive == 0L || ncontrol == 0L) {
    stop("Both classes must be present in Outcome (need 0s and 1s).")
  }

  # storage: keep your requested object names
  sel_mat <- matrix(0L, nrow = p, ncol = iteration, dimnames = list(predictors, NULL))

  # fixed per-iteration alpha, matching your original logic
  alpha <- 0.05

  # ---- bootstrap iterations (logic preserved) ----
  for (k in seq_len(iteration)) {
    # reproducibility behavior matches your original script:
    set.seed(k)

    # sample all cases with replacement, and the same number of controls with replacement
    sample_positive <- sample.int(n = npositive, size = npositive, replace = TRUE)
    data_positive_sample <- data_positive[sample_positive, , drop = FALSE]

    sample_control <- sample.int(n = ncontrol, size = npositive, replace = TRUE)
    data_control_sample <- data_control[sample_control, , drop = FALSE]

    finaldata1 <- rbind(data_positive_sample, data_control_sample)

    Outcome_sample1 <- as.factor(finaldata1[, 1])
    predictor_list_sample1 <- finaldata1[, -1, drop = FALSE]

    # limma expects features in rows → transpose predictor matrix
    mod1 <- stats::model.matrix(~ Outcome_sample1)
    fit1 <- limma::eBayes(limma::lmFit(t(predictor_list_sample1), mod1, method = "ls"))$p.value

    # keep the variable name you requested: col_ast
    col_ast <- grep("Outcome_sample1", colnames(fit1), value = TRUE)
    if (length(col_ast) < 1L) stop("Could not locate the Outcome coefficient column in limma p-values.")
    if (length(col_ast) > 1L) col_ast <- col_ast[1L]

    fit1_select <- as.integer(fit1[, col_ast] < alpha)

    sel_mat[, k] <- fit1_select
  }

  counts <- rowSums(sel_mat)
  selected <- which(counts >= cut)
  final_selection <- predictors[selected]

  # return object with the exact names you requested
  list(
    final_selection = final_selection,
    counts = counts,
    sel_mat = sel_mat,
    selected = selected
  )
}
