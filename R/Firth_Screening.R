#' ttScreening with Firth bias-reduced logistic regression + multiple testing
#'
#' Two-stage screening:
#' (1) Pre-screen each predictor with bias-reduced logistic regression (brglm2).
#' (2) For predictors passing pre-screen, run repeated train/test splits and count
#'     how often the test p-value < 0.05. Also report Bonferroni and FDR sets
#'     from the pre-screen p-values.
#'
#' @param predictor_list Numeric matrix of dimension n x p (rows = subjects, cols = predictors).
#'   Column names are treated as predictor IDs; if missing, they will be generated.
#' @param Outcome Integer or logical vector of length n with values in {0,1}.
#' @param iteration Integer, number of train/test resampling iterations (e.g., 100).
#' @param train_frac Numeric in (0,1); fraction of each class used for training (default 0.67).
#' @param alpha Numeric in (0,1) p-value cutoff used inside splits (default 0.05).
#' @param firth_count_threshold Integer, minimum number of "successes"
#'   to call a predictor selected by the ttScreening path (default 50).
#' @param verbose Logical; if TRUE prints a short summary (default TRUE).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{firth_tt}: character vector of predictor IDs selected by repeated ttScreening
#'         (count \code{>= firth_count_threshold})
#'   \item \code{Bonferroni}: character vector selected by Bonferroni-adjusted p < 0.05 (pre-screen)
#'   \item \code{FDR}: character vector selected by FDR-adjusted p < 0.05 (pre-screen)
#'   \item \code{prescreen_p}: named numeric vector of pre-screen p-values
#'   \item \code{counts}: named integer vector of ttScreening success counts
#' }
#'
#' @examples
#' \dontrun{
#' res <- firth_screening(
#'   predictor_list, Outcome,
#'   iteration = 100, train_frac = 0.67, alpha = 0.05, firth_count_threshold = 50
#' )
#' length(res$firth_tt); head(res$firth_tt)
#' }
#'
#' @importFrom stats glm p.adjust
#' @importFrom brglm2 brglmFit
#' @export
firth_screening <- function(
    predictor_list,
    Outcome,
    iteration,
    train_frac = 0.67,
    alpha = 0.05,
    firth_count_threshold = 50,
    verbose = TRUE
) {
  # ---- checks (minimal, non-invasive) ----
  if (!is.matrix(predictor_list)) stop("predictor_list must be a numeric matrix (n x p).")
  n <- nrow(predictor_list); p <- ncol(predictor_list)
  if (length(Outcome) != n) stop("Outcome length must equal nrow(predictor_list).")
  if (!all(Outcome %in% c(0, 1))) stop("Outcome must be binary (0/1).")
  if (iteration < 1L) stop("iteration must be >= 1.")
  if (!(train_frac > 0 && train_frac < 1)) stop("train_frac must be in (0,1).")
  if (!(is.numeric(alpha) && alpha > 0 && alpha < 1)) stop("alpha must be in (0,1).")
  if (firth_count_threshold < 0L) stop("firth_count_threshold must be >= 0.")

  # ensure predictor IDs
  if (is.null(colnames(predictor_list))) {
    colnames(predictor_list) <- sprintf("predictor_%d", seq_len(p))
  }
  predictors <- colnames(predictor_list)

  # prepare data with requested names
  Outcome_fac <- factor(as.integer(Outcome))
  data <- cbind(Outcome = as.integer(Outcome), predictor_list)

  # ---- pre-screen with bias-reduced GLM (brglm2) ----
  prescreen_p <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    fit <- stats::glm(Outcome_fac ~ predictor_list[, j],
                      family = "binomial", method = brglm2::brglmFit)
    s <- summary(fit)$coefficients
    prescreen_p[j] <- s[2, 4]  # p-value for predictor coefficient
  }
  names(prescreen_p) <- predictors

  p_FDR <- stats::p.adjust(prescreen_p, method = "fdr")
  p_BON <- stats::p.adjust(prescreen_p, method = "bonferroni")

  # ---- repeated train/test splits for features passing p<alpha pre-screen ----
  pass <- which(prescreen_p < alpha)
  counts <- integer(p); names(counts) <- predictors

  if (length(pass) > 0L) {
    pos_idx  <- which(data[, "Outcome"] == 1L)
    ctrl_idx <- which(data[, "Outcome"] == 0L)
    npos <- length(pos_idx); nctrl <- length(ctrl_idx)
    if (npos == 0L || nctrl == 0L) stop("Both classes must be present (Outcome 0s and 1s).")

    for (k in seq_len(iteration)) {
      # match your original HPC script's reproducibility behavior
      set.seed(k)

      ntrain_pos  <- max(1L, floor(train_frac * npos))
      ntrain_ctrl <- max(1L, floor(train_frac * nctrl))

      tr_pos  <- sample(pos_idx,  ntrain_pos,  replace = FALSE)
      tr_ctrl <- sample(ctrl_idx, ntrain_ctrl, replace = FALSE)
      te_pos  <- setdiff(pos_idx,  tr_pos)
      te_ctrl <- setdiff(ctrl_idx, tr_ctrl)

      train <- rbind(data[tr_pos, , drop = FALSE], data[tr_ctrl, , drop = FALSE])
      test  <- rbind(data[te_pos, , drop = FALSE], data[te_ctrl, , drop = FALSE])

      Outcome_tr <- factor(train[, 1]); X_tr <- train[, -1, drop = FALSE]
      Outcome_te <- factor(test[, 1]);  X_te <- test[, -1, drop = FALSE]

      for (j in pass) {
        # train
        s_tr <- summary(stats::glm(Outcome_tr ~ X_tr[, j],
                                   family = "binomial", method = brglm2::brglmFit))$coefficients
        p_tr <- s_tr[2, 4]
        if (!is.na(p_tr) && p_tr < alpha) {
          # test
          s_te <- summary(stats::glm(Outcome_te ~ X_te[, j],
                                     family = "binomial", method = brglm2::brglmFit))$coefficients
          p_te <- s_te[2, 4]
          if (!is.na(p_te) && p_te < alpha) counts[j] <- counts[j] + 1L
        }
      }
    }
  }

  # ---- Final selections ----
  CpG_firth       <- predictors[counts >= firth_count_threshold]
  CpG_bonferroni  <- predictors[p_BON < 0.05]
  CpG_fdr         <- predictors[p_FDR < 0.05]

  if (verbose) {
    msg <- sprintf(
      "[firth-tt] iterations=%d, train_frac=%.2f, alpha=%.3f, count>=%d -> firth_tt=%d, FDR=%d, Bonf=%d",
      iteration, train_frac, alpha, firth_count_threshold,
      length(CpG_firth), length(CpG_fdr), length(CpG_bonferroni)
    )
    message(msg)
  }

  list(
    firth_tt    = CpG_firth,
    Bonferroni  = CpG_bonferroni,
    FDR         = CpG_fdr,
    prescreen_p = prescreen_p,
    counts      = counts
  )
}
