\name{firth_screening}
\alias{firth_screening}
\title{ttScreening with Firth bias-reduced logistic regression + multiple testing}
\usage{
firth_screening(
  predictor_list,
  Outcome,
  iteration,
  train_frac = 0.67,
  alpha = 0.05,
  firth_count_threshold = 50,
  verbose = TRUE
)
}
\arguments{
\item{predictor_list}{Numeric matrix of dimension n x p (rows = subjects, cols = predictors).
Column names are treated as predictor IDs; if missing, they will be generated.}

\item{Outcome}{Integer or logical vector of length n with values in \{0,1\}.}

\item{iteration}{Integer, number of train/test resampling iterations (e.g., 100).}

\item{train_frac}{Numeric in (0,1); fraction of each class used for training (default 0.67).}

\item{alpha}{Numeric in (0,1) p-value cutoff used inside splits (default 0.05).}

\item{firth_count_threshold}{Integer, minimum number of "successes"
to call a predictor selected by the ttScreening path (default 50).}

\item{verbose}{Logical; if TRUE prints a short summary (default TRUE).}
}
\value{
A list with:
\itemize{
\item \code{firth_tt}: character vector of predictor IDs selected by repeated ttScreening
(count \code{>= firth_count_threshold})
\item \code{Bonferroni}: character vector selected by Bonferroni-adjusted p < 0.05 (pre-screen)
\item \code{FDR}: character vector selected by FDR-adjusted p < 0.05 (pre-screen)
\item \code{prescreen_p}: named numeric vector of pre-screen p-values
\item \code{counts}: named integer vector of ttScreening success counts
}
}
\description{
firth_screening: a two-stage approach combining ttScreening with Firth’s bias-reduced logistic regression. This method is designed to handle small sample sizes and rare outcomes by reducing bias in logistic regression coefficients.

Two-stage screening:
(1) Pre-screen each predictor with bias-reduced logistic regression (brglm2).
(2) For predictors passing pre-screen, run repeated train/test splits and count
how often the test p-value < 0.05. Also report Bonferroni and FDR sets
from the pre-screen p-values.

Dr. Abrar maintains and is the author of this functions, for questions contact him at abrar0708@gmail.com.
}
\examples{
\dontrun{
res <- firth_screening(
  predictor_list, Outcome,
  iteration = 100, train_frac = 0.67, alpha = 0.05, firth_count_threshold = 50
)
length(res$firth_tt); head(res$firth_tt)
}

}
