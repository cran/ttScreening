
\name{Rare_Screening}
\alias{Rare_Screening}
\title{Rare_Screening: resampling-based screening with limma}
\usage{
Rare_Screening(predictor_list, Outcome, iteration, cut)
}
\arguments{
\item{predictor_list}{Numeric matrix of size n x p (rows = subjects, columns = predictors).
Column names are treated as predictor IDs; if missing, they will be generated.}

\item{Outcome}{Integer or logical vector of length n with values in \{0,1\}.}

\item{iteration}{Integer, number of bootstrap iterations (e.g., 100).}

\item{cut}{Integer, selection count threshold (e.g., 70) used to define final selection.}
}
\value{
A list with:
\itemize{
\item \code{final_selection}: character vector of selected predictor IDs (selection count >= cut)
\item \code{counts}: integer vector of selection counts for all predictors (names = predictors)
\item \code{sel_mat}: p x iteration matrix of 0/1 selections per iteration
\item \code{selected}: integer indices of the predictors that met \code{cut}
}
}
\description{
Rare_Screening: a resampling-based method that repeatedly samples cases and controls with replacement and applies bayes moderated linear modeling (limma) to identify predictors that are consistently significant. The function runs a rare-event resampling procedure on user provided data and returns the
selected predictors along with iteration-wise selection counts.

Dr. Abrar maintains and is the author of this functions, for questions contact him at abrar0708@gmail.com.
}
\examples{
\dontrun{
# predictor_list: n x p matrix; Outcome: 0/1 vector of length n
res <- Rare_Screening(predictor_list, Outcome, iteration = 100, cut = 70)
head(res$final_selection)
}

}
