\name{ttScreening-package}
\alias{ttScreening-package}
\docType{package}
\title{
Genome-Wide DNA Methylation Sites Screening by Use of Training and Testing Samples
}
\description{
ttScreening is a screening process to filter out non-informative DNA methylation sites by applying (ordinary or robust) linear regressions to training data, and the results are further examined using testing samples. Surrogate variables are included to account for unknown factors. This package also contains two additional functions, Rare_Screening and firth_screening, for predictor screening in high-dimensional settings when the outcome is binary and is rare (binary with low prevalence). These methods were developed to handle challenges while screening predictors such as separation and small event counts and high-dimensionality, which are common in genomic and epigenomic studies.

}
\details{
\tabular{ll}{
Package: \tab ttScreening\cr
Type: \tab Package\cr
Version: \tab 1.8\cr
Date: \tab 2025-12-10\cr
License: \tab Artistic-2.0\cr
}
This package utilizes training and testing samples to filter out uninformative DNA methylation
sites for continuous and rare binary outcomes. Surrogate variables (SVs) of DNA methylation are included in the filtering process
to explain unknown factor effects.
}
\author{
Meredith Ray, Mohammad Abrar, Yu Jiang, Xin Tong, Hongmei Zhang

Maintainer: Meredith Ray <maray@memphis.edu>
}
\references{
Ray MA, Tong X, Lockett GA, Zhang H, and Karmaus WJJ. (2016) ``An Efficient Approach to Screening Epigenome-Wide Data'', BioMed Research International.

Leek JT and Storey JD. (2007) ``Capturing heterogeneity in gene expression studies by `Surrogate Variable Analysis'.'' PLoS Genetics, 3: e161.
}


\keyword{ ttScreening }

