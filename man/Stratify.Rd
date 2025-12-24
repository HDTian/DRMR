\name{Stratify}
\alias{Stratify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratification
}
\description{
Do doubly-ranked or/and residual stratification on a data
}
\usage{
Stratify(dat, onExposure = TRUE, Ns = 10, SoP = NA, seed = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a dataset to be stratified. It is important to ensure that the column names match the required naming conventions for the function to work properly:

The instrument: Z

The exposure: X

The outcome: Y

(if applicable) the covariate to stratify: M

(if applicable) the covariate to adjust: C1, C2, C3, ...

Note that the stratification is only about the instrument and exposure (or covariate if you wish to stratigy on covariate), so you do not really need other variales on the stratification stage.
}
  \item{onExposure}{
Logic value to determine whether stratification is done on exposure or not. If \code{onExposure=FALSE}, the covariate rather than the exposure will be stratified on. Default value is \code{onExposure==TRUE}.
}
  \item{Ns}{
a value indicates the number of strata to be built.
}
  \item{SoP}{
a value indicates the size of pre-strata. Default size of pre-stratum is equal to the number of strata.
}
 \item{seed}{
a value of seed for reproducibility. Default value is \code{seed=NA}.
}
}
\details{
\code{Stratify} does not require the pre-stratum or stratum has the equal size, so no need to drop individuals. Even if \code{nrow(dat)/Ns} or \code{nrow(dat)/SoP} is not an integer, stratification still works. One can let \code{SoP>Ns}, which will help to increase the stability of the doubly-ranked stratification results.

Note that the \code{Stratify} function will not induce internal randomness autonomously, in the sense that the same inputted data will always return the same output, so no worries about the reproducibility even you forget to fix seed. The argument \code{seed} help to induce internal randomness when the data contains the fixed instrument or exposure values.
}
\value{
\code{Stratify} returns the same data as \code{dat} with augmented stratificaiton information. The new columns \code{Rstratum}, \code{pre_stratum} and \code{DRstratum} represents the residual stratification results, the pre-stratum, and the doubly-ranked stratification results rexpectively.
}
\references{
\cite{Burgess, S., Davies, N. M., & Thompson, S. G. (2014). "Instrumental variable analysis with a nonlinear exposure-outcome relationship". Epidemiology (Cambridge, Mass.), 25(6), 877. (Residual stratification)}

\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2022). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". bioRxiv, 2022-06. (Doubly-ranked stratification)}
}

\examples{
dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='2' )
rdat<-Stratify(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
