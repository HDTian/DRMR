\name{Stratify}
\alias{Stratify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratification
}
\description{
Do doubly-ranked and residual stratification on a data
}
\usage{
Stratify(dat, onExposure = TRUE, Ns = 10, SoP = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a dataset to be stratified.
}
  \item{onExposure}{
Logic value to determine whether stratification is done on exposure or not. If `onExposure=FALSE`, the covariate rather than the exposure will be stratified on. Default value is `onExposure==TRUE`.
}
  \item{Ns}{
a value indicates the number of stratum to be built.
}
  \item{SoP}{
a value indicates the size of pre-stratum. Default size of pre-stratum is equal to the number stratum.
}
}
\details{
`Stratify` does not require the pre-stratum or stratum has the equal size, so no need to drop individuals. Even if 'nrow(dat)/Ns' or `nrow(dat)/SoP` is not an integer, stratification still works. One can let `SoP>Ns`, which will help to increase the stability of the doubly-ranked stratification results.
}
\value{
`Stratify` returns the same data as `dat` with augmented stratificaiton information. The new columns `Rstratum`, `pre_stratum` and `DRstratum` represents the residual stratification results, the pre-stratum, and the doubly-ranked stratification results rexpectively.
}
\references{
\cite{Burgess, S., Davies, N. M., & Thompson, S. G. (2014). "Instrumental variable analysis with a nonlinear exposure-outcome relationship". Epidemiology (Cambridge, Mass.), 25(6), 877. (Residual stratification)}

\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2022). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". bioRxiv, 2022-06. (Doubly-ranked stratification)}
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
