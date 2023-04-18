\name{getGRstats}
\alias{getGRstats}
\title{
Get Gelman-Rubin statistics
}
\description{
Get the Gelman-Rubin uniformity statistics for each stratum. This is used to check the degree of coarsenness when the exposure is coarsened.
}
\usage{
getGRstats(rdat, Nc = 2, roundnum = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdat}{
a data containing the stratification information. \code{rdat} is the result of \code{Stratify}.
}
  \item{Nc}{
integer value indicates how many chain used for calculating the GR statistic. No absolutely optimal choice. Default value is \code{Nc=2}
}
  \item{roundnum}{
a digit indicates how many decimal places the result shoul be retained
}
}
\details{
See supplementary Text S1 of the original paper for more details.
}
\value{
\code{getGRstats} gives the GR statistic values for each stratum. Small values (<1.02, the heuristic threshold) indicates a good degree of the coarsenness.
}
\references{
\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2022). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". bioRxiv, 2022-06.}
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
dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='1' ) #get a toy data
rdat<-Stratify(dat)  #Do stratification on the data
getGRstats(rdat,Nc=2,roundnum=3)
getGRstats(rdat,Nc=5,roundnum=3)
getGRstats(rdat,Nc=10,roundnum=3)
getGRstats(rdat,Nc=100,roundnum=3)
}

