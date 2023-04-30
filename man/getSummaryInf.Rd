\name{getSummaryInf}
\alias{getSummaryInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the summary information
}
\description{
Return the summary information of each stratum
}
\usage{
getSummaryInf(rdat,
              covariate=FALSE,
              target = FALSE,
              XYmodel = "1",
              bxthre = 1e-05,
              getHeterQ = TRUE,
              onlyDR = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdat}{
a data containing the stratification information. \code{rdat} is the result of \code{Stratify}.
}
\item{covariate}{
logic value indicates whether or not adjust covariates.
}
  \item{target}{
logic value indicates whether or not calculate the target effect values fo each stratum. Target effects can only be known when the true causal effect is known. If you have real data, you may not need this argument. The default value is \code{target = FALSE}.
}
  \item{XYmodel}{
a character indicating the exposure-outcome model on which the target effect is calculated and based. Only applicable when \code{target=TRUE}.
}
  \item{bxthre}{
a threshold value for the instrument-exposure assocations. When absolute instrument-exposure association is less than the threshold value, the MR ests will not be calculated for this stratum to avoid extreme MR results.
}
  \item{getHeterQ}{
logic value indicating whether generates the heterogeneity testing results. The default value is \code{getHeterQ=TRUE}.
}
  \item{onlyDR}{
logic vlaue indicating whether only get the summary information for the doubly-ranked stratification. Default value is \code{onlyDR=FALSE}
}
}
\details{
The heterogeneity statistic can be found in Supplementary Text S2. The details of the target effects can be found in Supplementary Text S3.

If covariates are to be adjusted, they should be named as \code{C1}, \code{C2}, etc in \code{rdat}.
}
\value{
\code{getSummaryInf} returns a list consisting of the following elements: \code{Rres}, \code{RHeterQ}, \code{DRres}, \code{DRHeterQ}. \code{Rres} and \code{DRres} contain the summry information for each stratum by the residual and doubly-ranked stratification method. \code{RHeterQ} and \code{DRHeterQ} contains the heterogeneity information.


The summary information for each stratum include (all variables are stratum-specific)
\item{size}{sample size}
\item{min}{minimal exposure (or covariate) value}
\item{1q}{the 1st quartile exposure (or covariate) value}
\item{3q}{the 3rd quartile exposure (or covariate) value}
\item{max}{maximal exposure (or covariate) value}
\item{Fvalue}{F statistic value of the instrument (ie the instrument strength)}
\item{bx}{instrument-exposure associations}
\item{bxse}{standard error of \code{bx}}
\item{by}{instrument-outcome associations}
\item{byse}{standard error of \code{by}}
\item{est}{MR estimate}
\item{se}{standard error of \code{est}}
\item{target}{target effect values}

Note that if the stratificaiton is for the covariate, the values of \code{min}, \code{1q}, \code{3q} and \code{max} are for the covariate as well. If you are confused about the variable, these values are for \code{rdat$M} (generally \code{rdat$M==rdat$X}).

The heterogeneity test information includes
\item{Q statistic}{heterogeneity statistic value}
\item{df}{the degree of freedom}
\item{p-value}{the p-value}

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
RES<-getSummaryInf( rdat, target=TRUE, bxthre=1e-5, XYmodel='1',getHeterQ=TRUE)
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
