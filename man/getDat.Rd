\name{getDat}
\alias{getDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get a Toy Data
}
\description{
Creat a simulated data set for stratification
}
\usage{
getDat<-function(N=10000,
                 IVtype='cont',
                 ZXmodel='A',
                 XYmodel='1',
                 printRR=FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
    a value indicates the sample size
}
\item{IVtype}{
    a character indicates which type of IV is used. Three IV types are available: the binary/dichotomous instrument (\code{IVtype='bi'}), the continuous instrument (\code{IVtype='cont'}), and the high-dimensional instruments (\code{IVtype='high-dim'})
}
\item{ZXmodel}{
    a character indicates which instrument-exposure model is used. The available chioces include: \code{'A'}, \code{'B'}, \code{'C'}, \code{'D'}, \code{'E'}, \code{'F'}, \code{'G'} and \code{'H'}. See the reference paper for the model details
}
\item{XYmodel}{
    a character indicates which exposure-outcome model is used. The available chioces include: \code{'A'}, \code{'B'}, \code{'C'}, \code{'D'}, \code{'E'}, \code{'F'}, \code{'G'} and \code{'H'}. See the reference paper for the model details
}
\item{printRR}{Logic value indicates whether the coefficient of determination values for the instrument strength to be returned.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{getDat} returns a data frame consisting of the instrument, the exposure, the outcome information.
}
\references{
\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2022). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". bioRxiv, 2022-06. (Doubly-ranked stratification)}
}
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
dat<-getDat( IVtype='cont', ZXmodel='A',XYmodel='1' ) #get a toy data
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
