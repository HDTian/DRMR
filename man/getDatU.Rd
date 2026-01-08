\name{getDatU}
\alias{getDatU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get a Toy Data
}
\description{
Creat a simulated data set for stratification. \code{getDatU} contains different scenarios of the confounder effects on the outcome.
}
\usage{
getDatU(N = 10000,
        IVtype = "cont",
        UYmodel = "U1",
        XYmodel = "1",
        printRR = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
a value indicates the sample size
}
  \item{IVtype}{
 a character indicates which type of IV is used. Two IV types are available: the binary/dichotomous instrument (\code{IVtype='bi'}) and the continuous instrument (\code{IVtype='cont'}).
}
  \item{UYmodel}{
a character indicates which confounder-outcome model is used. The available chioces include: \code{'U1'}, \code{'U2'} and \code{'U3'}. See the reference paper for the model details
}
  \item{XYmodel}{
a character indicates which exposure-outcome model is used. The available chioces include: \code{'A'}, \code{'B'}, \code{'C'}, \code{'D'}. See the reference paper for the model details
}
  \item{printRR}{
Logic value indicates whether the coefficient of determination values for the instrument strength to be returned.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{getDatU} returns a data frame consisting of the instrument, the exposure, the outcome information.
}
\references{
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
dat<-getDatU(IVtype='cont', XYmodel='1',UYmodel='U1',printRR=TRUE )
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
