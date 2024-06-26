\name{Smooth}
\alias{Smooth}

\title{
Smoothing stratification results
}
\description{
\code{Smooth} helps to smooth the stratification results based on the stratum-specifc results.
}
\usage{
Smooth(RES,
       StraMet='DR',
       Rall = NA,
       baseline = NA,
       splinestyle = 'Bspline',
       Norder = 1,
       XYmodel = '0',
       Knots = NA,
       Lambda = 0,
       random_effect = TRUE,
       getHeterQ = TRUE,
       Plot = FALSE,
       ylim_used = NA)
}
\arguments{
  \item{RES}{
the summary information data derived by \code{getSummaryInf()}. One example is \code{RES<-getSummaryInf( rdat,XYmodel='2')}
}
\item{StraMet}{
a character variable indicates which stratification method result will be used for further smoothing. Default is \code{'DR'}, representing the doubly-ranked stratification results.
}
  \item{Rall}{
a vector indicates the exposure range for smoothing and visualization. If \code{Rall} is not defined by user, the default exposure range is the mean exposure of the first stratum minus one and the mean exposure of the end stratum plus one.
}
  \item{baseline}{
a value representing the baseline exposure value for visualizing the causal effect shape. If \code{baseline} is not defined by user, the default baseline value is the mean exposure of the first stratum.
}
\item{splinestyle}{
a character variable indicates the spline style used for smoothing. Default is \code{'Bspline'}, representing the B-spline style.
}
  \item{Norder}{
a positive integer value indicates the order to be used for smoothing. The default value is \code{Norder=1}.
}
\item{XYmodel}{
a character variable indicates the index of the true underlying causal model.
}
  \item{Knots}{
either a vector or a character indicating the internal knots used for smoothing. If \code{Knots} is a vector, the vector is the internal knots. If \code{Knots} is a number character, the unform internal knots with number \code{Knots} over \code{Rall} will be used. If \code{Knots} is not defined by user, no internal knots will be considered.
}
  \item{Lambda}{
a positive value indicates the tuning parameter for smoothing roughness. Default value is \code{lambda=0}, which means no roughness penalty for smoothing.
}
  \item{random_effect}{
a logic value indicates whether to use random-effect for smoothing models. Default value is \code{random_effect=TRUE}
}
  \item{getHeterQ}{
a logic value indicates whether generates the heterogeneity testing results for the visulization results. The default value is \code{getHeterQ=TRUE}.
}
  \item{Plot}{
a logic value indicates whether print the visulization results immediately. Default is \code{Plot=FALSE}
}
\item{ylim_used}{
a vector indicates the y axis limit for visulization of h'(x).
}
}
\details{
The smoothing is based on B-spline system. Relevant papers will be appear to provide further details. The present smoothing methods for stratification results are fractional polynomial method and the piecewise linear method (see Reference below). Both of them can be achieved by \code{Smooth} (see Examples below).
}
\value{
\code{Smooth()} gives a list containing the following elements:
\item{thetahat}{the estimated parameters for the B-spline basisfunctions}
\item{var.matrix}{variance-covariance matrix of \code{thetahat}}
\item{summary}{summary information of the smoothing results, consisting of the estimates, standard errors and p-values for the B-spline basisfunctions.}
\item{p}{visulization plot for the derivatives of the causal effect shape}
\item{hp}{visulization plot for the causal effect shape}
}
\references{
\cite{Staley J R, Burgess S. Semiparametric methods for estimation of a nonlinear exposure-outcome relationship using instrumental variables with application to Mendelian randomization[J]. Genetic epidemiology, 2017, 41(4): 341-352.}
}

\examples{
dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='2' )
rdat<-Stratify(dat)
RES<-getSummaryInf( rdat,target=FALSE)

library(metafor)


#the fractional polynomial method (e.g. with degree 2)
smooth_res<-Smooth(RES,Norder=3,baseline=0)

#the piecewise linear method
cutting_values<-(RES$DRres$mean[-1] + head( RES$DRres$mean,-1) )/2
smooth_res<-Smooth(RES,Norder=1,baseline=0,Knots=cutting_values)

}

