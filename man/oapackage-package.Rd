\name{OApackage}
\alias{OApackage}
\docType{package}
\title{
Orthogonal Array package
}
\description{
This package creates 2-level D-optimal designs with a user specified optimization function.
}
\details{
The main function in the package is \code{\link{Doptimize}}. This function generates a single D-optimal design.
The optimization function is specified using three parameters
\eqn{latex}{\alpha_1, \alpha_2, \alpha_3}. The following
function is then optimized:

\deqn{latex}{\alpha_1 D + \alpha_2 D_s + \alpha_3 D_1}

Here \eqn{latex}{D} is the \eqn{latex}{D}-efficiency of the design.
The \eqn{latex}{D_s} and \eqn{latex}{D_1} efficiencies are described
in detail in the paper Eendebak (2015).

}
\author{
P.T. Eendebak <pieter.eendebak@gmail.com>
}
\references{
    Eendebak, P.T. and Schoen, E.D. (2015) 
    \emph{Two-level designs to estimate all main effects and
two-factor interactions},
    to be submitted
    
}
\keyword{ array }
\keyword{ design }
\keyword{ D-optimal }
\seealso{
The documententation for \code{\link{Doptimize}}. }
\examples{
    ### For an example see the file example_Doptimize.R or example_design_usage.R included in the package.    
}
