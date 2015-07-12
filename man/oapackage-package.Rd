\docType{package}
\name{oapackage-package}
\alias{oapackage-package}
\alias{oapackage}
\title{Orthogonal Array Package}
\description{
This package creates 2-level D-optimal designs with a user specified optimization function.

A design with non-zero D-efficiency can be used to estimate a model with all the main effects and all the two-factor interactions. The effect hierarchy assumption suggests that main effect estimation should be given more prominence than the estimation of two-factor interactions. The optimization functions in this package allow to user to create D-optimal designs that favor the main effect estimation
over the estimation of two-factor interactions.
}
\details{
The main function in the package is \code{\link{Doptimize}}. This function generates a single D-optimal design. The optimization function is specified using three parameters
\eqn{\alpha_1, \alpha_2, \alpha_3}{\alpha_1, \alpha_2, \alpha_3}. The following
function is then optimized:

\deqn{F = \alpha_1 D + \alpha_2 D_s + \alpha_3 D_1}{F = \alpha_1 D + \alpha_2 Ds + \alpha_3 D1}

Here \eqn{D} is the \eqn{D}-efficiency of the design. The  \eqn{D_s}{Ds}- and \eqn{D_1}{D1}-efficiency are defined further below.

When specifying \eqn{\alpha_1=1}{alpha_1=1}, \eqn{\alpha_2=\alpha_3=0}{alpha_2=alpha_3=0} the function generates design that optimizes the D-efficiency. For values of \eqn{\alpha_2>0}{alpha_2>0} the main effects are given more weight. This allows the user to create designs with good estimation of the main effects, while still allowing to estimate all 2-factor interactions.

The \eqn{D_s}{Ds}- and \eqn{D_1}{D1}-efficiency are defined as follows.
For a design \eqn{D} the model matrix
\eqn{X}{X} can be split into \eqn{[I X_1 X_{2}]}{[I X_1 X_2]}
with \eqn{I} the intercept, \eqn{X_1}{X_1} the main effects and \eqn{X_2}{X_2} the second order effects.
We let \eqn{X_{02}=[I X_2]} and define
\eqn{D1=(|X_1^T X_1|)^{1/(k+1)}}, 
\eqn{Ds=(|X^T X|/|X_{02}^T X_{02} |)^{1/k}} with \eqn{k} the number of factors of the design.

More details of the method and results of the the generation of optimal designs can be found the paper Eendebak (2015).
}
\author{
P.T. Eendebak <pieter.eendebak@gmail.com>

Alan Vazquez
}
\references{
    Eendebak, P.T. and Schoen, E.D. (2015) 
    \emph{Two-level designs to estimate all main effects and
two-factor interactions},
    submitted to Technometrics
    
}
\keyword{ array }
\keyword{ design }
\keyword{ D-optimal }
\seealso{
The documententation for \code{\link{Doptimize}} and \code{\link{Defficiencies}}.

For example code see the file \code{example_Doptimize.R} or \code{example_design_usage.R} included in the tests of the package.   
}
\examples{
p = Doptimize(N=32, k=7, nrestarts=20, alpha1=1, alpha2=1, alpha3=0)
}
