\name{ADAM.generateNullModel}
\alias{ADAM.generateNullModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate null profile of number of fitness genes across fixed numbers of cell lines and cumulative sums.
}
\description{
This function randomly perturbs the binary dependency matrix to generate a null distribution of profiles of fitness genes across fixed number of cell lines, and corresponding null distribution of cumulative sums.
}
\usage{
ADAM.generateNullModle(depMat,ntrials=100,display=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{depMat}{
Binary dependency matrix, rows are genes and columns are samples. 1 in position \emp{[i,j]} indicates that inactivation of the \emp{i}th gene exerts a significant loss of fitness in the \emp{j}th sample, 0 otherwise.
}
  \item{ntrials}{
  Integer, default =100. How many times to randomly perturb dependency matrix to generate the null distributions.
  }
  \item{display}{
  Boolean, default is TRUE. Should bar plots of the null profiles be plotted
  }
}
\details{
For a number of trials specified in (\code{ntrials}) the inputted binary dependency matrix is randomised, keeping its column marginal sums. The profiles of fitness genes across fixed number of cell lines, and corresponding cumulative sums, are returned for each random perturbation.
}
\value{
A list with the following two matrices
\item{nullProf}{Matrix of number of fitness genes for fixed number of cell lines from. Each rows of matrix corresponds to a random trial.}
\item{nullCumSum}{Matrix of profile of cumulative number of fitness genes in fixed number of cell lines. Each row of matrix is one random trial.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Pacini & F. Iorio
}
\note{
%%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(exampleDepMat)
pprofile <- ADAM.generateNullModel(depMat = exampleDepMat,ntrials=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")

