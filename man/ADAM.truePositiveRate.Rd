\name{ADAM.truePositiveRate}
\alias{ADAM.truePositiveRate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Profile of True Positive Rates
}
\description{
This function calculates a profile of True Positive Rates for fitness genes in at least n cell lines, with positive cases
from a reference set of essential genes.
}
\usage{
ADAM.truePositiveRate(depMat,essentialGeneSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{depMat}{
Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.
}
  \item{essentialGeneSet}{
  Reference set of predefined essential genes. This is used to define positive cases.
    }
}
\details{
This function calculates true positive rates for fitness genes in at least n cell lines (for each n). First, this function calculates the number of cell lines each gene is a fitness gene. Second, for a given number of cell lines, the set of genes that are fitness genes in at least that number of cell lines is determined. Finally, this set of genes is then compared to the reference set of essential genes to calculate a true positive rate.
}
\value{
A list of the following vectors:
\item{P}{Vector of number of genes that are depleted for a number of cell lines.}
\item{TP}{Vector of number of genes in sets of P are true positives, i.e. in the essentialGeneSet.}
\item{TPR}{TP divided by number of genes in set essentialGeneSet to give the true positive rate.}
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

pprofile<-ADAM.panessprofile(depMat=exampleDepMat)

nullmodel<-ADAM.generateNullModel(depMat=exampleDepMat,ntrials = 1000)

data(curated_BAGEL_essential)

EO<-ADAM.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )

TPR<-ADAM.truePositiveRate(exampleDepMat,curated_BAGEL_essential)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")

