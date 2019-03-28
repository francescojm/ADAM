\name{ADAM.randomisedepMat}
\alias{ADAM.randomisedepMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary matrix randomisation preserving column totals
}
\description{
This function takes in input a matrix and shuffles its entries column wisely. If the matrix is binary then then matrix resulting from this shuffling will have the same column marginal totals of the inpputted one.
}
\usage{
ADAM.randomisedepMat(depMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{depMat}{
  A numerical matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The matrix given in input with entries shuffled column wisely.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Pacini & F. Iorio
}
\examples{
    data(exampleDepMat)
    rnd_exampleDepMat<-ADAM.randomisedepMat(exampleDepMat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")

