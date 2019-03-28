\name{ADAM.panessprofile}
\alias{ADAM.panessprofile}
\title{
Calculate profile of number of fitness genes across fixed numbers of cell lines and cumulative sums.
}
\description{
This function calculates the numbers (and cumulative numbers) of genes whose inactivation exerts a fitness effect in \emph{n} cell lines, varying \emph{n} from 1 to the number of cell lines in the dataset in input.
}
\usage{
ADAM.panessprofile(depMat,
                  display=TRUE,
                  main_suffix='fitness genes in at least 1 cell line',
                  xlab='n. dependent cell lines')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{depMat}{
Binary dependency matrix, rows are genes and columns are samples. 1 in position \emp{[i,j]} indicates that inactivation of the \emp{i}th gene exerts a significant loss of fitness in the \emp{j}th sample, 0 otherwise.
}
  \item{display}{
  Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
  }
  \item{main_suffix}{
  	If display=TRUE, title suffix to give to plot of number of genes depleted in a give number of cell lines, default is 'genes depleted in at least 1 cell line'
  }
  \item{xlab}{
  	If display=TRUE, label to give to x-axis of the plots, default is 'n. cell lines'
  }
}
\details{
}
\value{
A list with the following two named vectors:
\item{panessprof}{Number of genes that are depleted for a number of cell lines}
\item{CUMsums}{Cumulative number of genes depleted in at least x cell lines}
}
\references{
}
\author{
C. Pacini & F. Iorio
}
\note{
%%  ~~further notes~~
}
\seealso{
}
\examples{
data(exampleDepMat)
pprofile <- ADAM.panessprofile(depMat = exampleDepMat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")

