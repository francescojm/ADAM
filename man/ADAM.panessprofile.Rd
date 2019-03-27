\name{ADAM.panessprofile}
\alias{ADAM.panessprofile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate profile of number of fitness genes across fixed numbers of cell lines.
}
\description{
This function calculates the numbers (and cumulative numbers) of fitness genes in \emph{n} cell lines, varying \emph{n} from 1 to the number of cell lines in the dataset in input.
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
Binary dependency matrix, rows are genes and columns are samples. 1 indicates a gene was signficantly depleted for that sample, 0 otherwise.
}
  \item{display}{
  Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
  }
  \item{main_suffix}{
  	If display=TRUE, title to give to plot of number of genes depleted in a give number of cell lines, default is 'genes depleted in at least 1 cell line'
  }
  \item{xlab}{
  	If display=TRUE, label to give to x-axis, default is 'n. cell lines'
  }
}
\details{

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{panessprof}{Number of genes that are depleted for a number of cell lines}
\item{CUMsums}{Cumulative number of genes depleted in at least x cell lines}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
F. Iorio
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	load(depMat)
pprofile <- ADAM.panessprofile(depMat = depMat)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
