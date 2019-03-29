\name{ADAM.tradeoffEO_TPR}
\alias{ADAM.tradeoffEO_TPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate ADaM threshold
}
\description{
This function finds the minimum number of cell lines in which a gene needs to be fitness in order to be called core-fitness. This is defined as the \emph{n} providing the best trade-off between i) coverage of priori-known essential genes in the resulting set of predicted core-fitness genes, i.e. fitness in at least \emph{n} cell lines, and ii) deviance from expectation of the number of fitness genes in \emph{n} cell lines.
}
\usage{
ADAM.tradeoffEO_TPR(EO,TPR,test_set_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EO}{
Profile of empirical odds values. Computed with the \code{ADAM.empiricalOdds} function.
}
  \item{TPR}{
  Profile of True positive rates for across number of cell line. Computed with the \code{ADAM.truePositiveRate} function.
  }
  \item{test_set_name}{
  Name to give to the analysis, used for plotting titles.
  }
}
\details{
Compare and plot the log10 odds ratios with the true positive rates to find the cross over point where the true positive rate falls below the odds ratio.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{point}{Number of cell lines for which a gene needs to be a fitness gene in order to be predicted as core-fitness gene.}
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

\seealso{
\code{\link{ADAM.empiricalOdds}},
\code{\link{ADAM.truePositiveRate}}
}
\examples{
data(exampleDepMat)

pprofile<-ADAM.panessprofile(depMat=exampleDepMat)

nullmodel<-ADAM.generateNullModel(depMat=exampleDepMat,ntrials = 1000)

data(curated_BAGEL_essential)

EO<-ADAM.empiricalOdds(observedCumSum = pprofile$CUMsums,
                       simulatedCumSum =nullmodel$nullCumSUM)

TPR<-ADAM.truePositiveRate(exampleDepMat,curated_BAGEL_essential)

crossoverpoint<-ADAM.tradeoffEO_TPR(EO,TPR$TPR,
                                    test_set_name = 'curated BAGEL essential')
}
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
