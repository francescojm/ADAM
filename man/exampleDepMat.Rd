\name{exampleDepMat}
\alias{exampleDepMat}
\docType{data}
\title{
Example dependency matrix data object
}
\description{
A binary matrix summarising the fitness effect status (1 = cellular fitness reduced upon gene inactivation via CRISPR-Cas9 targeting) of all the genes (at a genome scale) across 32 human colorectal cancer cell lines, derived from [1].
}
\usage{data("exampleDepMat")}
\format{
  The format is:
 num [1:17995, 1:32] 0 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:17995] "A1BG" "A1CF" "A2M" "A2ML1" ...
  ..$ : chr [1:32] "Cell_line_1" "Cell_line_2" "Cell_line_3" "Cell_line_4" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
[1] Behan FM & Iorio F & Picco G et al., \emph{Prioritisation of cancer therapeutic targets using CRISPR-Cas9 screens}. \strong{Nature}, In press.

}
\examples{
data(exampleDepMat)
head(exampleDepMat)
}
\keyword{datasets}
