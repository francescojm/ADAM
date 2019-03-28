\name{curated_BAGEL_essential}
\alias{curated_BAGEL_essential}
\docType{data}
\title{
Referecence set of core-fitness essential genes
}
\description{
Reference set of predefined core-fitness essential genes used in [1], derived from [2] and further curated as follows. In order to avoid their status (essential/non-essential) being defined a priori, a set of high-confidence cancer driver genes from [3] were filtered out.
}
\usage{data(curated_BAGEL_essential)}
\format{
    A vector of string with a gene symbol in each entry
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
[1] Behan FM & Iorio F & Picco G et al., Prioritisation of cancer therapeutic targets using CRISPR-Cas9 screens. Nature, In press.

[2] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, et al. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015;163:1515–26.

[3] Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, Schubert M, et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell. 2016;166:740–54.
}
\examples{
data(curated_BAGEL_essential)
head(curated_BAGEL_essential)
}
\keyword{datasets}
