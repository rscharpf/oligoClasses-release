\name{oligoSet}
\alias{oligoSet}
\docType{data}
\title{An example instance of oligoSnpSet class}
\description{
	An example instance of the \code{oligoSnpSet} class
}
\usage{data(oligoSetExample)}
\source{
	Created from the simulated locusLevelData provided in this package.
}
\seealso{\code{\link{locusLevelData}}}
\examples{
	\dontrun{
		## 'oligoSetExample' created by the following
		data(locusLevelData)
		oligoSet <- new("oligoSnpSet",
				copyNumber=integerMatrix(log2(locusLevelData[["copynumber"]]/100), 100),
				call=locusLevelData[["genotypes"]],
				callProbability=locusLevelData[["crlmmConfidence"]],
				annotation=locusLevelData[["platform"]],
				genome="hg19")
		oligoSet <- oligoSet[!is.na(chromosome(oligoSet)), ]
		oligoSet <- oligoSet[chromosome(oligoSet) < 3, ]
	}
	data(oligoSetExample)
	oligoSet
}
\keyword{datasets}
