\name{list.celfiles}
\alias{list.celfiles}

\title{List CEL files.}
\description{
  Function used to get a list of CEL files.
}
\usage{
	list.celfiles(..., listGzipped=FALSE)
}

\arguments{
  \item{\dots}{Passed to \code{\link{list.files}}}
  \item{listGzipped}{Logical. List .CEL.gz files?}
}
\value{
  Character vector with filenames.
}

\note{
  Quite often users want to use this function to pass filenames to other
  methods. In this situations, it is safer to use the argument 'full.names=TRUE'.
}
\seealso{\code{\link{list.files}}}
\examples{
if (require(hapmapsnp5)){
  path <- system.file("celFiles", package="hapmapsnp5")

  ## only the filenames
  list.celfiles(path)

  ## the filenames with full path...
  ## very useful when genotyping samples not in the working directory
  list.celfiles(path, full.names=TRUE)
}else{
  ## this won't return anything
  ## if in the working directory there isn't any CEL
  list.celfiles(getwd())
}
}
\keyword{IO}
\keyword{utilities}
