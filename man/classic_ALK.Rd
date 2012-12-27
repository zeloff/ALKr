\name{classic_ALK}
\alias{classic_ALK}
\title{Classic Age-Length Key}
\usage{
  classic_ALK(x, age_classes = colnames(x),
    length_classes = rownames(x))
}
\arguments{
  \item{x}{A \eqn{i \times j} matrix with the count of
  individuals of age \eqn{i} and length \eqn{j}.}

  \item{age_classes}{A vector with the name of each age
  class. Defaults to the column names of \code{x}.}

  \item{length_classes}{A vector with the name of each age
  class. Defaults to the row names of \code{x}.}
}
\value{
  An \code{ALKr} object, containing a \eqn{i \times j}
  matrix with the probability of an individual of length
  \eqn{i} having age \eqn{j}, i.e. \eqn{P(i|j)}, the
  vectors of age and length classes, and information about
  the method used to generate the key.
}
\description{
  \code{classicALK} returns an Age-Length Key calculated
  from a matrix with the count of individuals per age- and
  length-class, as described by Fridriksson (1934).
}
\examples{
data(hom)
classic_ALK(hom$otoliths[[1]])
}
\references{
  Fridriksson, A. (1934). On the calculation of
  age-distribution within a stock of cod by means of
  relatively few age determinations as a key to
  measurements on a large scale. Rapp. P.-V. CIEM, 86, 1-5.
}

