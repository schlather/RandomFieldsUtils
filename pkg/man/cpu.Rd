\name{Cores}
\alias{cores}
\alias{cores1}
\alias{cpus}
\alias{threads}
\title{CPU information}
\description{
  Returns some information of CPU and cores on linux systems
}
\usage{
cores()
cores1()
threads()
cpus()
}
\value{
  integer.

  \describe{
    \item{cores}{returns the number of cores on linux systems and
      \code{NA} otherwise}
    \item{cores1}{returns the number of cores on linux systems and
      \code{1} otherwise}
    \item{threads}{returns the number of threads (two if cores allow for
    hyperthreading)}
    \item{cpus}{number of cpus. i.e. the product \code{cores() * threads()}
    }
  }
}
\me

\examples{
cores()
}

\keyword{utilities}
\keyword{misc}
