\name{host}
\alias{hostname}
\alias{pid}

\title{System calls}
\description{
  The functions \code{hostname} and \code{pid} return
  the host name and the PID, respectively.

  \command{cpus}, \command{threads}, \command{cores}, and
  \command{cores1}
  give information about the CPU.
}
\usage{
hostname()

pid()

cpus()

threads()

cores()

cores1()

}

\details{
  If R runs on a unix platform the host name and the PID
  are returned, otherwise the empty string and naught, respectively.

  The number of cpus ususally equals the number of cores times the
  number of threads. When hyperthreading is available, threads equals 2,
  usually. 
}
\value{

  \item{hostname}{returns a string}

  \item{pid}{returns an unsigned integer}

  \item{cpus, threads, cores}{return a positive integer or \code{NA}}

  \item{cores1}{returns a positive integer; if the number of cores is
  unknown it returns 1}
}

\me


\examples{
cat("The name of your computer is '", hostname(),
    "'. Your R program has current pid ", pid(), ".\n", sep="")
}


\keyword{sysdata}
\keyword{utilities}

%  LocalWords:  pid unix Schlather url
