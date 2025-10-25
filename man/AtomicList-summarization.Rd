\name{AtomicList-summarization}

\alias{AtomicList-summarization}
\alias{AtomicList_summarization}

\alias{Summary,AtomicList-method}
\alias{Summary,CompressedRleList-method}
\alias{sum,CompressedIntegerList-method}
\alias{sum,CompressedLogicalList-method}
\alias{sum,CompressedNumericList-method}

\alias{all,CompressedRleList-method}
\alias{all,CompressedAtomicList-method}
\alias{any,CompressedAtomicList-method}
\alias{anyNA,CompressedAtomicList-method}

\alias{mean,AtomicList-method}
\alias{mean,CompressedLogicalList-method}
\alias{mean,CompressedIntegerList-method}
\alias{mean,CompressedNumericList-method}
\alias{mean,CompressedRleList-method}
\alias{var,AtomicList,missing-method}
\alias{var,AtomicList,AtomicList-method}
\alias{cov,AtomicList,AtomicList-method}
\alias{cor,AtomicList,AtomicList-method}
\alias{sd,AtomicList-method}
\alias{median,AtomicList-method}
\alias{median,CompressedAtomicList-method}
\alias{quantile,AtomicList-method}
\alias{mad,AtomicList-method}
\alias{IQR,AtomicList-method}
\alias{range,CompressedLogicalList-method}
\alias{range,CompressedIntegerList-method}
\alias{range,CompressedNumericList-method}
\alias{range,CompressedRleList-method}

\title{AtomicList summarization methods}

\description{
  Operations from the \code{Summary} group generic and other summarization
  operations on \link{AtomicList} objects.
}

\details{
  The following summarization operations are supported on \link{AtomicList}
  objects:
  \itemize{
    \item \code{any}, \code{all}, \code{min}, \code{max}, \code{range},
          \code{sum}, \code{prod}
    \item \code{mean}, \code{var}, \code{cov}, \code{cor}, \code{sd},
          \code{median}, \code{quantile}, \code{mad}, \code{IQR}
    \item \code{which.max}, \code{which.min}
  }
  See \link[methods]{S4groupGeneric} for more details.

  Note that the \code{which.min} and \code{which.max} functions have an
  extra argument, \code{global=FALSE}, which controls whether the returned
  subscripts are global (compatible with the unlisted form of the input)
  or local (compatible with the corresponding list element).
}

\seealso{
  \itemize{
    \item \link{AtomicList_utils} for common operations on
          \link{AtomicList} objects.

    \item \link{AtomicList} objects.
  }
}

\examples{
x1 <- c(1L, 2L, 3L, 5L, 2L, 8L)
x2 <- c(15L, -5L, 80L, -1L, 15L, 79L, 20L, 5L)
il <- IntegerList(one=x1, x2)

sum(il)        # equivalent to 'sapply(il, sum)'
mean(il)       # equivalent to 'sapply(il, mean)'
range(il)      # equivalent to 'sapply(il, range)'
which.max(il)  # equivalent to 'sapply(il, which.max)'
}
\keyword{methods}
