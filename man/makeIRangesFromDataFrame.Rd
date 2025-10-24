\name{makeIRangesFromDataFrame}

\alias{makeIRangesFromDataFrame}
\alias{coerce,data.frame,IRanges-method}
\alias{coerce,DataFrame,IRanges-method}

\title{Make a IRanges object from a data.frame or DataFrame}

\description{
  \code{makeIRangesFromDataFrame} takes a data-frame-like object as input
  and try to automatically find the columns that describe integer ranges.
  If successful, it returns them in an object.

  The function is also the workhorse behind the coercion method from
  data.frame (or \link[S4Vectors]{DataFrame}) to \link{IRanges}.
}

\usage{
makeIRangesFromDataFrame(df, keep.extra.columns=FALSE,
        start.field="start",
        end.field=c("end", "stop"),
        starts.in.df.are.0based=FALSE,
        na.rm=FALSE)
}

\arguments{
  \item{df}{
    A data.frame or \link[S4Vectors]{DataFrame} object. If not, then
    the function first tries to turn \code{df} into a data frame with
    \code{as.data.frame(df)}.
  }
  \item{keep.extra.columns}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the columns in \code{df} that are not used to
    form the integer ranges of the returned \link{IRanges} are returned
    as metadata columns on the object. Otherwise, they are ignored.

    Note that if \code{df} has a \code{width} column, then
    \code{makeIRangesFromDataFrame} will always ignore it.
  }
  \item{start.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the start values of the integer ranges.
    Only the first name in \code{start.field} that is found
    in \code{colnames(df)} is used.
    If no one is found, then an error is raised.
  }
  \item{end.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the end values of the integer ranges.
    Only the first name in \code{start.field} that is found
    in \code{colnames(df)} is used.
    If no one is found, then an error is raised.
  }
  \item{starts.in.df.are.0based}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the start values of the integer ranges in
    \code{df} are considered to be \emph{0-based} and are converted to
    \emph{1-based} in the returned \link{IRanges} object.
    This feature is intended to make it more convenient to handle input
    that contains data obtained from resources using the "0-based
    start" convention. A notorious example of such resource is the UCSC
    Table Browser (\url{http://genome.ucsc.edu/cgi-bin/hgTables}).
  }
  \item{na.rm}{
    \code{TRUE} or \code{FALSE} (the default).

    If \code{TRUE}, then rows in the \code{df} with missing start or end
    values (i.e. the value is \code{NA}) are ignored.
    Otherwise, they trigger an error.
  }
}

\value{
  An \link{IRanges} object.

  If \code{na.rm} is set to \code{FALSE} (the default), then the 
  returned object is guaranteed to have one element per row in the input.
  However, if \code{na.rm} is set to \code{TRUE}, then the length of the
  returned object can be less than \code{nrow(df)}.

  If \code{df} has non-automatic row names (i.e. \code{rownames(df)} is
  not \code{NULL} and is not \code{seq_len(nrow(df))}), then they will be
  used to set names on the returned \link{IRanges} object.
}

\note{
  Coercing a data.frame or \link[S4Vectors]{DataFrame} \code{df}
  to \link{IRanges} (with \code{as(df, "IRanges")}), or calling
  \code{IRanges(df)}, are both equivalent to calling
  \code{makeIRangesFromDataFrame(df, keep.extra.columns=TRUE)}.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link[GenomicRanges]{makeGRangesFromDataFrame}} in the
          \pkg{GenomicRanges} package to make a \link[GenomicRanges]{GRanges}
	  object from a data.frame or \link[S4Vectors]{DataFrame}.

    \item \link{IRanges} objects.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.
  }
}

\examples{
df <- data.frame(ID=letters[1:5],
                 locus_stART=11:15,
                 locus_STop=12:16,
                 score=1:5)

## makeIRangesFromDataFrame() tries hard to figure out what columns
## to use for the start and end values of the ranges.
makeIRangesFromDataFrame(df)

makeIRangesFromDataFrame(df, keep.extra.columns=TRUE)

as(df, "IRanges")  # equivalent to the above
}

\keyword{manip}
