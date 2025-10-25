\name{AtomicList-utils}

\alias{AtomicList-utils}
\alias{AtomicList_utils}

\alias{Ops,AtomicList,AtomicList-method}
\alias{Ops,CompressedAtomicList,CompressedAtomicList-method}
\alias{Ops,SimpleAtomicList,CompressedAtomicList-method}
\alias{Ops,CompressedAtomicList,SimpleAtomicList-method}
\alias{Ops,AtomicList,atomic-method}
\alias{Ops,atomic,AtomicList-method}
\alias{Ops,AtomicList,missing-method}
\alias{Ops,CompressedAtomicList,atomic-method}
\alias{Ops,atomic,CompressedAtomicList-method}

\alias{Math,AtomicList-method}
\alias{Math,CompressedAtomicList-method}
\alias{Math2,AtomicList-method}
\alias{Math2,CompressedAtomicList-method}
\alias{Complex,AtomicList-method}
\alias{Complex,CompressedAtomicList-method}

\alias{which,CompressedLogicalList-method}
\alias{which,SimpleLogicalList-method}
\alias{which,CompressedRleList-method}
\alias{which,SimpleRleList-method}
\alias{which.max,CompressedRleList-method}
\alias{which.max,RleList-method}
\alias{which.max,IntegerList-method}
\alias{which.max,NumericList-method}
\alias{which.min,CompressedRleList-method}
\alias{which.min,RleList-method}
\alias{which.min,IntegerList-method}
\alias{which.min,NumericList-method}

\alias{diff.AtomicList}
\alias{diff,CompressedAtomicList-method}

\alias{pmax,IntegerList-method}
\alias{pmax,NumericList-method}
\alias{pmax,RleList-method}
\alias{pmin,IntegerList-method}
\alias{pmin,NumericList-method}
\alias{pmin,RleList-method}
\alias{pmax.int,IntegerList-method}
\alias{pmax.int,NumericList-method}
\alias{pmax.int,RleList-method}
\alias{pmin.int,IntegerList-method}
\alias{pmin.int,NumericList-method}
\alias{pmin.int,RleList-method}

\alias{cumsum,CompressedAtomicList-method}
\alias{cumprod,CompressedAtomicList-method}
\alias{cummin,CompressedAtomicList-method}
\alias{cummax,CompressedAtomicList-method}

\alias{smoothEnds,CompressedIntegerList-method}
\alias{smoothEnds,SimpleIntegerList-method}
\alias{smoothEnds,NumericList-method}
\alias{smoothEnds,RleList-method}
\alias{runmed,CompressedIntegerList-method}
\alias{runmed,SimpleIntegerList-method}
\alias{runmed,NumericList-method}
\alias{runmed,RleList-method}
\alias{runmean,RleList-method}
\alias{runsum,RleList-method}
\alias{runwtsum,RleList-method}
\alias{runq,RleList-method}

\alias{nchar,CompressedCharacterList-method}
\alias{nchar,SimpleCharacterList-method}
\alias{nchar,CompressedRleList-method}
\alias{nchar,SimpleRleList-method}
\alias{chartr,ANY,ANY,CompressedCharacterList-method}
\alias{chartr,ANY,ANY,SimpleCharacterList-method}
\alias{chartr,ANY,ANY,CompressedRleList-method}
\alias{chartr,ANY,ANY,SimpleRleList-method}
\alias{tolower,CompressedCharacterList-method}
\alias{tolower,SimpleCharacterList-method}
\alias{tolower,CompressedRleList-method}
\alias{tolower,SimpleRleList-method}
\alias{toupper,CompressedCharacterList-method}
\alias{toupper,SimpleCharacterList-method}
\alias{toupper,CompressedRleList-method}
\alias{toupper,SimpleRleList-method}
\alias{sub,ANY,ANY,CompressedCharacterList-method}
\alias{sub,ANY,ANY,SimpleCharacterList-method}
\alias{sub,ANY,ANY,CompressedRleList-method}
\alias{sub,ANY,ANY,SimpleRleList-method}
\alias{gsub,ANY,ANY,CompressedCharacterList-method}
\alias{gsub,ANY,ANY,SimpleCharacterList-method}
\alias{gsub,ANY,ANY,CompressedRleList-method}
\alias{gsub,ANY,ANY,SimpleRleList-method}
\alias{unstrsplit,CharacterList-method}
\alias{unstrsplit,RleList-method}
\alias{paste,CompressedAtomicList-method}
\alias{startsWith,CharacterList,ANY-method}
\alias{startsWith,RleList,ANY-method}
\alias{endsWith,CharacterList,ANY-method}
\alias{endsWith,RleList,ANY-method}

\alias{rank,CompressedAtomicList-method}
\alias{order,CompressedAtomicList-method}
\alias{selfmatch,CompressedAtomicList-method}
\alias{intersect,CompressedAtomicList,CompressedAtomicList-method}

\alias{ifelse2}
\alias{ifelse2,ANY,ANY,List-method}
\alias{ifelse2,ANY,List,ANY-method}
\alias{ifelse2,List,ANY,ANY-method}
\alias{ifelse2,CompressedLogicalList,ANY,ANY-method}
\alias{ifelse2,CompressedLogicalList,List,ANY-method}
\alias{ifelse2,CompressedLogicalList,ANY,List-method}
\alias{ifelse2,CompressedLogicalList,List,ANY-method}
\alias{ifelse2,CompressedLogicalList,List,List-method}
\alias{ifelse2,SimpleLogicalList,ANY,ANY-method}
\alias{ifelse2,SimpleLogicalList,List,ANY-method}
\alias{ifelse2,SimpleLogicalList,ANY,List-method}
\alias{ifelse2,SimpleLogicalList,List,ANY-method}
\alias{ifelse2,SimpleLogicalList,List,List-method}

\title{Common operations on AtomicList objects}

\description{
  Common operations on \link{AtomicList} objects.
}

\section{Group Generics}{
  \link{AtomicList} objects have support for S4 group generic functionality
  to operate within elements across objects:
  \describe{
    \item{\code{Arith}}{\code{"+"}, \code{"-"}, \code{"*"}, \code{"^"},
      \code{"\%\%"}, \code{"\%/\%"}, \code{"/"}}
    \item{\code{Compare}}{\code{"=="}, \code{">"}, \code{"<"}, \code{"!="}, 
      \code{"<="}, \code{">="}}
    \item{\code{Logic}}{\code{"&"}, \code{"|"}}
    \item{\code{Ops}}{\code{"Arith"}, \code{"Compare"}, \code{"Logic"}}
    \item{\code{Math}}{\code{"abs"}, \code{"sign"}, \code{"sqrt"},
      \code{"ceiling"}, \code{"floor"}, \code{"trunc"}, \code{"cummax"},
      \code{"cummin"}, \code{"cumprod"}, \code{"cumsum"}, \code{"log"},
      \code{"log10"}, \code{"log2"}, \code{"log1p"}, \code{"acos"},
      \code{"acosh"}, \code{"asin"}, \code{"asinh"}, \code{"atan"},
      \code{"atanh"}, \code{"exp"}, \code{"expm1"}, \code{"cos"},
      \code{"cosh"}, \code{"sin"}, \code{"sinh"}, \code{"tan"}, \code{"tanh"},
      \code{"gamma"}, \code{"lgamma"}, \code{"digamma"}, \code{"trigamma"}}
    \item{\code{Math2}}{\code{"round"}, \code{"signif"}}
    \item{\code{Complex}}{\code{"Arg"}, \code{"Conj"}, \code{"Im"},
      \code{"Mod"}, \code{"Re"}}
  }

  See \link[methods]{S4groupGeneric} for more details.

  See \link{AtomicList_summarization} for operations from the \code{Summary}
  group generic and other summarization operations on \link{AtomicList} objects.
}

\section{Other Methods}{
  \link{AtomicList} objects also support a large number of basic methods.
  Like the group generics above, these methods perform the corresponding
  operation on each element of the list separately. The methods are: 

  \describe{
    \item{Logical}{\code{!}, \code{which}}
    \item{Numeric}{\code{diff},
                   \code{pmax}, \code{pmax.int}, \code{pmin}, \code{pmin.int}}
    \item{Running Window}{\code{smoothEnds}, \code{runmed}. \code{runmean},
                          \code{runsum}, \code{runwtsum}, \code{runq}}
    \item{Character}{\code{nchar}, \code{chartr}, \code{tolower},
                     \code{toupper}, \code{sub}, \code{gsub},
                     \code{startsWith}, \code{endsWith}}
  }

  The \code{rank} method only supports tie methods \dQuote{average},
  \dQuote{first}, \dQuote{min} and \dQuote{max}.

  Since \code{\link{ifelse}} relies on non-standard evaluation for
  arguments that need to be in the generic signature, we provide
  \code{ifelse2}, which has eager but otherwise equivalent semantics.
}

\section{Specialized Methods}{
  \describe{
    \item{\code{unstrsplit(x, sep="")}:}{}: A fast
    \code{sapply(x, paste0, collapse=sep)}. See \code{?\link{unstrsplit}}
    for the details.
  }
}

\author{P. Aboyoun}

\seealso{
  \itemize{
    \item \link{AtomicList_summarization} for operations from the
          \code{Summary} group generic and other summarization
          operations on \link{AtomicList} objects.

    \item \link{AtomicList} objects.
  }
}

\examples{
## group generics
int1 <- c(1L,2L,3L,5L,2L,8L)
int2 <- c(15L,45L,20L,1L,15L,100L,80L,5L)
col1 <- IntegerList(one = int1, int2)
2 * col1
col1 + col1
col1 > 2
}
\keyword{methods}
