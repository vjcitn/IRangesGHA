### =========================================================================
### Summarization methods for CompressedAtomicList objects
### -------------------------------------------------------------------------
###
### Summarization methods:
###   - anyNA()
###   - 'Summary' group: any(), all(), min(), max(), range(), sum(), prod()
###   - mean()
###   - median()
###   - which.min(), which.max()


setMethod("any", "CompressedAtomicList", function(x, na.rm = FALSE) {
              stopifnot(isTRUEorFALSE(na.rm))
              ans <- sum(x, na.rm=TRUE) > 0L
              if (!na.rm) {
                  ans[!ans & any(is.na(x), na.rm=TRUE)] <- NA
              }
              ans
          })

setMethod("all", "CompressedAtomicList", function(x, na.rm = FALSE) {
              stopifnot(isTRUEorFALSE(na.rm))
              ans <- !any(!x, na.rm=TRUE)
              if (!na.rm) {
                  ans[ans & any(is.na(x), na.rm=TRUE)] <- NA
              }
              ans
          })

setMethod("anyNA", "CompressedAtomicList", function(x, recursive=FALSE) {
    callNextMethod(x, recursive=FALSE) ## recursion will just slow us down
})

setCompressedListSummaryMethod <- function(fun, where=topenv(parent.frame()))
{
    setCompressedNumericalListMethod(fun, function(x, na.rm = FALSE) {
        stopifnot(isTRUEorFALSE(na.rm))
        .Call2(C_fun, x, na.rm, PACKAGE="IRanges")
    }, where)
}

setCompressedListSummaryMethod("sum")
setCompressedListSummaryMethod("prod")
setCompressedListSummaryMethod("min")
setCompressedListSummaryMethod("max")

setMethods("range",
           list("CompressedLogicalList",
                "CompressedIntegerList",
                "CompressedNumericList",
                "CompressedRleList"),
           function(x, na.rm=FALSE) {
               stopifnot(isTRUEorFALSE(na.rm))
               cbind(min(x, na.rm=na.rm), max(x, na.rm=na.rm))
           })

setMethod("Summary", "CompressedRleList",
          function(x, ..., na.rm = FALSE) {
            toViewFun <- list(max = viewMaxs, min = viewMins, sum = viewSums)
            if (!is.null(viewFun <- toViewFun[[.Generic]])) {
              ans <- viewFun(as(x, "RleViews"), na.rm = na.rm)
              names(ans) <- names(x)
              ans
            } else if (.Generic %in% c("any", "all"))
                callNextMethod()
            else sapply(x, .Generic, na.rm = na.rm)
          })

setMethod("all", "CompressedRleList", function(x, ..., na.rm = FALSE) {
  args <- list(...)
  if (length(args) > 0L)
    stop("Only a single argument in '...' is supported for now")
  if (!isTRUEorFALSE(na.rm))
    stop("'na.rm' must be TRUE or FALSE")
  rv <- runValue(x)
  if (na.rm)
    rv <- rv[!is.na(rv)]
  rv_eltNROWS <- elementNROWS(rv)
  ans <- rv_eltNROWS == 0L
  singletons <- rv_eltNROWS == 1L
  ans[singletons] <- unlist(rv, use.names = FALSE)[singletons[togroup(PartitioningByWidth(rv))]]
  ans
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Numerical methods
###

setCompressedListWhichSummaryMethod <-
    function(fun, where=topenv(parent.frame()))
    {
        def <- function(x, global = FALSE) {
            stopifnot(isTRUEorFALSE(global))
            ans <- .Call2(C_fun, x, PACKAGE="IRanges")
            if (global) {
                ans <- toglobal(ans, x)
            }
            ans
        }
        setCompressedNumericalListMethod(fun, def, where)
    }
setCompressedListWhichSummaryMethod("which.min")
setCompressedListWhichSummaryMethod("which.max")

setMethods("mean",
           list("CompressedLogicalList",
                "CompressedIntegerList",
                "CompressedNumericList",
                "CompressedRleList"),
           function(x, trim = 0, na.rm = FALSE) {
               stopifnot(isTRUEorFALSE(na.rm))
               stopifnot(isSingleNumber(trim))
               if (trim > 0) {
                   return(callNextMethod())
               }
               x_eltNROWS <- if (na.rm) sum(!is.na(x)) else elementNROWS(x)
               sum(x, na.rm=na.rm) / x_eltNROWS
           })

setMethod("median", "CompressedAtomicList", function(x, na.rm=FALSE) {
    stopifnot(isTRUEorFALSE(na.rm))
    sx <- sort(x)
    n <- lengths(sx)
    half <- (n + 1L)%/%2L
    even <- n%%2L != 1L
    ind <- IRanges(half, width=1L+even)
    NAs <- half == 0L
    ind <- relist(ind[!NAs], PartitioningByWidth(as.integer(!NAs)))
    ## ind <- as(half, "IntegerList")
    ## ind[even] <- ind[even] + as(0:1, "IntegerList")
    ans <- mean(sx[ind])
    if (!na.rm) {
        NAs <- NAs | anyNA(x)
    }
    if (any(NAs)) {
        ans[NAs] <- as(NA, elementType(x))
    }
    ans
})

