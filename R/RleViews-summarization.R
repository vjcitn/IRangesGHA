### =========================================================================
### View summarization methods for RleViews objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### viewApply() method for RleViews objects
###

setMethod("viewApply", "RleViews",
          function(X, FUN, ..., simplify = TRUE) {
              X <- trim(X)
              ans <-
                aggregate(subject(X), start = structure(start(X), names = names(X)),
                          end = end(X), FUN = FUN, ..., simplify = simplify)
              if (!simplify) {
                  ans <- S4Vectors:::new_SimpleList_from_list("SimpleList",
                                                  ans,
                                                  metadata=metadata(X),
                                                  mcols=mcols(X))
              }
              ans
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### viewMins(), viewMaxs(), viewSums(), viewMeans(), viewWhichMins(), and
### viewWhichMaxs() methods for RleViews objects
###

### Coerce 'x' to integer if all its values are within the integer range.
.double_to_integer_if_in_range <- function(x)
{
    stopifnot(is.double(x))
    x_min <- suppressWarnings(min(x, na.rm=TRUE))
    x_max <- suppressWarnings(max(x, na.rm=TRUE))
    int_max_plus_one <- .Machine$integer.max + 1
    if (-int_max_plus_one < x_min && x_max < int_max_plus_one)
        x <- setNames(as.integer(x), names(x))
    x
}

.C_summarize_RleViews <- function(op, x_subject, x_ranges, na.rm)
{
    C_ans <- .Call2("C_summarize_RleViews", op, x_subject, x_ranges, na.rm,
                                            PACKAGE="IRanges")
    if (op == "sum") {
        run_vals <- runValue(x_subject)
        if (is.integer(run_vals) || is.logical(run_vals))
            C_ans <- .double_to_integer_if_in_range(C_ans)
    }
    C_ans
}

.summarize_RleViews <- function(op, x, na.rm=FALSE)
{
    stopifnot(isSingleString(op), is(x, "RleViews"), isTRUEorFALSE(na.rm))
    x <- trim(x)
    x_ranges <- ranges(x)
    x_subject <- subject(x)
    run_vals <- runValue(x_subject)
    if (!is.complex(run_vals))
        return(.C_summarize_RleViews(op, x_subject, x_ranges, na.rm))
    if (op %in% c("min", "max", "which.min", "which.max"))
        stop(wmsg("operation not supported when the subject ",
                  "is an Rle that contains \"complex\" values"))
    ans_names <- names(x_ranges)
    x_ranges <- unname(x_ranges)
    ans_r <- .C_summarize_RleViews(op, Re(x_subject), x_ranges, na.rm)
    ans_i <- .C_summarize_RleViews(op, Im(x_subject), x_ranges, na.rm)
    setNames(complex(real=ans_r, imaginary=ans_i), ans_names)
}

setMethod("viewMins", "RleViews",
    function(x, na.rm=FALSE) .summarize_RleViews("min", x, na.rm=na.rm)
)

setMethod("viewMaxs", "RleViews",
    function(x, na.rm=FALSE) .summarize_RleViews("max", x, na.rm=na.rm)
)

setMethod("viewSums", "RleViews",
    function(x, na.rm=FALSE) .summarize_RleViews("sum", x, na.rm=na.rm)
)

setMethod("viewMeans", "RleViews",
    function(x, na.rm=FALSE) .summarize_RleViews("mean", x, na.rm=na.rm)
)

### Even though base::which(), base::which.min(), and base::which.max() always
### ignore NAs and don't have an 'na.rm' argument to control that, somehow
### someone felt it could be a good idea to give the viewWhichMins() and
### viewWhichMaxs() generics an 'na.rm' argument, not sure why. Plus, the
### old implementations (from IRanges < 2.41.3) of the viewWhichMins() and
### viewWhichMaxs() methods for RleViews objects (which were based on old .Call
### entry points C_viewWhichMins_RleViews and C_viewWhichMaxs_RleViews) were
### apparently trying to support the 'na.rm' argument but that support was
### broken.
setMethod("viewWhichMins", "RleViews",
    function(x, na.rm=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning(wmsg("the viewWhichMins() method for RleViews objects ",
                         "always ignores NAs (the 'na.rm' argument has ",
                         "no effect)"))
        .summarize_RleViews("which.min", x)
    }
)

setMethod("viewWhichMaxs", "RleViews",
    function(x, na.rm=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning(wmsg("the viewWhichMaxs() method for RleViews objects ",
                         "always ignores NAs (the 'na.rm' argument has ",
                         "no effect)"))
        .summarize_RleViews("which.max", x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### viewRangeMins() and viewRangeMaxs() methods for RleViews objects
###

setMethod("viewRangeMins", "RleViews",
    function(x, na.rm=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning(wmsg("the viewRangeMins() method for RleViews objects ",
                         "always ignores NAs (the 'na.rm' argument has ",
                         "no effect)"))
        mins <- viewWhichMins(x)
        pintersect(findRange(mins, subject(x)), trim(x))
    }
)

setMethod("viewRangeMaxs", "RleViews",
    function(x, na.rm=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning(wmsg("the viewRangeMaxs() method for RleViews objects ",
                         "always ignores NAs (the 'na.rm' argument has ",
                         "no effect)"))
        maxs <- viewWhichMaxs(x)
        pintersect(findRange(maxs, subject(x)), trim(x))
    }
)

