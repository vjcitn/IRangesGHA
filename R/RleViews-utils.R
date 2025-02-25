### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "viewApply", "viewMins", "viewMaxs", and "viewSums" generics and
### methods.
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

setMethod("viewMins", "RleViews",
          function(x, na.rm = FALSE)
          .Call2("C_viewMins_RleViews", trim(x), na.rm, PACKAGE="IRangesGHA"))

setMethod("viewMaxs", "RleViews",
          function(x, na.rm = FALSE)
          .Call2("C_viewMaxs_RleViews", trim(x), na.rm, PACKAGE="IRangesGHA"))

setMethod("viewSums", "RleViews",
          function(x, na.rm = FALSE)
          .Call2("C_viewSums_RleViews", trim(x), na.rm, PACKAGE="IRangesGHA"))

setMethod("viewMeans", "RleViews",
          function(x, na.rm = FALSE)
          .Call2("C_viewMeans_RleViews", trim(x), na.rm, PACKAGE="IRangesGHA"))

setMethod("viewWhichMins", "RleViews",
          function(x, na.rm = FALSE)
          .Call2("C_viewWhichMins_RleViews", trim(x), na.rm, PACKAGE="IRangesGHA"))

setMethod("viewWhichMaxs", "RleViews",
          function(x, na.rm = FALSE)
          .Call2("C_viewWhichMaxs_RleViews", trim(x), na.rm, PACKAGE="IRangesGHA"))

setMethod("viewRangeMaxs", "RleViews",
          function(x, na.rm = FALSE) {
              maxs <- viewWhichMaxs(trim(x), na.rm = na.rm)
              if (S4Vectors:::anyMissing(maxs))
                  stop("missing values present, set 'na.rm = TRUE'")
              findRange(maxs, subject(x))
          })

setMethod("viewRangeMins", "RleViews",
          function(x, na.rm = FALSE) {
              mins <- viewWhichMins(trim(x), na.rm = na.rm)
              if (S4Vectors:::anyMissing(mins))
                  stop("missing values present, set 'na.rm = TRUE'")
              findRange(mins, subject(x))
          })
