### =========================================================================
### makeIRangesFromDataFrame()
### -------------------------------------------------------------------------


### NOT exported but used in the GenomicRanges package.
normarg_field <- function(field, what)
{
    if (!is.character(field) || any(is.na(field)))
        stop("'", what, ".field' must be a character vector with no NAs")
    tolower(field)
}

.collect_prefixes <- function(x, suffixes)
{
    stopifnot(is.character(x), is.character(suffixes))
    all_prefixes <- lapply(suffixes,
        function(suffix) {
            ok <- S4Vectors:::has_suffix(x, suffix) & x != suffix
            x2 <- x[ok]
            prefix_nc <- nchar(x2) - nchar(suffix)
            substr(x2, 1L, prefix_nc)
        })
    unique(unlist(all_prefixes, use.names=FALSE))
}

### NOT exported but used in the GenomicRanges package.
get_data_frame_col_as_numeric <- function(df, colidx)
{
    stopifnot(isSingleInteger(colidx))
    col <- df[[colidx]]
    if (is(col, "Rle"))
        col <- S4Vectors:::decodeRle(col)
    if (is.numeric(col))
        return(col)
    if (is.factor(col)) {
        col <- as.character(col)
    } else if (!is.vector(col)) {
        stop(wmsg("the \"", names(df)[[colidx]], "\" column is not ",
                  "an atomic vector, list, factor, or Rle object"))
    }
    ## as.numeric() will generate a warning if the coercion introduces NAs.
    old_warn <- getOption("warn")
    options(warn=2)
    on.exit(options(warn=old_warn))
    as.numeric(col)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_core_IRanges_cols()
###

### NOT exported but used in the GenomicRanges package.
find_start_end_cols <- function(df_colnames, start.field, end.field,
                                startend_prefix="")
{
    colidx1 <- which(df_colnames %in% paste0(startend_prefix, start.field))
    colidx2 <- which(df_colnames %in% paste0(startend_prefix, end.field))
    if (length(colidx1) == 1L && length(colidx2) == 1L)
        return(list(c(start=colidx1, end=colidx2), startend_prefix))
    if (length(colidx1) != 0L || length(colidx2) != 0L ||
        startend_prefix != "")
    {
        stop(wmsg("cannnot determine start/end columns"))
    }
    prefixes1 <- .collect_prefixes(df_colnames, start.field)
    prefixes2 <- .collect_prefixes(df_colnames, end.field)
    if (length(prefixes1) != 1L || length(prefixes2) != 1L ||
        prefixes1 != prefixes2)
    {
        stop(wmsg("cannnot determine start/end columns"))
    }
    startend_prefix <- prefixes1
    find_start_end_cols(df_colnames, start.field, end.field, startend_prefix)
}

### NOT exported but used in the GenomicRanges package.
find_width_col <- function(df_colnames, width.field, startend_prefix)
{
    colidx <- which(df_colnames %in% paste0(startend_prefix, width.field))
    if (length(colidx) == 0L)
        colidx <- which(df_colnames %in% width.field)
    if (length(colidx) == 0L)
        return(NA_integer_)
    if (length(colidx) >= 2L) {
        warning("cannnot determine width column unambiguously")
        return(colidx[[1L]])
    }
    colidx
}

### NOT exported but used in the unit tests.
find_core_IRanges_cols <-
    function(df_colnames, start.field="start", end.field=c("end", "stop"))
{
    ## The heuristic we use to find the core GRanges columns is
    ## case insensitive.
    df_colnames0 <- tolower(df_colnames)
    start.field0 <- normarg_field(start.field, "start")
    end.field0 <- normarg_field(end.field, "end")

    start_end_cols <- find_start_end_cols(df_colnames0,
                                          start.field0, end.field0)
    startend_prefix <- start_end_cols[[2L]]
    ## Name of "width" field is not under user control for now (until we
    ## need that).
    width_col <- find_width_col(df_colnames0, "width", startend_prefix)

    c(start_end_cols[[1L]], width=width_col)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeIRangesFromDataFrame()
###

### NOT exported but used in the GenomicRanges package.
drop_rows_with_na_start_end <- function(df, corecol_map, na.rm)
{
    start_col <- get_data_frame_col_as_numeric(df, corecol_map[["start"]])
    end_col <- get_data_frame_col_as_numeric(df, corecol_map[["end"]])
    is_na <- is.na(start_col) | is.na(end_col)
    if (!any(is_na))
        return(df)
    if (na.rm) {
        keep_idx <- which(!is_na)
        ans <- S4Vectors:::extract_data_frame_rows(df, keep_idx)
        df_rownames <- rownames(df)
        if (!identical(df_rownames, as.character(seq_len(nrow(df)))))
            rownames(ans) <- df_rownames[keep_idx]
        return(ans)
    }
    start_colname <- names(df)[[corecol_map[["start"]]]]
    end_colname <- names(df)[[corecol_map[["end"]]]]
    where <- c("\"", start_colname, "\" and/or \"", end_colname, "\" columns")
    stop(wmsg(
        "The ", where, " contain NAs. Use 'na.rm=TRUE' to ignore ",
        "input rows with NAs in the ", where, "."
    ))
}

### 'df' must be a data.frame (or tibble) or DataFrame object.
makeIRangesFromDataFrame <- function(df,
                                     keep.extra.columns=FALSE,
                                     start.field="start",
                                     end.field=c("end", "stop"),
                                     starts.in.df.are.0based=FALSE,
                                     na.rm=FALSE)
{
    ## Check args.
    if (is.character(df))  # for people that provide the path to a file
        stop("'df' must be a data.frame or DataFrame object")
    if (!(is.data.frame(df) || is(df, "DataFrame")))
        df <- as.data.frame(df)
    if (!isTRUEorFALSE(keep.extra.columns))
        stop("'keep.extra.columns' must be TRUE or FALSE")
    if (!isTRUEorFALSE(starts.in.df.are.0based))
        stop("'starts.in.df.are.0based' must be TRUE or FALSE")
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    corecol_map <- find_core_IRanges_cols(names(df),
                                          start.field=start.field,
                                          end.field=end.field)
    df <- drop_rows_with_na_start_end(df, corecol_map, na.rm)

    ## Prepare 'ans_start' and 'ans_end'.
    ans_start <- get_data_frame_col_as_numeric(df, corecol_map[["start"]])
    ans_end <- get_data_frame_col_as_numeric(df, corecol_map[["end"]])
    if (starts.in.df.are.0based)
        ans_start <- ans_start + 1L

    ans <- new_IRanges(ans_start, ans_end)

    ## Add names.
    ans_names <- rownames(df)
    if (!identical(ans_names, as.character(seq_len(nrow(df)))))
        names(ans) <- ans_names

    ## Add metadata columns.
    if (keep.extra.columns) {
        drop_colidx <- c(corecol_map[["start"]], corecol_map[["end"]])
        if (!is.na(corecol_map[["width"]]))
            drop_colidx <- c(drop_colidx, corecol_map[["width"]])
        ans_mcols <- df[-drop_colidx]
	if (length(ans_mcols) != 0L)
            mcols(ans) <- ans_mcols
    }

    ans
}

setAs("data.frame", "IRanges",
    function(from) makeIRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

setAs("DataFrame", "IRanges",
    function(from) makeIRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

