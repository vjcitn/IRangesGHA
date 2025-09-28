### =========================================================================
### OpenMP thread control
### -------------------------------------------------------------------------
###


.normarg_nthread <- function(nthread)
{
    if (!isSingleNumber(nthread))
        stop(wmsg("'nthread' must be a single number"))
    if (!is.integer(nthread))
        nthread <- as.integer(nthread)
    nthread
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### R wrappers to OpenMP thread control functions
###

### Wrapper to omp_get_num_procs().
### Returns the same as parallel::detectCores() or 0 if OpenMP is not
### available (e.g. on macOS).
get_num_procs <- function()
    .Call2("C_get_num_procs", PACKAGE="IRanges")

### Wrapper to omp_get_max_threads().
### Returns 0 if OpenMP is not available (e.g. on macOS).
### Note that the "initial" value returned by get_max_threads(), that is,
### the value it returns **before** omp_set_num_threads() gets called, is
### controlled by environment variable OMP_NUM_THREADS.
get_max_threads <- function()
    .Call2("C_get_max_threads", PACKAGE="IRanges")

### Wrapper to omp_set_num_threads().
### No-op if OpenMP is not available (e.g. on macOS).
### Returns previous omp_get_max_threads() value.
set_max_threads <- function(nthread)
{
    nthread <- .normarg_nthread(nthread)
    .Call2("C_set_max_threads", nthread, PACKAGE="IRanges")
}

