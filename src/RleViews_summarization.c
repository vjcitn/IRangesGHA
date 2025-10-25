/****************************************************************************
 *             View summarization methods for RleViews objects              *
 ****************************************************************************/
#include "IRanges.h"
#include "S4Vectors_interface.h"

#include <limits.h>  /* for INT_MAX */


/****************************************************************************
 * Low-level C_summarize_RleViews() helpers
 */

static void viewMins_int_RleViews(
		const int *run_vals,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		int nranges, int narm, int *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		int min;  /* uninitialized */
		int min_is_initialized = 0;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			int v = run_vals[k];
			if (v == NA_INTEGER) {
				if (narm)
					continue;
				min = NA_INTEGER;
				min_is_initialized = 1;
				break;
			}
			if (!min_is_initialized || v < min) {
				min = v;
				min_is_initialized = 1;
			}
		}
		out[i] = min_is_initialized ? min : NA_INTEGER;
	}
	return;
}

static void viewMaxs_int_RleViews(const int *run_vals,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		int nranges, int narm, int *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		int max;  /* uninitialized */
		int max_is_initialized = 0;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			int v = run_vals[k];
			if (v == NA_INTEGER) {
				if (narm)
					continue;
				max = NA_INTEGER;
				max_is_initialized = 1;
				break;
			}
			if (!max_is_initialized || v > max) {
				max = v;
				max_is_initialized = 1;
			}
		}
		out[i] = max_is_initialized ? max : NA_INTEGER;
	}
	return;
}

static void viewMins_double_RleViews(
		const double *run_vals,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		int nranges, int narm, double *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		double min = R_PosInf;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			double v = run_vals[k];
			if (ISNAN(v)) {
				if (narm)
					continue;
				min = NA_REAL;
				break;
			}
			if (v < min)
				min = v;
		}
		out[i] = min;
	}
	return;
}

static void viewMaxs_double_RleViews(const double *run_vals,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		int nranges, int narm, double *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		double max = R_NegInf;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			double v = run_vals[k];
			if (ISNAN(v)) {
				if (narm)
					continue;
				max = NA_REAL;
				break;
			}
			if (v > max)
				max = v;
		}
		out[i] = max;
	}
	return;
}

static void viewSums_int_RleViews(
		const int *run_vals, const int *run_lens,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		const int *mapped_range_Rtrim,
		int nranges, int narm, double *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		double sum = 0.0;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			int v = run_vals[k];
			if (v == NA_INTEGER) {
				if (narm)
					continue;
				sum = NA_REAL;
				break;
			}
			int n = run_lens[k];
			if (k == k1)
				n -= mapped_range_Ltrim[i];
			if (k == k2)
				n -= mapped_range_Rtrim[i];
			sum += (double) v * n;
		}
		out[i] = sum;
	}
	return;
}

static void viewSums_double_RleViews(
		const double *run_vals, const int *run_lens,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		const int *mapped_range_Rtrim,
		int nranges, int narm, double *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		double sum = 0.0;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			double v = run_vals[k];
			if (ISNAN(v)) {
				if (narm)
					continue;
				sum = NA_REAL;
				break;
			}
			int n = run_lens[k];
			if (k == k1)
				n -= mapped_range_Ltrim[i];
			if (k == k2)
				n -= mapped_range_Rtrim[i];
			sum += v * n;
		}
		out[i] = sum;
	}
	return;
}

static void viewMeans_int_RleViews(
		const int *run_vals, const int *run_lens,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		const int *mapped_range_Rtrim,
		int nranges, int narm, double *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		long long int nval = 0;
		double sum = 0.0;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			int v = run_vals[k];
			if (v == NA_INTEGER) {
				if (narm)
					continue;
				nval = 1;  /* any nonzero val works */
				sum = NA_REAL;
				break;
			}
			int n = run_lens[k];
			if (k == k1)
				n -= mapped_range_Ltrim[i];
			if (k == k2)
				n -= mapped_range_Rtrim[i];
			nval += n;
			sum += (double) v * n;
		}
		out[i] = nval == 0 ? R_NaN : (sum / nval);
	}
	return;
}

static void viewMeans_double_RleViews(
		const double *run_vals, const int *run_lens,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		const int *mapped_range_Rtrim,
		int nranges, int narm, double *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		long long int nval = 0;
		double sum = 0.0;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			double v = run_vals[k];
			if (ISNAN(v)) {
				if (narm)
					continue;
				nval = 1;  /* any nonzero val works */
				sum = NA_REAL;
				break;
			}
			int n = run_lens[k];
			if (k == k1)
				n -= mapped_range_Ltrim[i];
			if (k == k2)
				n -= mapped_range_Rtrim[i];
			nval += n;
			sum += v * n;
		}
		out[i] = nval == 0 ? R_NaN : (sum / nval);
	}
	return;
}

static void viewWhichMins_int_RleViews(
		const int *run_vals, const int *run_starts,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		int nranges, int *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		int min;  /* uninitialized */
		int min_is_initialized = 0;
		int which_min = NA_INTEGER;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			int v = run_vals[k];
			if (v == NA_INTEGER)
				continue;
			if (!min_is_initialized || v < min) {
				min = v;
				min_is_initialized = 1;
				which_min = run_starts[k];
				if (k == k1)
					which_min += mapped_range_Ltrim[i];
			}
		}
		out[i] = which_min;
	}
	return;
}

static void viewWhichMaxs_int_RleViews(
		const int *run_vals, const int *run_starts,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		int nranges, int *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		int max;  /* uninitialized */
		int max_is_initialized = 0;
		int which_max = NA_INTEGER;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			int v = run_vals[k];
			if (v == NA_INTEGER)
				continue;
			if (!max_is_initialized || v > max) {
				max = v;
				max_is_initialized = 1;
				which_max = run_starts[k];
				if (k == k1)
					which_max += mapped_range_Ltrim[i];
			}
		}
		out[i] = which_max;
	}
	return;
}

static void viewWhichMins_double_RleViews(
		const double *run_vals, const int *run_starts,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		int nranges, int *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		double min = R_PosInf;
		int which_min = NA_INTEGER;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			double v = run_vals[k];
			if (ISNAN(v))
				continue;
			if (v < min) {
				min = v;
				which_min = run_starts[k];
				if (k == k1)
					which_min += mapped_range_Ltrim[i];
			}
		}
		out[i] = which_min;
	}
	return;
}

static void viewWhichMaxs_double_RleViews(
		const double *run_vals, const int *run_starts,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		int nranges, int *out)
{
	for (int i = 0; i < nranges; i++) {
		if (i % 100000 == 99999)
			R_CheckUserInterrupt();
		double max = R_NegInf;
		int which_max = NA_INTEGER;
		int k1 = mapped_range_offset[i];
		int k2 = k1 + mapped_range_span[i] - 1;
		for (int k = k1; k <= k2; k++) {
			double v = run_vals[k];
			if (ISNAN(v))
				continue;
			if (v > max) {
				max = v;
				which_max = run_starts[k];
				if (k == k1)
					which_max += mapped_range_Ltrim[i];
			}
		}
		out[i] = which_max;
	}
	return;
}


/****************************************************************************
 * High-level C_summarize_RleViews() helpers
 */

static SEXP viewMins_RleViews(SEXP run_vals,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		int nranges, int narm)
{
	SEXP ans;
	SEXPTYPE inputRtype = TYPEOF(run_vals);
	switch (inputRtype) {
	    case INTSXP: case LGLSXP:
		ans = PROTECT(NEW_INTEGER(nranges));
		viewMins_int_RleViews(INTEGER(run_vals),
				mapped_range_offset,
				mapped_range_span,
				nranges, narm, INTEGER(ans));
	    break;
	    case REALSXP:
		ans = PROTECT(NEW_NUMERIC(nranges));
		viewMins_double_RleViews(REAL(run_vals),
				mapped_range_offset,
				mapped_range_span,
				nranges, narm, REAL(ans));
	    break;
	    default:
		error("viewMins_RleViews() is not implemented for Rle "
		      "subjects with \"%s\" values", type2char(inputRtype));
	}
	UNPROTECT(1);
	return ans;
}

static SEXP viewMaxs_RleViews(SEXP run_vals,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		int nranges, int narm)
{
	SEXP ans;
	SEXPTYPE inputRtype = TYPEOF(run_vals);
	switch (inputRtype) {
	    case INTSXP: case LGLSXP:
		ans = PROTECT(NEW_INTEGER(nranges));
		viewMaxs_int_RleViews(INTEGER(run_vals),
				mapped_range_offset,
				mapped_range_span,
				nranges, narm, INTEGER(ans));
	    break;
	    case REALSXP:
		ans = PROTECT(NEW_NUMERIC(nranges));
		viewMaxs_double_RleViews(REAL(run_vals),
				mapped_range_offset,
				mapped_range_span,
				nranges, narm, REAL(ans));
	    break;
	    default:
		error("viewMaxs_RleViews() is not implemented for Rle "
		      "subjects with \"%s\" values", type2char(inputRtype));
	}
	UNPROTECT(1);
	return ans;
}

static SEXP viewSums_RleViews(SEXP run_vals, const int *run_lens,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		const int *mapped_range_Rtrim,
		int nranges, int narm)
{
	SEXP ans;
	SEXPTYPE inputRtype = TYPEOF(run_vals);
	switch (inputRtype) {
	    case INTSXP: case LGLSXP:
		ans = PROTECT(NEW_NUMERIC(nranges));
		viewSums_int_RleViews(INTEGER(run_vals), run_lens,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				mapped_range_Rtrim,
				nranges, narm, REAL(ans));
	    break;
	    case REALSXP:
		ans = PROTECT(NEW_NUMERIC(nranges));
		viewSums_double_RleViews(REAL(run_vals), run_lens,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				mapped_range_Rtrim,
				nranges, narm, REAL(ans));
	    break;
	    default:
		error("viewSums_RleViews() is not implemented for Rle "
		      "subjects with \"%s\" values", type2char(inputRtype));
	}
	UNPROTECT(1);
	return ans;
}

static SEXP viewMeans_RleViews(SEXP run_vals, const int *run_lens,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		const int *mapped_range_Rtrim,
		int nranges, int narm)
{
	SEXP ans;
	SEXPTYPE inputRtype = TYPEOF(run_vals);
	switch (inputRtype) {
	    case INTSXP: case LGLSXP:
		ans = PROTECT(NEW_NUMERIC(nranges));
		viewMeans_int_RleViews(INTEGER(run_vals), run_lens,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				mapped_range_Rtrim,
				nranges, narm, REAL(ans));
	    break;
	    case REALSXP:
		ans = PROTECT(NEW_NUMERIC(nranges));
		viewMeans_double_RleViews(REAL(run_vals), run_lens,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				mapped_range_Rtrim,
				nranges, narm, REAL(ans));
	    break;
	    default:
		error("viewMeans_RleViews() is not implemented for Rle "
		      "subjects with \"%s\" values", type2char(inputRtype));
	}
	UNPROTECT(1);
	return ans;
}

static SEXP viewWhichMins_RleViews(SEXP run_vals, const int *run_starts,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		int nranges)
{
	SEXP ans;
	SEXPTYPE inputRtype = TYPEOF(run_vals);
	switch (inputRtype) {
	    case INTSXP: case LGLSXP:
		ans = PROTECT(NEW_INTEGER(nranges));
		viewWhichMins_int_RleViews(INTEGER(run_vals), run_starts,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				nranges, INTEGER(ans));
	    break;
	    case REALSXP:
		ans = PROTECT(NEW_INTEGER(nranges));
		viewWhichMins_double_RleViews(REAL(run_vals), run_starts,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				nranges, INTEGER(ans));
	    break;
	    default:
		error("viewWhichMins_RleViews() is not implemented for Rle "
		      "subjects with \"%s\" values", type2char(inputRtype));
	}
	UNPROTECT(1);
	return ans;
}

static SEXP viewWhichMaxs_RleViews(SEXP run_vals, const int *run_starts,
		const int *mapped_range_offset,
		const int *mapped_range_span,
		const int *mapped_range_Ltrim,
		int nranges)
{
	SEXP ans;
	SEXPTYPE inputRtype = TYPEOF(run_vals);
	switch (inputRtype) {
	    case INTSXP: case LGLSXP:
		ans = PROTECT(NEW_INTEGER(nranges));
		viewWhichMaxs_int_RleViews(INTEGER(run_vals), run_starts,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				nranges, INTEGER(ans));
	    break;
	    case REALSXP:
		ans = PROTECT(NEW_INTEGER(nranges));
		viewWhichMaxs_double_RleViews(REAL(run_vals), run_starts,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				nranges, INTEGER(ans));
	    break;
	    default:
		error("viewWhichMaxs_RleViews() is not implemented for Rle "
		      "subjects with \"%s\" values", type2char(inputRtype));
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_summarize_RleViews()
 */

static void compute_run_starts(const int *run_lens, int nrun, int *run_starts)
{
	unsigned int start = 1;
	for (int k = 0; k < nrun; k++) {
		run_starts[k] = start;
		start += run_lens[k];
		if (start > (unsigned int) INT_MAX + 1)
			error("operation not supported when "
			      "the subject is a **long** Rle");
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_summarize_RleViews(SEXP op, SEXP x_subject, SEXP x_ranges, SEXP na_rm)
{
	if (!(IS_CHARACTER(op) && LENGTH(op) == 1))
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	const char *op_string = CHAR(op);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	int narm = LOGICAL(na_rm)[0];
	if (narm == NA_INTEGER)
		error("'na.rm' cannot be NA");

	SEXP run_vals = GET_SLOT(x_subject, install("values"));
	SEXP run_lens = GET_SLOT(x_subject, install("lengths"));
	int nrun = LENGTH(run_lens);

	SEXP x_start = _get_IRanges_start(x_ranges);
	SEXP x_width = _get_IRanges_width(x_ranges);
	const int *start_p, *width_p;
	int nranges = check_integer_pairs(x_start, x_width,
					  &start_p, &width_p,
					  "start", "width");

	int *mapped_range_offset = (int *) R_alloc(nranges, sizeof(int));
	int *mapped_range_span   = (int *) R_alloc(nranges, sizeof(int));
	int *mapped_range_Ltrim  = (int *) R_alloc(nranges, sizeof(int));
	int *mapped_range_Rtrim  = (int *) R_alloc(nranges, sizeof(int));
	const char *errmsg = ranges_mapper(INTEGER(run_lens), nrun,
					start_p, width_p, nranges,
					mapped_range_offset,
					mapped_range_span,
					mapped_range_Ltrim,
					mapped_range_Rtrim, 0);
	if (errmsg != NULL)
		error("%s", errmsg);
	SEXP ans = R_NilValue;
	if (strcmp(op_string, "min") == 0) {
		ans = viewMins_RleViews(run_vals,
				mapped_range_offset,
				mapped_range_span,
				nranges, narm);
	} else if (strcmp(op_string, "max") == 0) {
		ans = viewMaxs_RleViews(run_vals,
				mapped_range_offset,
				mapped_range_span,
				nranges, narm);
	} else if (strcmp(op_string, "sum") == 0) {
		ans = viewSums_RleViews(run_vals, INTEGER(run_lens),
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				mapped_range_Rtrim,
				nranges, narm);
	} else if (strcmp(op_string, "mean") == 0) {
		ans = viewMeans_RleViews(run_vals, INTEGER(run_lens),
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				mapped_range_Rtrim,
				nranges, narm);
	} else if (strcmp(op_string, "which.min") == 0) {
		int *run_starts = (int *) R_alloc(nrun, sizeof(int));
		compute_run_starts(INTEGER(run_lens), nrun, run_starts);
		ans = viewWhichMins_RleViews(run_vals, run_starts,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				nranges);
	} else if (strcmp(op_string, "which.max") == 0) {
		int *run_starts = (int *) R_alloc(nrun, sizeof(int));
		compute_run_starts(INTEGER(run_lens), nrun, run_starts);
		ans = viewWhichMaxs_RleViews(run_vals, run_starts,
				mapped_range_offset,
				mapped_range_span,
				mapped_range_Ltrim,
				nranges);
	}
	if (ans == R_NilValue)
		error("'op == \"%s\"' not ready yet", op_string);
	PROTECT(ans);
	SEXP names = PROTECT(duplicate(_get_IRanges_names(x_ranges)));
	SET_NAMES(ans, names);
	UNPROTECT(2);
	return ans;
}

