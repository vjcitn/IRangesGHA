#include "IRanges.h"


static int debug = 0;

SEXP debug_copy_byteblocks()
{
#ifdef DEBUG_IRANGES
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}


/*
 * Memory comparison
 */
int _compare_byteblocks(const char *a, int ia, const char *b, int ib,
		int n, size_t blocksize)
{
#ifdef DEBUG_IRANGES
	if (debug) {
		Rprintf("[DEBUG] _compare_byteblocks(): ");
		Rprintf("a=%p ia=%d b=%p ib=%d n=%d blocksize=%d\n",
			a, ia, b, ib, n, blocksize);
	}
#endif
	a += ia * blocksize;
	b += ib * blocksize;
	/* memcmp() doesn't try to be smart by checking if a == b */
	return a == b ? 0 : memcmp(a, b, n * blocksize);
}


/****************************************************************************
 All the functions below are performing cyclic copy i.e. copy with recycling
 either at the destination ('dest') or at the source ('src'). In this file,
 "Ocopy" is an abbreviation for "cyclic copy".
 ****************************************************************************/

/*
 * Performs (in short):
 *   dest[(i-i1) % dest_nblocks] <- src[i] for i1 <= i <= i2
 * Details:
 *   - Reads the linear subset of blocks from 'src' defined by 'i1', 'i2'.
 *   - Writing is recycled in 'dest': it starts at its first block
 *     and comes back to it after it reaches its last block.
 *   - Don't do anything if i1 > i2.
 */
void _Ocopy_byteblocks_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nblocks,
		const char *src, size_t src_nblocks, size_t blocksize)
{
	const char *b;
	int i2next, i1max, q;
	size_t dest_size;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= src_nblocks)
		error("subscript out of bounds");
	if (dest_nblocks == 0)
		error("no destination to copy to");
	i2next = i2 + 1;
	i1max = i2next - dest_nblocks;
	b = src + i1 * blocksize;
	dest_size = dest_nblocks * blocksize;
	while (i1 <= i1max) {
		memcpy(dest, b, dest_size);
		b += dest_size;
		i1 += dest_nblocks;
	}
	q = i2next - i1;
	if (q > 0) {
		/* Safe because q is always < dest_nblocks */
		memcpy(dest, b, q * blocksize);
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return;
}

/*
 * Performs (in short):
 *   dest[k % dest_nblocks] <- src[subset[k] - 1] for 0 <= k <= n
 * Details:
 *   - Reads the blocks from 'src' that have the 1-based offsets passed
 *     in 'subset'.
 *   - Writing is recycled in 'dest': it starts at its first block
 *     and comes back to it after it reaches its last block.
 */
void _Ocopy_byteblocks_from_subset(const int *subset, int n,
		char *dest, size_t dest_nblocks,
		const char *src, size_t src_nblocks, size_t blocksize)
{
	char *a;
	const char *b;
	int i, j, k, z;

	if (dest_nblocks == 0 && n != 0)
		error("no destination to copy to");
	a = dest;
	for (i = k = 0; k < n; i++, k++) {
		j = subset[k] - 1;
		if (j < 0 || j >= src_nblocks)
			error("subscript out of bounds");
		if (i >= dest_nblocks) {
			i = 0; /* recycle */
			a = dest;
		}
		b = src + j * blocksize;
		for (z = 0; z < blocksize; z++) {
			*(a++) = *(b++);
		}
	}
	if (i != dest_nblocks)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[i] <- src[(i-i1) % src_nblocks] for i1 <= i <= i2
 * Details:
 *   - Writes to the linear subset of blocks in 'dest' defined by 'i1', 'i2'.
 *   - Reading is recycled in 'src': it starts at its first block
 *     and comes back to it after it reaches its last block.
 *   - Don't do anything if i1 > i2.
 */
void _Ocopy_byteblocks_to_i1i2(int i1, int i2,
		char *dest, size_t dest_nblocks,
		const char *src, size_t src_nblocks, size_t blocksize)
{
	char *a;
	int i2next, i1max, q;
	size_t src_size;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= dest_nblocks)
		error("subscript out of bounds");
	if (src_nblocks == 0)
		error("no value provided");
	i2next = i2 + 1;
	i1max = i2next - src_nblocks;
	a = dest + i1 * blocksize;
	src_size = src_nblocks * blocksize;
	while (i1 <= i1max) {
		memcpy(a, src, src_size);
		a += src_size;
		i1 += src_nblocks;
	}
	q = i2next - i1;
	if (q > 0) {
		/* Safe because q is always < src_nblocks */
		memcpy(a, src, q * blocksize);
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return;
}

/*
 * Performs (in short):
 *   dest[subset[k] - 1] <- src[k % src_nblocks] for 0 <= k <= n
 * Details:
 *   - Writes the blocks in 'dest' that have the 1-based offsets passed
 *     in 'subset'.
 *   - Reading is recycled in 'src': it starts at its first block
 *     and comes back to it after it reaches its last block.
 */
void _Ocopy_byteblocks_to_subset(const int *subset, int n,
		char *dest, size_t dest_nblocks,
		const char *src, size_t src_nblocks, size_t blocksize)
{
	char *a;
	const char *b;
	int i, j, k, z;

	if (src_nblocks == 0 && n != 0)
		error("no value provided");
	b = src;
	for (j = k = 0; k < n; j++, k++) {
		i = subset[k] - 1;
		if (i < 0 || i >= dest_nblocks)
			error("subscript out of bounds");
		if (j >= src_nblocks) {
			j = 0; /* recycle */
			b = src;
		}
		a = dest + i * blocksize;
		for (z = 0; z < blocksize; z++) {
			*(a++) = *(b++);
		}
	}
	if (j != src_nblocks)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[(i-i1) % dest_nbytes] <- tr(src[i]) for i1 <= i <= i2
 * Note: tr() stands for translation.
 * Details:
 *   - Reads the linear subset of bytes from 'src' defined by 'i1', 'i2'.
 *   - Writing is recycled in 'dest': it starts at its first byte
 *     and comes back to it after it reaches its last byte.
 *   - Don't do anything if i1 > i2.
 */
void _Ocopy_bytes_from_i1i2_with_lkup(int i1, int i2,
		char *dest, int dest_nbytes,
		const char *src, int src_nbytes,
		const int *lkup, int lkup_length)
{
	const char *b;
	char src_val;
	int i, j, lkup_key, lkup_val;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= src_nbytes)
		error("subscript out of bounds");
	if (dest_nbytes == 0)
		error("no destination to copy to");
	b = src + i1;
	for (i = i1, j = 0; i <= i2; i++, j++) {
		if (j >= dest_nbytes) { /* recycle */
			j = 0;
		}
		src_val = *(b++);
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = (char) lkup_val;
	}
	if (j < dest_nbytes)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[k % dest_nbytes] <- tr(src[subset[k] - 1]) for 0 <= k <= n
 * Note: tr() stands for translation.
 * Details:
 *   - Reads the bytes from 'src' that have the 1-based offsets passed
 *     in 'subset'.
 *   - Writing is recycled in 'dest': it starts at its first byte
 *     and comes back to it after it reaches its last byte.
 */
void _Ocopy_bytes_from_subset_with_lkup(const int *subset, int n,
		char *dest, int dest_nbytes,
		const char *src, int src_nbytes,
		const int *lkup, int lkup_length)
{
	char src_val;
	int i, j, k, lkup_key, lkup_val;

	if (dest_nbytes == 0 && n != 0)
		error("no destination to copy to");
	for (k = j = 0; k < n; k++, j++) {
		if (j >= dest_nbytes) { /* recycle */
			j = 0;
		}
		i = subset[k] - 1;
		if (i < 0 || i >= src_nbytes)
			error("subscript out of bounds");
		src_val = src[i];
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = (char) lkup_val;
	}
	if (j < dest_nbytes)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[i] <- tr(src[(i-i1) % src_nbytes]) for i1 <= i <= i2
 * Note: tr() stands for translation.
 * Details:
 *   - Writes to the linear subset of bytes in 'dest' defined by 'i1', 'i2'.
 *   - Reading is recycled in 'src': it starts at its first byte
 *     and comes back to it after it reaches its last byte.
 *   - Don't do anything if i1 > i2.
 */
void _Ocopy_bytes_to_i1i2_with_lkup(int i1, int i2,
		char *dest, int dest_nbytes,
		const char *src, int src_nbytes,
		const int *lkup, int lkup_length)
{
	char *a, src_val;
	int i, j, lkup_key, lkup_val;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= dest_nbytes)
		error("subscript out of bounds");
	if (src_nbytes == 0)
		error("no value provided");
	a = dest + i1;
	for (i = i1, j = 0; i <= i2; i++, j++) {
		if (j >= src_nbytes) { /* recycle */
			j = 0;
		}
		src_val = src[j];
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		*(a++) = (char) lkup_val;
	}
	if (j < src_nbytes)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[subset[k] - 1] <- tr(src[k % src_nbytes]) for 0 <= k <= n
 * Note: tr() stands for translation.
 * Details:
 *   - Writes the bytes in 'dest' that have the 1-based offsets passed
 *     in 'subset'.
 *   - Reading is recycled in 'src': it starts at its first byte
 *     and comes back to it after it reaches its last byte.
 */
void _Ocopy_bytes_to_subset_with_lkup(const int *subset, int n,
		char *dest, int dest_nbytes,
		const char *src, int src_nbytes,
		const int *lkup, int lkup_length)
{
	char src_val;
	int i, j, k, lkup_key, lkup_val;

	if (src_nbytes == 0 && n != 0)
		error("no value provided");
	for (k = j = 0; k < n; k++, j++) {
		if (j >= src_nbytes) { /* recycle */
			j = 0;
		}
		i = subset[k] - 1;
		if (i < 0 || i >= dest_nbytes)
			error("subscript out of bounds");
		src_val = src[j];
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[i] = (char) lkup_val;
	}
	if (j < src_nbytes)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[(dest_nblocks-1-(i-i1)) % dest_nblocks] <- src[i] for i1 <= i <= i2
 * Note: the order of the blocks is reversed during the copy.
 * Details:
 *   - Reads the linear subset of blocks from 'src' defined by 'i1', 'i2'.
 *   - Writing is recycled in 'dest': it starts at its last block
 *     and comes back to it after it reaches its first block.
 *   - Don't do anything if i1 > i2.
 */
void _Orevcopy_byteblocks_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nblocks,
		const char *src, size_t src_nblocks, size_t blocksize)
{
	char *a;
	const char *b;
	int i, j, z;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= src_nblocks)
		error("subscript out of bounds");
	if (dest_nblocks == 0)
		error("no destination to copy to");
	b = src + i1 * blocksize;
	for (i = i1, j = dest_nblocks - 1; i <= i2; i++, j--) {
		if (j < 0) { /* recycle */
			j = dest_nblocks - 1;
		}
		a = dest + j * blocksize;
		for (z = 0; z < blocksize; z++) {
			*(a++) = *(b++);
		}
	}
	if (j >= 0)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[(dest_nbytes-1-(i-i1)) % dest_nbytes] <- tr(src[i]) for i1 <= i <= i2
 * Notes: - tr() stands for translation.
 *        - the order of the bytes is reversed during the copy.
 * Details:
 *   - Reads the linear subset of bytes from 'src' defined by 'i1', 'i2'.
 *   - Writing is recycled in 'dest': it starts at its last byte
 *     and comes back to it after it reaches its first byte.
 *   - Don't do anything if i1 > i2.
 */
void _Orevcopy_bytes_from_i1i2_with_lkup(int i1, int i2,
		char *dest, int dest_nbytes,
		const char *src, int src_nbytes,
		const int *lkup, int lkup_length)
{
	const char *b;
	char src_val;
	int i, j, lkup_key, lkup_val;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= src_nbytes)
		error("subscript out of bounds");
	if (dest_nbytes == 0)
		error("no destination to copy to");
	b = src + i1;
	for (i = i1, j = dest_nbytes - 1; i <= i2; i++, j--) {
		if (j < 0) { /* recycle */
			j = dest_nbytes - 1;
		}
		src_val = *(b++);
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = (char) lkup_val;
	}
	if (j >= 0)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

/*
 * Performs (in short):
 *   dest[(i-i1) % dest_nbytes] <- toComplex(src[i]) for i1 <= i <= i2
 * Note: toComplex() stands for conversion to complex values.
 * Details:
 *   - Reads the linear subset of bytes from 'src' defined by 'i1', 'i2'.
 *   - Writing is recycled in 'dest': it starts at its first element
 *     and comes back to it after it reaches its last element.
 *   - Don't do anything if i1 > i2.
 */
void _Ocopy_bytes_from_i1i2_to_complex(int i1, int i2,
		Rcomplex *dest, int dest_nbytes,
		const char *src, int src_nbytes,
		const Rcomplex *lkup, int lkup_length)
{
	const char *b;
	char src_val;
	int i, j, lkup_key;
	Rcomplex lkup_val;

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= src_nbytes)
		error("subscript out of bounds");
	if (dest_nbytes == 0)
		error("no destination to copy to");
	b = src + i1;
	for (i = i1, j = 0; i <= i2; i++, j++) {
		if (j >= dest_nbytes) { /* recycle */
			j = 0;
		}
		src_val = *(b++);
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length
		 || ISNA((lkup_val = lkup[lkup_key]).r)
		 || ISNA(lkup_val.i)) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = lkup_val;
	}
	if (j < dest_nbytes)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

