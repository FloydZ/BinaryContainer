#ifndef BINARY_ROW_H
#define BINARY_ROW_H

#include <cstdint>      // needed for 'uint64_t'
#include <cstddef>      // needed for 'size_t'
#include <array>
#include <vector>

#include <iostream>

#include "random.h"

#if defined(USE_AVX) || defined(USE_AVX2)
#include <immintrin.h>
#include <emmintrin.h>
#endif

// performance helpers
#if defined(USE_BRANCH_PREDICTION) && (!defined(DEBUG))
#   define likely(x)       __builtin_expect(!!(x), 1)
#   define unlikely(x)     __builtin_expect(!!(x), 0)
#else
#   define likely(x)       x
#   define unlikely(x)     x
#endif

#ifdef USE_PREFETCH
/*
 * The value of addr is the address of the memory to prefetch. There are two optional arguments, rw and locality.
 * The value of rw is a compile-time constant one or zero; one means that the prefetch is preparing for a write to the
 * memory address and zero the default, means that the prefetch is preparing for a read. The value locality must be a
 * compile-time constant integer between zero and three. A value of zero means that the data has no temporal locality,
 * so it need not be left in the cache after the access. A value of three means that the data has a high degree of
 * temporal locality and should be left in all levels of cache possible. Values of one and two mean, respectively,
 * a low or moderate degree of temporal locality. The default is three.
 */
#   define prefetch(m, x, y) __builtin_prefetch(m, x, y)
#   define cryptanalysislib_prefetch(address)  __builtin_prefetch((const void *)(address), 0, 0)
#   define cryptanalysislib_prefetchw(address) __builtin_prefetch((const void *)(address), 1, 0)
constexpr ptrdiff_t prefetch_distance = 32*4;     //Prefetch amount in bytes
#else
#   define prefetch(m, x, y)
#   define cryptanalysislib_prefetch(address)
#   define cryptanalysislib_prefetchw(address)
constexpr ptrdiff_t prefetch_distance = 0;
#endif

#ifdef USE_LOOP_UNROLL
#define LOOP_UNROLL()                      \
    _Pragma(STRINGIFY(clang loop unroll(full)))
#else
#define LOOP_UNROLL()
#endif

#ifdef DEBUG
#include <cassert>
#   define ASSERT(x) assert(x)
#else
#   define ASSERT(x)
#endif

#ifdef DEBUG
#   define INLINE
#else
#   define INLINE inline
#endif

// helper macros
#define BINARYCONTAINER_COMPARE(limb1, limb2, op1, op2) \
if (limb1 op1 limb2)                                    \
	return 1;                                           \
else if(limb1 op2 limb2)                                \
	return 0;

#define BINARYCONTAINER_COMPARE_MASKED(limb1, limb2, mask, op1, op2)\
if ((limb1&mask) op1 (limb2&mask))                                  \
	return 1;                                                       \
else if((limb1&mask) op2 (limb2&mask))                              \
	return 0;                                                       \


template<unsigned int length=1024, typename LT=uint64_t>
class BinaryContainer {
public:
	typedef LT   LimbType;
	typedef bool DataType;

private:
	constexpr static LimbType ONE = LimbType(-1);

	// TODO replace __builtin_popcountll with arbitrary correct function for different types
	INLINE constexpr static uint64_t popcount(LimbType x) { return __builtin_popcountll(x); }
	INLINE constexpr static uint64_t popcount(LimbType x, LimbType mask) { return __builtin_popcountll(x & mask); }

public:
	// how many limbs to we need and how wide are they.
	INLINE constexpr static uint16_t limb_bits_width() { return limb_bytes_width() * 8; };
	INLINE constexpr static uint16_t limb_bytes_width() { return sizeof(LimbType); };

private:
	// DO NOT CALL THIS FUNCTION. Use 'limbs()'.
	INLINE constexpr static uint16_t compute_limbs() {
#ifdef BINARY_CONTAINER_ALIGNMENT
		return (alignment()+limb_bits_width()-1)/limb_bits_width();
#else
		return (length+limb_bits_width()-1)/limb_bits_width();
#endif
	};

public:

	/// default constructor
	constexpr BinaryContainer(): __data() { ASSERT(length > 0); }

	/// Copy Constructor
	constexpr BinaryContainer(const BinaryContainer& a) : __data(a.__data) {}


	// round a given amount of 'in' bits to the nearest limb excluding the the lowest overflowing bits
	// eg 13 -> 64
	INLINE constexpr static uint16_t round_up(uint16_t in) { return round_up_to_limb(in) * limb_bits_width(); }
	INLINE constexpr static uint16_t round_up_to_limb(uint16_t in) {return (in/limb_bits_width())+1; }

	// the same as above only rounding down
	// 13 -> 0
	INLINE constexpr static uint16_t round_down(uint16_t in) { return round_down_to_limb(in) * limb_bits_width(); }
	INLINE constexpr static uint16_t round_down_to_limb(uint16_t in) { return (in/limb_bits_width()); }

	// calculate from a bit position 'i' the mask to set it.
	INLINE constexpr static LimbType mask(uint16_t i) {
		ASSERT(i <= length && "wrong access index");
		LimbType u = i%limb_bits_width();
		return (LimbType(1) << u);
	}

	// same as the function belowe, but catches the special case when i == 0 %64.
	INLINE constexpr static LimbType lower_mask2(const uint16_t i) {
		ASSERT(i <= length);
		LimbType u = i%limb_bits_width();
		if (u == 0) return LimbType(-1);
		return ((LimbType(1) << u) - 1);
	}

	// given the i-th bit this function will return a bits mask where the lower 'i' bits are set. Everything will be
	// realigned to limb_bits_width().
	INLINE constexpr static LimbType lower_mask(const uint16_t i) {
		ASSERT(i <= length);
		return ((LimbType(1) << (i%limb_bits_width())) - 1);
	}

	// given the i-th bit this function will return a bits mask where the higher (n-i)bits are set.
	INLINE constexpr static LimbType higher_mask(const uint16_t i) {
		ASSERT(i <= length);
		return (~((LimbType(1) << (i%limb_bits_width())) - 1));
	}

	// given the i-th bit this function will return a bits mask where the lower 'n-i' bits are set. Everything will be
	// realigned to limb_bits_width().
	INLINE constexpr static LimbType lower_mask_inverse(const uint16_t i) {
		ASSERT(i <= length && "wrong access index");
		LimbType u = i%limb_bits_width();

		if (u == 0)
			return -1;

		auto b = (LimbType(1) << (limb_bits_width()-u));
		return b - LimbType(1);
	}

	// given the i-th bit this function will return a bits mask where the higher (i) bits are set.
	INLINE constexpr static LimbType higher_mask_inverse(const uint16_t i) {
		ASSERT(i <= length && "wrong access index");
		return (~lower_mask_inverse(i));
	}

	// not shifted.
	INLINE constexpr LimbType get_bit(const uint16_t i) const {
		return __data[round_down_to_limb(i)] & mask(i);
	}

	// shifted.
	INLINE constexpr bool get_bit_shifted(const uint16_t i) const {
		return (__data[round_down_to_limb(i)] & mask(i)) >> i;
	}

	// return the bits [i,..., j) in one limb
	INLINE LimbType get_bits(const uint16_t i, const uint16_t j) const {
		ASSERT(j > i && j-i <= limb_bits_width() && j <= length);

/*		//TODO not fully correct if both are in the same limb
		const LimbType lmask = higher_mask(i);
		const LimbType rmask = lower_mask2(j);
		const int64_t lower_limb = i/limb_bits_width();
		const int64_t higher_limb = (j-1)/limb_bits_width();
		const uint64_t shift = i%limb_bits_width();
		__uint128_t tmp  = __data[lower_limb] & lmask;
		const uint64_t shift2 = (higher_limb - lower_limb) << 6;
		tmp ^= __uint128_t(__data[higher_limb] & rmask) << (shift2);
		return uint64_t(tmp >> shift);
*/

		const LimbType lmask = higher_mask(i);
		const LimbType rmask = lower_mask2(j);
		const int64_t lower_limb = i/limb_bits_width();
		const int64_t higher_limb = (j-1)/limb_bits_width();

		const uint64_t shift = i%limb_bits_width();
		if (lower_limb == higher_limb) {
			return (__data[lower_limb] & lmask & rmask) >> (shift);
		} else {
			const LimbType a = __data[lower_limb] & lmask;
			const LimbType b = __data[higher_limb] & rmask;

			auto c = (a >> shift);
			auto d = (b << (limb_bits_width()-shift));
			auto r = c ^ d;
			return  r ;
		}
	}

	/// call like this
	/// const LimbType lmask = higher_mask(i);
	/// const LimbType rmask = lower_mask2(j);
	/// const int64_t lower_limb = i / limb_bits_width();
	/// const int64_t higher_limb = (j - 1) / limb_bits_width();
	/// const uint64_t shift = i % limb_bits_width();
	inline LimbType get_bits(const uint64_t llimb, const uint64_t rlimb, const uint64_t lmask, const uint64_t rmask, const uint64_t shift) const {
		//TODO ASSERT(j > i && j - i <= limb_bits_width() && j <= length);
		/*__uint128_t tmp = __data[llimb] & lmask;
		tmp ^= (__uint128_t(__data[rlimb] & rmask) << ((llimb != rlimb) << 6));
		return uint64_t(tmp >> shift);*/
		if (llimb == rlimb) {
			return (__data[llimb] & lmask & rmask) >> (shift);
		} else {
			const LimbType a = __data[llimb] & lmask;
			const LimbType b = __data[rlimb] & rmask;

			auto c = (a >> shift);
			auto d = (b << (limb_bits_width()-shift));
			auto r = c ^ d;
			return r;
		}
	}

	INLINE constexpr void __write_bit(const uint16_t pos, const uint16_t spot, const uint8_t bit) {
		__data[pos] &= ((ONE << spot) | (LT(bit) << (spot)));
	}

	INLINE constexpr void write_bit(const uint16_t pos, const uint8_t bit) {
		__write_bit(pos/limb_bits_width(), pos%limb_bits_width(), bit);
	}

	INLINE constexpr void set_bit(const uint16_t pos) {
		__data[round_down_to_limb(pos)] |= (LimbType(1) << (pos));
	}

	INLINE constexpr void flip_bit(const uint16_t pos) {
		__data[round_down_to_limb(pos)] ^= (LimbType(1) << (pos));
	}

	INLINE constexpr void clear_bit(const uint16_t pos) {
		__data[round_down_to_limb(pos)] &= ~(LimbType(1) << (pos));
	}

	/// zero the complete data vector
	constexpr void zero() {
		LOOP_UNROLL();
		for (unsigned int i = 0; i < limbs(); ++i) {
			__data[i] = 0;
		}
	}

	constexpr void zero(const uint16_t k_lower, const uint16_t k_upper) {
		const uint16_t lower = round_down_to_limb(k_lower);
		const uint16_t upper = round_down_to_limb(k_upper);

		const LimbType lmask = higher_mask(k_lower);
		const LimbType umask = k_upper%limb_bits_width() == 0
								? LimbType(0)
								: lower_mask(k_upper);

		if (lower == upper) {
			const LimbType mask = ~(lmask&umask);
			__data[lower] &= mask;
			return ;
		}

		__data[lower] &= ~lmask;
		__data[upper] &= ~umask;

		LOOP_UNROLL();
		for (uint16_t i = lower+1; i < upper; ++i) {
			__data[i] = 0;
		}
	}

	// seth the whole array to 'fff...fff'
	void one() {
		for (int i = 0; i < limbs(); ++i) {
			__data[i] = ~(__data[i] & 0);
		}
	}

	void one(const uint16_t k_lower, const uint16_t k_upper) {
		const uint16_t lower = round_down_to_limb(k_lower);
		const uint16_t upper = round_down_to_limb(k_upper);

		const LimbType lmask = higher_mask(k_lower);
		const LimbType umask = k_upper%limb_bits_width() == 0
									? LimbType(0)
									: lower_mask(k_upper);

		if (lower == upper) {
			const LimbType mask = (lmask&umask);
			__data[lower] |= ONE & mask;
			return ;
		}

		__data[lower] |= ONE & lmask;
		__data[upper] |= ONE & umask;

		for (int i = lower+1; i < upper; ++i) {
			__data[i] = ~(__data[i] & 0);
		}
	}

	LimbType static random_limb() {
		return fastrandombytes_uint64();
	}

	LimbType static random_limb(const uint16_t k_lower, const uint16_t k_upper) {
		ASSERT(k_upper > k_lower && k_upper - k_lower < limb_bits_width());

		const LimbType lmask = higher_mask(k_lower);
		const LimbType umask = k_upper%limb_bits_width() == 0
									? LimbType(-1) :
									lower_mask(k_upper%limb_bits_width());

		return fastrandombytes_uint64() & lmask & umask;
	}

	/// split the full length BinaryContainer into `k` windows. Inject in every window weight `w` on random positions.
	void random_with_weight_per_windows(const uint16_t w, const uint16_t k) {
		std::vector<uint64_t> buckets_windows{};

		// this stupid approach needs to be done, because if w is not dividing n the last bits would be unused.
		buckets_windows.resize(k+1);
		for (uint16_t i = 0; i < k; ++i) {
			buckets_windows[i] = i*length/k;
		}
		buckets_windows[k] = length;

		// clear everything.
		zero();

		// for every window.
		for (int i = 0; i < k; ++i) {
			uint64_t cur_offset = buckets_windows[i];
			uint64_t windows_length = buckets_windows[i+1] - buckets_windows[i];

			for (int j = 0; j < w; ++j) {
				write_bit(cur_offset + j, true);
			}

			// now permute
			for (int l = 0; l < windows_length; ++l) {
				uint64_t pos = random_limb() % (windows_length - l);
				auto t = get_bit_shifted(cur_offset + l);
				write_bit(cur_offset + l, get_bit_shifted(cur_offset+l+pos));
				write_bit(cur_offset+l+pos, t);
			}
		}
	}

	void random_with_weight(const uint64_t w){
		zero();

		for (int i = 0; i < w; ++i) {
			write_bit(i, true);
		}

		// now permute
		for (int i = 0; i < length; ++i) {
			uint64_t pos = random_limb() % (length - i);
			bool t = get_bit_shifted(i);
			write_bit(i, get_bit_shifted(i+pos));
			write_bit(i+pos, t);
		}
	}

	/// set the whole data array on random data.
	INLINE void random() {
		constexpr uint64_t apply_mask = length%limb_bits_width()==0 ? lower_mask(length)-1 : lower_mask(length);

		if constexpr (length < limb_bits_width()) {
			__data[0] = fastrandombytes_uint64() & apply_mask;
		} else {
			for (int i = 0; i < limbs()-1; ++i) {
				__data[i] = fastrandombytes_uint64();
			}
			__data[limbs()-1] = fastrandombytes_uint64() & apply_mask;
		}
	}

	bool is_zero() const {
		for (int i = 0; i < limbs(); ++i) {
			if(__data[i] != 0)
				return false;
		}

		return true;
	}

	// returns the position in which bits are set.
	void get_bits_set(uint32_t *P, const uint16_t p=1) {
		uint16_t ctr = 0;
		for (uint16_t i = 0; i < length; ++i) {
			if (get_bit(i)) {
				P[ctr++] = i;
				if(ctr == p)
					return;
			}
		}
	}

#ifdef M4RI_M4RI_H
	// M4RI (method Of The 4 Russians) glue code.
	// export/import function
	void to_m4ri(word *p){
		p = __data.data();
	}
	void column_from_m4ri(const mzd_t *H, const uint32_t col){
		ASSERT(H->ncols > col);
		for (int i = 0; i < H->nrows; ++i) {
			write_bit(i, mzd_read_bit(H, i, col));
		}
	}
	void from_m4ri(const word *p){
		for (int i = 0; i < limbs(); ++i) {
			__data[i] = p[i];
		}
	}
	void to_m4ri(mzd_t *p){
		ASSERT(p->nrows == 1 && p->ncols > 0);
		to_m4ri(p->rows[0]);
	}
	void from_m4ri(const mzd_t *p){
		ASSERT(p->nrows == 1 && p->ncols > 0);
		from_m4ri(p->rows[0]);
	}
#endif

	// swap the two bits i, j
	INLINE constexpr void swap(const uint16_t i, const uint16_t j) {
		ASSERT(i < length && j < length);
		auto t = get_bit_shifted(i);
		write_bit(i, get_bit_shifted(j));
		write_bit(j, t);
	}

	INLINE constexpr void flip(const uint16_t i) {
		ASSERT(i < length);
		__data[round_down_to_limb(i)] ^= mask(i);
	}

	INLINE constexpr void neg(const uint16_t k_lower, const uint16_t k_upper) {
		// do nothing.
	}

	// full length addition.
	INLINE int add(BinaryContainer const &v) {
		return add(*this, *this, v);
	}

	// windowed addition.
	INLINE void add(BinaryContainer const &v, const uint16_t k_lower, const uint16_t k_upper) {
		add(*this, *this, v, k_lower, k_upper);
	}

	// TODO templating
	constexpr static void add(LimbType *v5, LimbType const *v1, LimbType const *v2, LimbType const *v3, LimbType const *v4,const uint32_t limbs) {
		int64_t i = 0;

#ifdef USE_AVX2
		LOOP_UNROLL()
		for (; i+4 <= limbs; i+=4) {
			// we need to access the memory unaligned
			__m256 x_avx = MM256_LOAD((float*)v1 + 2*i);
			__m256 y_avx = MM256_LOAD((float*)v2 + 2*i);
			__m256 y1_avx = MM256_LOAD((float*)v3 + 2*i);
			__m256 y2_avx = MM256_LOAD((float*)v4 + 2*i);
			__m256 z_avx = _mm256_xor_ps(x_avx, y_avx);
			_mm256_xor_ps(z_avx, y1_avx);
			_mm256_xor_ps(z_avx, y2_avx);

			MM256_STORE((float*)v5 + 2*i, z_avx);
		}
#endif

		for (; i < limbs; ++i) {
			v5[i] = v1[i] ^ v2[i] ^ v3[i] ^ v4[i];
		}
	}


	constexpr static void add(LimbType *v4, LimbType const *v1, LimbType const *v2, LimbType const *v3, const uint32_t limbs) {
		int64_t i = 0;

#ifdef USE_AVX2
		LOOP_UNROLL()
		for (; i+4 <= limbs; i+=4) {
			// we need to access the memory unaligned
			__m256 x_avx = MM256_LOAD((float*)v1 + 2*i);
			__m256 y_avx = MM256_LOAD((float*)v2 + 2*i);
			__m256 y1_avx = MM256_LOAD((float*)v3 + 2*i);
			__m256 z_avx = _mm256_xor_ps(x_avx, y_avx);
			_mm256_xor_ps(z_avx, y1_avx);
			MM256_STORE((float*)v4 + 2*i, z_avx);
		}
#endif

		for (; i < limbs; ++i) {
			v4[i] = v1[i] ^ v2[i] ^ v3[i];
		}
	}

	// full length addition
	constexpr static void add(LimbType *v3, LimbType const *v1, LimbType const *v2, const uint32_t limbs) {
		int64_t i = 0;

#ifdef USE_AVX2
		LOOP_UNROLL()
		for (; i+4 <= limbs; i+=4) {
			// we need to access the memory unaligned
			__m256 x_avx = MM256_LOAD((float*)v1 + 2*i);
			__m256 y_avx = MM256_LOAD((float*)v2 + 2*i);
			__m256 z_avx = _mm256_xor_ps(x_avx, y_avx);
			MM256_STORE((float*)v3 + 2*i, z_avx);
		}
#endif

		for (; i < limbs; ++i) {
			v3[i] = v1[i] ^ v2[i];
		}
	}

	// full length addition
	constexpr static void add(LimbType *v3, LimbType const *v1, LimbType const *v2) {
		int64_t i = 0;

#ifdef USE_AVX2
		LOOP_UNROLL()
		for (; i+4 < limbs(); i+=4) {
			// we need to access the memory unaligned
			__m256 x_avx = MM256_LOAD((float*)v1 + 2*i);
			__m256 y_avx = MM256_LOAD((float*)v2 + 2*i);
			__m256 z_avx = _mm256_xor_ps(x_avx, y_avx);
			MM256_STORE((float*)v3 + 2*i, z_avx);
		}
#endif

		for (; i < limbs(); ++i) {
			v3[i] = v1[i] ^ v2[i];
		}
	}

	// full length addition
	//  IMPORTANT: this function does a full length addition
	constexpr static int add(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2) {
		constexpr uint64_t upper = limbs();
		int64_t i = 0;
#ifdef USE_AVX2
		LOOP_UNROLL()
		for (; i+4 < upper; i+=4) {
			// we need to access the memory unaligned
			__m256 x_avx = MM256_LOAD((float*)v1.__data.data() + 2*i);
			__m256 y_avx = MM256_LOAD((float*)v2.__data.data() + 2*i);
			__m256 z_avx = _mm256_xor_ps(x_avx, y_avx);
			MM256_STORE((float*)v3.__data.data() + 2*i, z_avx);
		}
#endif
		LOOP_UNROLL()
		for (; i < upper; ++i) {
			v3.__data[i] = v1.__data[i] ^ v2.__data[i];
		}

		return 0;
	}

	// TODO Tests
	template<const uint32_t llimb, const uint32_t ulimb, const LimbType lmask, const LimbType rmask>
	constexpr static void add(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2) {
		if constexpr (llimb == ulimb) {
			constexpr LimbType mask = (lmask & rmask);
			LimbType tmp1 = (v3.__data[llimb] & ~(mask));
			LimbType tmp2 = (v1.__data[llimb] ^ v2.__data[llimb]) & mask;
			v3.__data[llimb] = tmp1 ^ tmp2;
			return;
		}

		LOOP_UNROLL();
		for (int64_t i = llimb+1; i < ulimb; ++i) {
			v3.__data[i] = v1.__data[i] ^ v2.__data[i];
		}

		LimbType tmp1 = (v1.__data[llimb] ^ v2.__data[llimb]) & lmask;
		LimbType tmp2 = (v1.__data[ulimb] ^ v2.__data[ulimb]) & rmask;
		LimbType tmp11 = (v3.__data[llimb] & ~(lmask));
		LimbType tmp21 = (v3.__data[ulimb] & ~(rmask));

		v3.__data[llimb] = tmp1^tmp11;
		v3.__data[ulimb] = tmp2^tmp21;
	}

	// TODO test
	template<const uint32_t llimb, const uint32_t ulimb, const LimbType lmask, const LimbType rmask>
	constexpr static void add(LimbType *v3, LimbType const *v1, LimbType const *v2) {
		if constexpr (llimb == ulimb) {
			constexpr LimbType mask = (lmask & rmask);
			LimbType tmp1 = (v3[llimb] & ~(mask));
			LimbType tmp2 = (v1[llimb] ^ v2[llimb]) & mask;
			v3[llimb] = tmp1 ^ tmp2;
			return;
		}

		LOOP_UNROLL();
		for (int64_t i = llimb+1; i < ulimb; ++i) {
			v3[i] = v1[i] ^ v2[i];
		}

		LimbType tmp1 = (v1[llimb] ^ v2[llimb]) & lmask;
		LimbType tmp2 = (v1[ulimb] ^ v2[ulimb]) & rmask;
		LimbType tmp11 = (v3[llimb] & ~(lmask));
		LimbType tmp21 = (v3[ulimb] & ~(rmask));

		v3[llimb] = tmp1^tmp11;
		v3[ulimb] = tmp2^tmp21;
	}


	constexpr static void add(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2,
	                          const uint16_t k_lower, const uint16_t k_upper) {
		ASSERT(k_upper <= length && k_lower < k_upper);
		const LimbType lmask        = higher_mask(k_lower%limb_bits_width());
		const LimbType rmask        = lower_mask2(k_upper%limb_bits_width());
		const int64_t lower_limb    = k_lower / limb_bits_width();
		const int64_t higher_limb   = (k_upper-1) / limb_bits_width();

		if (lower_limb == higher_limb) {
			const LimbType mask = k_upper%limb_bits_width() == 0 ? lmask : (lmask & rmask);
			LimbType tmp1 = (v3.__data[lower_limb] & ~(mask));
			LimbType tmp2 = (v1.__data[lower_limb] ^ v2.__data[lower_limb]) & mask;
			v3.__data[lower_limb] = tmp1 ^ tmp2;
			return;
		}

		LOOP_UNROLL();
		for (int64_t i = lower_limb+1; i < higher_limb; ++i) {
			v3.__data[i] = v1.__data[i] ^ v2.__data[i];
		}

		LimbType tmp1 = (v1.__data[lower_limb] ^ v2.__data[lower_limb]) & lmask;
		LimbType tmp2 = (v1.__data[higher_limb] ^ v2.__data[higher_limb]) & rmask;
		LimbType tmp11 = (v3.__data[lower_limb] & ~(lmask));
		LimbType tmp21 = (v3.__data[higher_limb] & ~(rmask));

		v3.__data[lower_limb] = tmp1^tmp11;
		v3.__data[higher_limb]= tmp2^tmp21;
	}

	constexpr static uint32_t add_weight(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2,
	                                const uint16_t k_lower, const uint16_t k_upper) {
		ASSERT(k_upper <= length && k_lower < k_upper && 0 < k_upper);

		uint32_t cnorm = 0;
		const LimbType lmask        = higher_mask(k_lower%limb_bits_width());
		const LimbType rmask        = lower_mask2(k_upper%limb_bits_width());
		const int64_t lower_limb    = k_lower/limb_bits_width();
		const int64_t higher_limb   = (k_upper-1)/limb_bits_width();

		if (lower_limb == higher_limb) {
			const LimbType mask = k_upper%limb_bits_width() == 0 ? lmask : (lmask & rmask);
			LimbType tmp1 = (v3.__data[lower_limb] & ~(mask));
			LimbType tmp2 = (v1.__data[lower_limb] ^ v2.__data[lower_limb]) & mask;
			v3.__data[lower_limb] = tmp1 ^ tmp2;
			auto b = popcount(tmp2);
			return b;
		}

		LOOP_UNROLL();
		for (int64_t i = lower_limb+1; i < higher_limb; ++i) {
			v3.__data[i] = v1.__data[i] ^ v2.__data[i];
			cnorm += popcount(v3.__data[i]);
		}

		LimbType tmp1 = (v1.__data[lower_limb] ^ v2.__data[lower_limb]) & lmask;
		LimbType tmp2 = (v1.__data[higher_limb] ^ v2.__data[higher_limb]) & rmask;
		LimbType tmp11 = (v3.__data[lower_limb] & ~(lmask));
		LimbType tmp21 = (v3.__data[higher_limb] & ~(rmask));

		v3.__data[lower_limb] = tmp1^tmp11;
		v3.__data[higher_limb]= tmp2^tmp21;

		cnorm += popcount(tmp1);
		cnorm += popcount(tmp2);

		return cnorm;
	}

	INLINE constexpr static bool add(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2,
	                                 const uint16_t k_lower, const uint16_t k_upper,
	                                 const uint32_t norm) {
		ASSERT(k_upper <= length && k_lower < k_upper);

		if (norm == uint32_t(-1)){
			// fallback to normal addition.
			add(v3, v1, v2, k_lower, k_upper);
			return false;
		} else {
			ASSERT((&v1 != &v3 && &v2 != &v3) || (norm == -1));
			uint32_t cnorm = add_weight(v3, v1, v2, k_lower, k_upper);
			if (cnorm >= norm)
				return true;
		}

		return false;
	}


	// full length subtraction=addition in F_2
	INLINE int sub(BinaryContainer const &v) {
		return this->add(v);
	}

	/// alias for add
	INLINE void sub(BinaryContainer const &v, const uint16_t k_lower, const uint16_t k_upper) {
		return add(v, k_lower, k_upper);
	}

	/// alias for add
	INLINE constexpr static int sub(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2) {
		return add(v3, v1, v2);
	}

	/// alias for add
	INLINE constexpr static bool sub(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2,
	                                 const uint16_t k_lower, const uint16_t k_upper) {
		add(v3, v1, v2, k_lower, k_upper);
		return false; //TODO not implemented
	}

	/// alias for add
	INLINE static bool sub(BinaryContainer &v3, BinaryContainer const &v1, BinaryContainer const &v2,
	                       const uint16_t k_lower, const uint16_t k_upper,
	                       const uint32_t norm) {
		return add(v3, v1, v2, k_lower, k_upper, norm);
	}

	INLINE constexpr static bool cmp(BinaryContainer const &v1, BinaryContainer const &v2) {
		return cmp(v1, v2, 0, ssize());
	}

	/// implements only a 2 way comparison. E.g. implements the `!=` operator.
	INLINE constexpr static bool cmp(BinaryContainer const &v1, BinaryContainer const &v2,
	                                 const uint16_t k_lower, const uint16_t k_upper) {
		ASSERT(k_upper <= length && k_lower < k_upper);
		int64_t lower = round_down_to_limb(k_lower);
		int64_t upper = round_down_to_limb(k_upper-1);
		const LimbType lmask = higher_mask(k_lower);
		const LimbType rmask = lower_mask2(k_upper);

		if (lower == upper) {   // the two offsets lay in the same limb.
			const LimbType mask = k_upper%limb_bits_width() == 0 ? lmask : (lmask & rmask);
			return cmp_simple2(v1, v2, lower, mask);
		} else {                // the two offsets lay in two different limbs
			// first check the highest limb with the mask
			return cmp_ext2(v1, v2, lower, upper, lmask, rmask);
		}
	}

	/// Important: lower != higher
	/// unrolled high speed implementation of a multi limb compare function
	template<const uint16_t lower, const uint16_t upper, const LimbType lmask, const LimbType umask>
	INLINE constexpr static bool cmp_ext(BinaryContainer const &v1, BinaryContainer const &v2) {
		ASSERT(lower != upper && lower < upper);

		// first check the highest limb with the mask
		if ((v1.__data[upper]&umask) != (v2.__data[upper]&umask))
			return false;

		// check all limbs in the middle
		LOOP_UNROLL()
		for(uint16_t i = upper-1; i > lower; i--) {
			if (v1.__data[i] != v2.__data[i])
				return false;
		}

		if ((v1.__data[lower]&lmask) != (v2.__data[lower]&lmask))
			return false;
		return true;
	}

	/// IMPORTANT; lower < upper so you have to compare at least two limbs.
	/// use this function if you have to compare a lot of different elements on the same coordinate. So you can precompute
	/// the mask and the limbs.
	INLINE constexpr static bool cmp_ext2(BinaryContainer const &v1, BinaryContainer const &v2,
	                                      const uint16_t lower, const uint16_t upper,
	                                      const LimbType lmask, const LimbType umask) {
		ASSERT(lower < upper && lmask != 0 && upper < limbs());
		// first check the highest limb with the mask.
		if ((v1.__data[upper]&umask) != (v2.__data[upper]&umask))
			return false;

		// check all limbs in the middle.
		for(uint64_t i = upper-1; i > lower; i--) {
			if (v1.__data[i] != v2.__data[i])
				return false;
		}

		// and at the end check the lowest limb.
		if ((v1.__data[lower]&lmask) != (v2.__data[lower]&lmask))
			return false;
		return true;
	}

	/// IMPORTANT: lower != higher => mask != 0. This is actually only a sanity check.
	/// high speed implementation of a same limb cmompare function
	template<const uint16_t limb, const LimbType mask>
	INLINE constexpr static bool cmp_simple(BinaryContainer const &v1, BinaryContainer const &v2) {
		ASSERT(limb != uint64_t(-1) && mask != 0);
		return ((v1.__data[limb]&mask) == (v2.__data[limb]&mask));
	}

	/// IMPORTANT: mask != 0.
	/// use this function if you have to compare a lot of different elements on the same coordinate. So you can precompute
	/// the mask and the limb.
	INLINE constexpr static bool cmp_simple2(BinaryContainer const &v1, BinaryContainer const &v2,
										  const uint16_t limb, const LimbType mask) {
		ASSERT(limb != uint64_t(-1) && mask != 0);
		return ((v1.__data[limb]&mask) == (v2.__data[limb]&mask));
	}

	INLINE constexpr static int cmp_ternary_simple2(BinaryContainer const &v1, BinaryContainer const &v2,
												 const uint16_t limb, const LimbType mask) {
		ASSERT(limb != uint64_t(-1));
		if ((v1.__data[limb]&mask) > (v2.__data[limb]&mask))
			return 1;
		else if ((v1.__data[limb]&mask) < (v2.__data[limb]&mask))
			return -1;

		return 0;
	}

	/// IMPORTANT: k_lower < k_upper is enforced.
	/// sets v1 = v2[k_lower, ..., k_upper].
	/// Does not change anything else.
	inline constexpr static void set(BinaryContainer &v1, BinaryContainer const &v2,
								  const uint16_t k_lower, const uint16_t k_upper) {
		ASSERT(k_upper <= length && k_lower < k_upper);
		const int64_t lower = round_down_to_limb(k_lower);
		const int64_t upper = round_down_to_limb(k_upper-1);
		const LimbType lmask = higher_mask(k_lower);
		const LimbType rmask = lower_mask2(k_upper);

		if (lower == upper) {   // the two offsets lay in the same limb.
			const LimbType mask = k_upper%limb_bits_width() == 0 ? lmask : (lmask & rmask);
			v1.__data[lower] = (v1.__data[lower] & ~mask) | (v2.__data[lower] & mask);
			return;
		} else {                // the two offsets lay in two different limbs
			v1.__data[upper] = (v1.__data[upper] & ~rmask) | (v2.__data[upper] & rmask);
			v1.__data[lower] = (v1.__data[lower] & ~lmask) | (v2.__data[lower] & lmask);
			for(uint64_t i = upper-1; i > lower; i--) {
				v1.data()[i] = v2.data()[i];
			}
		}
	}

	///  out[s: ] = in[0:s]
	INLINE static void shift_right(BinaryContainer &out, const BinaryContainer &in, const uint16_t s) {
		for (int j = 0; j < s; ++j) {
			out.write_bit(j+s, in.get_bit_shifted(j));
		}
	}

	/// checks whether this == obj on the interval [k_lower, ..., k_upper]
	/// the level of the calling 'list' object.
	/// \return
	INLINE bool is_equal(const BinaryContainer &obj, const uint16_t k_lower, const uint16_t k_upper) const {
		return cmp(*this, obj, k_lower, k_upper);
	}

	template<const uint16_t lower, const uint16_t upper, const LimbType lmask, const LimbType umask>
	INLINE bool is_equal_ext(const BinaryContainer &obj) const {
		return cmp_ext2<lower, upper, lmask, umask>(*this, obj);
	}

	INLINE bool is_equal_ext2(const BinaryContainer &obj, const uint16_t lower, const uint16_t upper,
	                          const LimbType lmask, const LimbType umask) const {
		return cmp_ext2(*this, obj, lower, upper, lmask, umask);
	}

	INLINE bool is_equal_simple2(const BinaryContainer &obj, const uint16_t limb, const LimbType mask) const {
		return cmp_simple2(*this, obj, limb, mask);
	}

	/// implements a strict comparison. Call this function if you dont know what to call. Its the most generic implementaion
	/// and it works for all input.s
	INLINE bool is_greater(BinaryContainer const &obj, const uint16_t k_lower, const uint16_t k_upper) const {
		ASSERT(k_upper <= length && k_lower < k_upper);
		int64_t lower = round_down_to_limb(k_lower);
		int64_t upper = round_down_to_limb(k_upper-1);
		const LimbType lmask = higher_mask(k_lower);
		const LimbType rmask = lower_mask2(k_upper);

		if (lower == upper) {   // the two offsets lay in the same limb.
			const LimbType mask = k_upper%limb_bits_width() == 0 ? lmask : (lmask & rmask);
			return this->is_greater_simple2(obj, lower, mask);
		} else {                // the two offsets lay in two different limbs
			return this->is_greater_ext2(obj, lower, upper, lmask, rmask);
		}
	}

	/// return *this > obj on the limbs [lower, upper]. `lmask` and `umask` are bitmask for the lowest and highest limbs.
	/// Technically this is a specially unrolled implementation of `BinaryContainer::is_greater` if you have to compare a lot
	/// of containers repeatedly on the same coordinates.
	INLINE bool is_greater_ext2(BinaryContainer const &obj, const uint16_t lower, const uint16_t upper,
	                            const LimbType lmask, const LimbType umask) const {
		ASSERT(lower < upper && lmask != 0 && upper < limbs());
		// umask is allowed to be zero. Otherwise cases like k_upper = 128 wouldn't make sense.

		BINARYCONTAINER_COMPARE_MASKED(__data[upper], obj.__data[upper], umask, >, <)
		// check all limbs in the middle
		for(uint64_t i = upper-1; i > lower; i--) {
			BINARYCONTAINER_COMPARE(__data[i], obj.__data[i], >, <)
		}

		BINARYCONTAINER_COMPARE_MASKED(__data[lower], obj.__data[lower], lmask, >, <)
		return false;
	}

	INLINE bool is_greater_equal_ext2(BinaryContainer const &obj, const uint16_t lower, const uint16_t upper,
	                                  const LimbType lmask, const LimbType umask) const {
		ASSERT(lower < upper && lmask != 0 && upper < limbs());
		// umask is allowed to be zero. Otherwise cases like k_upper = 128 wouldn't make sense.
		BINARYCONTAINER_COMPARE_MASKED(__data[upper], obj.__data[upper], umask, >=, <)
		// check all limbs in the middle
		for(uint64_t i = upper-1; i > lower; i--) {
			BINARYCONTAINER_COMPARE(__data[i], obj.__data[i], >=, <)
		}

		BINARYCONTAINER_COMPARE_MASKED(__data[lower], obj.__data[lower], lmask, >=, <)
		return false;
	}

	/// most simple type of comparison implemented for this class.
	/// returns *this < obj on bits specified by the parameter `limb` and `mask.`
	/// call like this:
	///		using BinaryContainerTest = BinaryContainer<64>;
	///		BinaryContainerTest b1, b2;
	///		uint64_t limb = 0;
	///		mask = BinaryContainerTest::higher_mask(k_lower) & BinaryContainerTest::lower_mask2(k_higher);
	///		b1.is_greater_simple2(b2, limb, mask);
	INLINE bool is_greater_simple2(BinaryContainer const &obj, const uint16_t limb, const LimbType mask) const {
		ASSERT(limb < limbs() && mask != 0);
		return ((__data[limb]&mask) > (obj.__data[limb]&mask));
	}

	/// not testet
	INLINE bool is_greater_equal_simple2(BinaryContainer const &obj, const uint16_t limb, const LimbType mask) const {
		ASSERT(limb < limbs() && mask != 0);
		return ((__data[limb]&mask) >= (obj.__data[limb]&mask));
	}

	/// main comparison function for the < operator. If you dont know what function to use, use this one. Its the most generic
	/// implementation and works for all inputs.
	INLINE bool is_lower(BinaryContainer const &obj, const uint16_t k_lower, const uint16_t k_upper) const {
		ASSERT(k_upper <= length && k_lower < k_upper);
		int64_t lower = round_down_to_limb(k_lower);
		int64_t upper = round_down_to_limb(k_upper-1);
		const LimbType lmask = higher_mask(k_lower);
		const LimbType rmask = lower_mask2(k_upper);

		if (lower == upper) {   // the two offsets lay in the same limb.
			const LimbType mask = k_upper%limb_bits_width() == 0
									? lmask
									: (lmask & rmask);
			return is_lower_simple2(obj, lower, mask);
		} else {                // the two offsets lay in two different limbs
			return is_lower_ext2(obj, lower, upper, lmask, rmask);
		}
	}

	/// return *this < obj on the limbs [lower, upper]. `lmask` and `umask` are bitmask for the lowest and highest limbs.
	/// Technically this is a specially unrolled implementation of `BinaryContainer::is_lower` if you have to compare a lot
	/// of containers repeatedly on the same coordinates.
	/// Example Code:
	///		using BinaryContainerTest = BinaryContainer<G_n>;
	///
	///		BinaryContainerTest2 b1, b2;
	///			FILL WITH DATA HERE.
	///		const uint64_t lower = BinaryContainerTest2::round_down_to_limb(k_higher);
	///		const uint64_t upper = BinaryContainerTest2::round_down_to_limb(b1.size()-1);
	///		const BinaryContainerTest2::LimbType lmask = BinaryContainerTest2::higher_mask(k_higher);
	///		const BinaryContainerTest2::LimbType umask = BinaryContainerTest2::lower_mask2(b1.size());
	///		is_lower_ext2(b2, lower, upper, lmask, umask))
	/// You __MUST__ be extremely carefully with the chose of `upper`.
	/// Normally you want to compare two elements between `k_lower` and `k_upper`. This can be done by:
	///		const uint64_t lower = BinaryContainerTest2::round_down_to_limb(k_lower);
	///		const uint64_t upper = BinaryContainerTest2::round_down_to_limb(k_higher);
	///							... (as above)
	/// Note that you dont have to pass the -1 to k_higher. The rule of thumb is that you __MUST__ add a -1 to the computation
	/// of the upper limb.
	INLINE bool is_lower_ext2(BinaryContainer const &obj, const uint16_t lower, const uint16_t upper,
	                          const LimbType lmask, const LimbType umask) const {
		ASSERT(lower < upper && lmask != 0 && upper < limbs());
		// umask is allowed to be zero. Otherwise cases like k_upper = 128 wouldnt make sense.

		BINARYCONTAINER_COMPARE_MASKED(__data[upper], obj.__data[upper], umask, <, >)
		// check all limbs in the middle
		for(uint64_t i = upper-1; i > lower; i--) {
			BINARYCONTAINER_COMPARE(__data[i], obj.__data[i], <, >)
		}

		BINARYCONTAINER_COMPARE_MASKED(__data[lower], obj.__data[lower], lmask, <, >)
		return false;
	}

	/// not testet
	INLINE bool is_lower_equal_ext2(BinaryContainer const &obj, const uint16_t lower, const uint16_t upper,
	                                const LimbType lmask, const LimbType umask) const {
		ASSERT(lower < upper && lmask != 0 && upper < limbs());
		// umask is allowed to be zero. Otherwise cases like k_upper = 128 wouldnt make sense.

		BINARYCONTAINER_COMPARE_MASKED(__data[upper], obj.__data[upper], umask, <=, >)
		// check all limbs in the middle
		for(uint64_t i = upper-1; i > lower; i--) {
			BINARYCONTAINER_COMPARE(__data[i], obj.__data[i], <=, >)
		}

		BINARYCONTAINER_COMPARE_MASKED(__data[lower], obj.__data[lower], lmask, <=, >)
		return false;
	}

	// efficient reimplementation of `is_lower` for the special case that `k_lower` and `k_upper` are in the same limb.
	/// call like this:
	///		using BinaryContainerTest = BinaryContainer<64>;
	///		BinaryContainerTest b1, b2;
	///			FILL WITH DATA
	///		uint64_t limb = 0;
	///		mask = BinaryContainerTest::higher_mask(k_lower) & BinaryContainerTest::lower_mask2(k_higher);
	///		b1.is_lower_simple2(b2, limb, mask);
	INLINE bool is_lower_simple2(BinaryContainer const &obj, const uint16_t limb, const LimbType mask) const {
		ASSERT(limb < limbs() < length && mask != 0);
		return ((__data[limb]&mask) < (obj.__data[limb]&mask));
	}

	INLINE bool is_lower_equal_simple2(BinaryContainer const &obj, const uint16_t limb, const LimbType mask) const {
		ASSERT(limb < limbs() < length && mask != 0);
		return ((__data[limb]&mask) <= (obj.__data[limb]&mask));
	}

	// calculates the weight = hamming weight of the data container.
	inline uint16_t weight() const {
		return weight(0, length);
	}

	uint32_t weight(const uint16_t k_lower, const uint16_t k_upper) const {
		ASSERT(k_upper <= length && k_lower < k_upper);
		const uint64_t lower = round_down_to_limb(k_lower);
		const uint64_t upper = round_down_to_limb(k_upper-1);

		const LimbType l_mask = higher_mask(k_lower);
		const LimbType u_mask = lower_mask2(k_upper);

		uint64_t weight;
		// if only one limb needs to be checked to check
		if(lower == upper){
			uint64_t b = (l_mask & u_mask);
			uint64_t c = uint64_t(__data[lower]);
			uint64_t d = uint64_t(b) & uint64_t(c);
			uint64_t w_ = popcount(d);
			return w_;
		}

		weight = popcount(l_mask&__data[lower]);
		weight += popcount(u_mask&__data[upper]);
		for(uint64_t i=lower+1;i<upper;++i)
			weight += popcount(__data[i]);

		return weight;
	}

	template<const uint16_t lower, const uint16_t upper, const LimbType l_mask, const LimbType u_mask>
	constexpr uint32_t weight() const {
		ASSERT(lower <= upper);
		uint32_t weight = 0;

		// if only one limb needs to be checked to check
		if constexpr (lower == upper){
			uint64_t b = (l_mask & u_mask);
			uint64_t c = uint64_t(__data[lower]);
			uint64_t d = uint64_t(b) & uint64_t(c);
			uint64_t w_ = popcount(d);
			return w_;
		}

		weight = popcount(l_mask&__data[lower]);
		for(uint32_t i=lower+1;i<upper;++i)
			weight += popcount(__data[i]);
		weight += popcount(u_mask&__data[upper]);
		return weight;
	}



	INLINE constexpr static void _not(BinaryContainer &out, const BinaryContainer &in);
	INLINE constexpr static void _not(BinaryContainer &out, const BinaryContainer &in,
									 const uint16_t lower, const uint16_t higher);

	constexpr void _not();
	constexpr void _not(const uint16_t lower, const uint16_t higher);


	constexpr static void _and(BinaryContainer &out, const BinaryContainer & in1, const BinaryContainer & in2);
	constexpr static void _and(BinaryContainer &out, const BinaryContainer & in1, const BinaryContainer & in2,
	                             const uint16_t lower, const uint16_t higher);

	void _and(BinaryContainer &out, const BinaryContainer &in);
	void _and(BinaryContainer &out, const BinaryContainer &in,
			const uint16_t lower, const uint16_t higher);


	constexpr static void _xor(BinaryContainer &out, const BinaryContainer & in1,const  BinaryContainer & in2);
	constexpr static void _xor(BinaryContainer &out, const BinaryContainer & in1,const  BinaryContainer & in2,
	                           const uint16_t lower, const uint16_t higher);

	void _xor(BinaryContainer &out, const BinaryContainer &in);
	void _xor(BinaryContainer &out, const BinaryContainer &in,
	          const uint16_t lower, const uint16_t higher);


	constexpr static void diff(BinaryContainer &out, const BinaryContainer & in1,const  BinaryContainer & in2);
	constexpr static void diff(BinaryContainer &out, const BinaryContainer & in1,const  BinaryContainer & in2,
	                           const uint16_t lower, const uint16_t higher);

	void diff(BinaryContainer &out, const BinaryContainer &in);
	void diff(BinaryContainer &out, const BinaryContainer &in,
	          const uint16_t lower, const uint16_t higher);

	constexpr static int setMaskRow(BinaryContainer &in, const size_t bitlen, const size_t len);
	int setMaskPositionRow(BinaryContainer &in, const size_t bitlen, const size_t pos, const size_t len);


	constexpr static void lshdRow(BinaryContainer &out, const BinaryContainer &a, const size_t x);
	void lshdRowInplace(BinaryContainer &a, const size_t x);

	constexpr static void lshRow(BinaryContainer &out, const BinaryContainer &c, const size_t b, const size_t len);
	void lshRowInplace(BinaryContainer &a, const size_t x, const size_t len);

	constexpr static void rshdRow(BinaryContainer &out, const BinaryContainer &a, const size_t x);
	void rshdRowInplace(BinaryContainer &a, const size_t x);

	constexpr static void rshRow(BinaryContainer &out, const BinaryContainer &c, const size_t b, const size_t len);
	void rshRowInplace(BinaryContainer &a, const size_t x, const size_t len);




	constexpr static int clearRow(BinaryContainer &row);
	int clearRowN(BinaryContainer &row, const size_t len);

	int setRowUInt32(BinaryContainer &row, const uint32_t n);
	int setRowUInt64(BinaryContainer &row, const uint64_t n);

	constexpr static uint32_t getFirstBitSet(const BinaryContainer &row, const size_t bitlen);
	constexpr static uint32_t getLastBitSet(const BinaryContainer &row, const size_t bitlen);


	// hack it like C++
	// Took from glibc++
	class reference {
		friend class BinaryContainer;

		// pointer to the limb
		LimbType     *wp;
		// bit position in the whole data array.
		const size_t 	        mask_pos;

		// left undefined
		reference();

	public:
		reference(BinaryContainer &b, size_t pos) : mask_pos(mask(pos)){
			wp = &b.data().data()[round_down_to_limb(pos)];
		}

#if __cplusplus >= 201103L
		reference(const reference&) = default;
#endif

		~reference() = default;

		// For b[i] = __x;
		reference& operator=(bool x) {
			if (x)
				*wp |= mask_pos;
			else
				*wp &= ~mask_pos;
			return *this;
		}

		// For b[i] = b[__j];
		reference& operator=(const reference& j) {
			if (*(j.wp) & j.mask_pos)
				*wp |= mask_pos;
			else
				*wp &= ~mask_pos;
			return *this;
		}

		// Flips the bit
		bool operator~() const { return (*(wp) & mask_pos) == 0; }

		// For __x = b[i];
		operator bool() const {
			return (*(wp) & mask_pos) != 0;
		}

		// For b[i].flip();
		reference& flip() {
			*wp ^= mask_pos;
			return *this;
		}

		unsigned int get_data() const { return bool(); }
		unsigned int data() const { return bool(); }
	};
	friend class reference;


	reference operator[](size_t pos) { return reference(*this, pos); }
	constexpr bool operator[](const size_t pos) const { return (__data[round_down_to_limb(pos)] & mask(pos)) != 0; }

	/// Assignment operator implementing copy assignment
	/// see https://en.cppreference.com/w/cpp/language/operators
	/// \param obj
	/// \return
	BinaryContainer& operator =(BinaryContainer const &obj) {
		if (this != &obj) { // self-assignment check expected
			std::copy(&obj.__data[0], &obj.__data[0] + obj.__data.size(), &this->__data[0]);
		}

		return *this;
	}

	/// Assignment operator implementing move assignment
	/// Alternative definition: Value& operator =(Value &&obj) = default;
	/// see https://en.cppreference.com/w/cpp/language/move_assignment
	/// \param obj
	/// \return
	BinaryContainer& operator =(BinaryContainer &&obj) noexcept {
		if (this != &obj) { // self-assignment check expected really?
			// move the data
			__data = std::move(obj.__data);
		}

		return *this;
	}

	void print(const uint64_t k_lower, const uint64_t k_upper) {
		ASSERT(k_lower < length && k_upper < length && k_lower < k_upper);
		for (uint64_t i = k_lower; i < k_upper; ++i) {
			std::cout << __data[i] << " ";
		}
		std::cout << "\n";
	}

	//LimbType data(uint64_t index) { ASSERT(index < length); return get_bit_shifted(index); }
	bool data(uint64_t index) const { ASSERT(index < length); return get_bit_shifted(index); }

	LimbType get_type() {return __data[0]; }

#ifdef BINARY_CONTAINER_ALIGNMENT
	static constexpr uint16_t alignment() {
		// Aligns to a multiple of 32 Bytes
		//constexpr uint16_t bytes = (length+7)/8;
		//constexpr uint16_t limbs = (((bytes-1) >> 5)+1)<<5;
		//return limbs*8;

		// Aligns to a multiple of 16 Bytes
		constexpr uint16_t bytes = (length+7)/8;
		constexpr uint16_t limbs = (((bytes-1) >> 4)+1)<<4;
		return limbs*8;
	}
#endif

	// returns size in bytes
	static constexpr uint64_t copyable_ssize() {
#ifdef BINARY_CONTAINER_ALIGNMENT
		return alignment()/8;
#else
		return (length+7)/8;
#endif
	}
	static constexpr uint64_t ssize() { return length; }
	inline constexpr uint64_t size() const { return length; }
	inline constexpr static uint64_t limbs() { return (length+limb_bits_width()-1)/limb_bits_width(); }
	inline constexpr uint64_t get_limbs() { return (length+limb_bits_width()-1)/limb_bits_width(); }

	LimbType* ptr() { return __data.data(); };

	std::array<LimbType, compute_limbs()>& data() {
		return __data;
	};
	const std::array<LimbType, compute_limbs()>& data() const {
		return __data;
	};

private:
	// actual data container.
	std::array<LimbType, compute_limbs()> __data;
};


template<unsigned int length=1024, typename LT=uint64_t>
std::ostream& operator<< (std::ostream &out, const BinaryContainer<length> &obj) {
	for (uint64_t i = 0; i < obj.size(); ++i) {
		out << obj[i];
	}
	return out;
}

#endif //BINARY_ROW_H
