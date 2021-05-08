#include "binary_container.h"

#include <limits>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <cstdint>

uint64_t get_random_uint64_t() {
	static std::random_device dev;
	static std::mt19937 rng(dev());
	static std::uniform_int_distribution<uint64_t> dist(1,  std::numeric_limits<uint64_t>::max());

	return dist(rng);
}

template<typename digit>
int binary_row<digit>::resizeRow(digit *in, const size_t new_bin_len, const size_t old_binLen){
    //TODO
    return MA_OK;
}


/// Set the last 'bitmask' Bits to 1;
// 'bitLen' is the actual binary len of the input. NOT the number of limbs
/// \tparam digit
/// \param in
/// \param bitmask
/// \param bitLen
/// \return
template<typename digit>
int binary_row<digit>::setMaskRow(digit *in, const size_t bitmask, const size_t bitLen){
    if (unlikely(bitmask >= bitLen))
        return MA_VAL;

    // we need to clear the complete row.
    clearRow(in);

    digit mask = 0u;
    mask = ~(mask & 0u);
    size_t limbs = bitmask / BITS_OF_DIGIT;
    digit bit = (((digit)(1u) << (bitmask % BITS_OF_DIGIT)) - 1);

    // set the first limbs completely to 1;
    for (size_t i = 0; i < limbs; ++i){
        in[i] = mask;
    }

    // Set the last limb
    in[limbs] = bit;

    return MA_OK;
}


/// Set the last 'bitmask' Bits to 1 starting from the bit 'pos';
/// \param in
/// \param bitmask      number of bts to 1 starting from `pos`
/// \param pos          starting bit of the mask
/// \param bitLen       actual bit len, not number of limbs
/// \return
template<typename digit>
int binary_row<digit>::setMaskPositionRow(digit *in, const size_t bitmask, const size_t pos, const size_t bitLen){
    if (unlikely((bitmask + pos) >= bitLen))
        return MA_VAL;

    // we need to clear the complete row.
    clearRow(in);


    // TODO optimize
    for (size_t i = pos; i < pos + bitmask; i++){
        set_bit(in, i, 1);
    }

    /*digit mask = 0u;
    mask = ~(mask & 0u);        // set everything to 1

    // get the number of bits set in the 'pos' limb
    digit startingMask = pos;
    for (; startingMask > 0 ; startingMask -= BITS_OF_DIGIT);

    // calculate the bit masks for the first an last limb
    startingMask = (mask << (BITS_OF_DIGIT - startingMask));
    digit endMask = startingMask; // TODO not correct

    // how many limbs between the first an last do we need to set completely to `1...1`
    size_t limbs = bitmask / BITS_OF_DIGIT;
    digit bit = (((digit)(1) << (bitmask % BITS_OF_DIGIT)) - 1);

    // set the limbs
    size_t n = pos/BITS_OF_DIGIT;
    in[n] = startingMask;
    for (size_t i = n + 1; i < n + limbs; ++i){
        in[i] = mask;
    }
    in[n+limbs] = endMask;*/

    return MA_OK;
}

///
///
/// \param out
/// \param in
template<typename digit>
void binary_row<digit>::notRow(digit *__restrict__ out, const digit *__restrict__ in) {
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(int i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(in + 32*i));
        STOREBYPASS256((__m256i *)(out + i*32), _mm256_not_si256(a));
    }
    // dont need to the rest, because we make sure, that number of limbs is always div by 32
#else
    LOOP_UNROLL()
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, !get_bit(in, i));
    }
#endif

}

///
/// \param out
template<typename digit>
void binary_row<digit>::notRowInplace(digit *out) {
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(size_t i = 0; i < NUMBER_OF_ROW_LIMBS_256; i += 1){
        __m256i a = LOADBYPASS256((__m256i *)(out + 32*i));
        STOREBYPASS256((__m256i *)(out + 32*i), _mm256_not_si256(a));
    }

    // dont need to the rest, because we make sure, that number of limbs is always div by 32
#else
    LOOP_UNROLL()
    for(size_t i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, !get_bit(out, i));
    }
#endif

}

///
/// \param out = !in
/// \param in
/// \param bitLen actual length of the in vector in bits. not limbs
template<typename digit>
void binary_row<digit>::notRowN(digit *__restrict__ out, const digit *__restrict__ in, const size_t bitLen) {
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen/256; i += 32){
        __m256i a = LOADBYPASS256((__m256i *)(in + i));
        STOREBYPASS256((__m256i *)(out + i), _mm256_not_si256(a));
    }

    // TODO this can be optimized for COMPRESSED MEM
    for(; i < bitLen; i++){
        set_bit(out, i, !get_bit(in, i));
    }
#else
    LOOP_UNROLL()
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, !get_bit(in, i));
    }
#endif

}

///
/// \param out
/// \param bitLen  actual length in bits. Not in Limbs
template<typename digit>
void binary_row<digit>::notRowNInplace(digit *out, const size_t bitLen) {
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen/256; i += 32){
        digit *p = out + i;
        __m256i a = LOADBYPASS256((__m256i *)p);
        STOREBYPASS256((__m256i *)p, _mm256_not_si256(a));
        //printf("AVX used\n");
    }

    //printf("i: %d; N: %d\n", i, bitLen);
    // TODO this can be optimized
    for(; i < bitLen; ++i){
        set_bit(out, i, !get_bit(out, i));
    }
#else
    LOOP_UNROLL()
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, !get_bit(out, i));
    }
#endif

}

///
/// \tparam digit
/// \param __restrict__
/// \param __restrict__
/// \param __restrict__
template<typename digit>
void binary_row<digit>::andRow(digit *__restrict__ out, const digit *__restrict__ in1, const digit *__restrict__ in2){
#ifdef USE_AVX2
    LOOP_UNROLL()
    for (int i = 0; i < BYTES_IN_ROW; i += 32) {
        __m256i a = LOADBYPASS256((__m256i*)&in1[i]);
        __m256i b = LOADBYPASS256((__m256i*)&in2[i]);

        //__m256i a = LOADBYPASS256((__m256i *)(in1 + 32*i));
        //__m256i b = LOADBYPASS256((__m256i *)(in2 + 32*i));
        __m256i c = _mm256_and_si256(a, b);

        //STOREBYPASS256((__m256i *)(out + i*32), c);
        STOREBYPASS256((__m256i*)&out[i], c);
    }
#else
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, get_bit(in1, i) & get_bit(in2, i));
    }
#endif
}

///
/// \tparam digit
/// \param __restrict__
/// \param __restrict__
template<typename digit>
void binary_row<digit>::andRowInplace(digit *__restrict__ out, const digit *__restrict__ in){
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(size_t i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(out + 32*i));
        __m256i b = LOADBYPASS256((__m256i *)(in + 32*i));
        __m256i c = _mm256_and_si256(a, b);

        STOREBYPASS256((__m256i *)(out + i*32), c);
    }

#else
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, get_bit(out, i) & get_bit(in, i));

    }
#endif
}


/// 'bitLen' is the actual binary len of the input which should be 'and'. NOT the number of limbs
/// \tparam digit
/// \param out
/// \param in1
/// \param in2
/// \param bitLen
template<typename digit>
void binary_row<digit>::andRowN(digit *__restrict__ out, const digit *__restrict__ in1, const digit *__restrict__ in2, const size_t bitLen){
#ifdef USE_AVX2
    size_t limit_byte = bitLen >> 3u ;
    size_t i = 0;

    LOOP_UNROLL()
    for(; i + BYTES_PER_LIMBS_256 <= limit_byte; i += BYTES_PER_LIMBS_256){
        __m256i a = LOADBYPASS256((__m256i *)(in1 + i));
        __m256i b = LOADBYPASS256((__m256i *)(in2 + i));

        STOREBYPASS256((__m256i *)(out + i), _mm256_and_si256(a, b));
    }

    //TODO here maybe here sse3

    LOOP_UNROLL()
    for(; i + BYTES_PER_LIMBS_64 <= limit_byte; i += BYTES_PER_LIMBS_64){
        uint64_t tmp = ((uint64_t)in1[i] & (uint64_t)in2[i]);
        ((uint64_t *)&out[i])[0] = tmp;
    }

    LOOP_UNROLL()
    for(; i + BYTES_PER_LIMBS_32 <= limit_byte; i += BYTES_PER_LIMBS_32){
        uint32_t tmp = ((uint32_t)(in1[i]) & (uint32_t)in2[i]);
        ((uint32_t *)&out[i])[0] = tmp;

        //printf("whup 32\n");
    }

    //printf("aaab %d %d %d\n", limit_byte,  i, bitLen);


    LOOP_UNROLL()
    for(; i + BYTES_PER_LIMBS_8 <= limit_byte; i += BYTES_PER_LIMBS_8){
        uint8_t tmp = ((uint8_t)in1[i] & (uint8_t)in2[i]);
        ((uint8_t *)&out[i])[0] = tmp;
    }

    //printf("aaa %d %d %d\n", limit_byte,  i, bitLen);

    i = i*8;

    //printf("aaa %d %d %d\n", limit_byte,  i, bitLen);

    for(; i < bitLen; i++){
        set_bit(out, i, get_bit(in1, i) & get_bit(in2, i));
    }
#else
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, get_bit(in1, i) & get_bit(in2, i));
    }
#endif
}


///
/// \tparam digit
/// \param __restrict__
/// \param __restrict__
/// \param bitLen  'bitLen' is the actual binary len of the input. NOT the number of limbs
template<typename digit>
void binary_row<digit>::andRowNInplace(digit *__restrict__ out, const digit *__restrict__ in, const size_t bitLen){
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen / 32; i += 32){
        __m256i a = LOADBYPASS256((__m256i *)(out + i));
        __m256i b = LOADBYPASS256((__m256i *)(in +  i));
        STOREBYPASS256((__m256i *)(out + i), _mm256_and_si256(a, b));
    }
        // TODO this can be optimized for COMPRESSED MEM
    for(; i < bitLen; i++){
        set_bit(out, i, get_bit(out, i) & get_bit(in, i));
    }
#else
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, get_bit(out, i) & get_bit(in, i));
    }
#endif
}

///
/// \tparam digit
/// \param out
/// \param __restrict__
/// \param __restrict__
template<typename digit>
void binary_row<digit>::xorRow(digit *out, const digit *__restrict__ in1, const digit *__restrict__ in2){
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(size_t i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(in1 + 32*i));
        __m256i b = LOADBYPASS256((__m256i *)(in2 + 32*i));
        __m256i c = _mm256_xor_si256(a, b);

        STOREBYPASS256((__m256i *)(out + i*32), c);
    }
#else
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, get_bit(in1, i) ^ get_bit(in2, i));
    }
#endif
}

template<typename digit>
void binary_row<digit>::xorRowInplace(digit *out, const digit *in){
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(size_t i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(in + 32*i));
        __m256i b = LOADBYPASS256((__m256i *)(out + 32*i));
        __m256i c = _mm256_xor_si256(a, b);

        STOREBYPASS256((__m256i *)(out + i*32), c);
    }
#else
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, get_bit(out, i) ^ get_bit(in, i));
    }
#endif
}

///
/// \tparam digit
/// \param out
/// \param __restrict__
/// \param __restrict__
/// \param bitLen bitLen' is the actual binary len of the input. NOT the number of limbs
template<typename digit>
void binary_row<digit>::xorRowN(digit *out, const digit *__restrict__ in1, const digit *__restrict__ in2, const size_t bitLen){
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen / 256; i++){
        __m256i *inn1 = ((__m256i *)in1) + i;
        __m256i *inn2 = ((__m256i *)in2) + i;
        __m256i *outt = ((__m256i *)out) + i;

        __m256i a = LOADBYPASS256(inn1);
        __m256i b = LOADBYPASS256(inn2);

        STOREBYPASS256(outt, _mm256_xor_si256(a, b));
    }

    // TODO this can be optimized for COMPRESSED MEM
    for(; i < bitLen; i++){
        set_bit(out, i, get_bit(in1, i) ^ get_bit(in2, i));
    }
#else
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, get_bit(in1, i) ^ get_bit(in2, i));
    }
#endif
}


///
/// \tparam digit
/// \param out
/// \param in
/// \param bitLen bitLen' is the actual binary len of the input. NOT the number of limbs
template<typename digit>
void binary_row<digit>::xorRowNInplace(digit *out, const digit *in, const size_t bitLen){
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen / 256; i++){
        __m256i *inn = ((__m256i *)in) + i;
        __m256i *outt = ((__m256i *)out) + i;

        __m256i a = LOADBYPASS256(inn);
        __m256i b = LOADBYPASS256(outt);

        STOREBYPASS256(outt, _mm256_xor_si256(a, b));
    }

    // TODO this can be optimized for COMPRESSED MEM
    for(; i < bitLen; i++){
        set_bit(out, i, get_bit(out, i) ^ get_bit(in, i));
    }
#else
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, get_bit(out, i) ^ get_bit(in, i));
    }
#endif
}

/// calculate the difference vector between 'in1' and 'in2'. So e.g.
///         out = in1 & ~in2
/// \param out
/// \param in1
/// \param in2
template<typename digit>
void binary_row<digit>::diffRow(digit *out, const digit *__restrict__ in1, const digit *__restrict__ in2){
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(size_t i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(in1 + 32*i));
        __m256i b = _mm256_not_si256(LOADBYPASS256((__m256i *)(in2 + 32*i)));
        __m256i c = _mm256_and_si256(a, b);

        STOREBYPASS256((__m256i *)(out + i*32), c);
    }
#else
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, get_bit(in1, i) &  (~get_bit(in2, i)));
    }
#endif
}

/// same as 'diffRow'
/// \param out
/// \param in
template<typename digit>
void binary_row<digit>::diffRowInplace(digit *out, const digit *in){
#ifdef USE_AVX2
    LOOP_UNROLL()
    for(size_t i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = _mm256_not_si256(LOADBYPASS256((__m256i *)(in + 32*i)));
        __m256i b = LOADBYPASS256((__m256i *)(out + 32*i));
        __m256i c = _mm256_and_si256(a, b);

        STOREBYPASS256((__m256i *)(out + i*32), c);
    }
#else
    for(int i = 0; i < BITS_IN_ROW; i++){
        set_bit(out, i, get_bit(out, i) & (~get_bit(in, i)));
    }
#endif
}

///
/// \param out
/// \param in1
/// \param in2
/// \param bitLen
template<typename digit>
void binary_row<digit>::diffRowN(digit *out, const digit *__restrict__ in1, const digit *__restrict__ in2, const size_t bitLen){
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen / 256; i++){
        __m256i *inn1 = ((__m256i *)in1) + i;
        __m256i *inn2 = ((__m256i *)in2) + i;
        __m256i *outt = ((__m256i *)out) + i;

        __m256i a = _mm256_not_si256(LOADBYPASS256(inn1));
        __m256i b = LOADBYPASS256(inn2);

        STOREBYPASS256(outt, _mm256_and_si256(a, b));
    }

    // TODO this can be optimized for COMPRESSED MEM
    for(; i < bitLen; i++){
        set_bit(out, i, get_bit(in1, i) &  (~get_bit(in2, i)));
    }
#else
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, get_bit(in1, i) &  (~get_bit(in2, i)));
    }
#endif
}

///
/// \param out
/// \param in
/// \param bitLen
template<typename digit>
void binary_row<digit>::diffRowNInplace(digit *out, const digit *in, const size_t bitLen){
#ifdef USE_AVX2
    size_t i = 0;

    LOOP_UNROLL()
    for(; i < bitLen / 256; i++){
        __m256i *inn = ((__m256i *)in) + i;
        __m256i *outt = ((__m256i *)out) + i;

        __m256i a = _mm256_not_si256(LOADBYPASS256(inn));
        __m256i b = LOADBYPASS256(outt);

        STOREBYPASS256(outt, _mm256_and_si256(a, b));
    }

    // TODO this can be optimized for COMPRESSED MEM
    for(; i < bitLen; i++){
        set_bit(out, i, get_bit(out, i) & (~get_bit(in, i)));
    }
#else
    for(int i = 0; i < bitLen; i++){
        set_bit(out, i, get_bit(out, i) & (~get_bit(in, i)));
    }
#endif
}


///
/// Shift x*BITS_OF_DIGIT bits: out = (a<<x*(BITS_OF_DIGIT))
/// 'bitLen' is not  parameter, because its quite usless.
///
/// \tparam digit
/// \param out
/// \param a
/// \param x
template<typename digit>
void binary_row<digit>::lshdRow(digit *out, const digit *a, const size_t x){
    if (unlikely(x > NUMBER_OF_ROW_LIMBS)){
	    printf("to many limbs\n");
        return;
    }

    if (x == 0)
        return cpyRow(out, a);

    int y = NUMBER_OF_ROW_LIMBS;
    for (; y >= x; y--) {
        out[y] = a[y-x];
    }
    for (; y >= 0; y--) {
        out[y] = 0;
    }

}


/// Shift x*BITS_OF_DIGIT bit
/// \tparam digit
/// \param a
/// \param x
template<typename digit>
void binary_row<digit>::lshdRowInplace(digit *a, const size_t x){
    if (unlikely(x > NUMBER_OF_ROW_LIMBS)){
        printf("to many limbs");
        return;
    }

    if (x == 0)
        return;

    int y = NUMBER_OF_ROW_LIMBS;
    for (; y >= x; y--) {
        a[y] = a[y-x];
    }

    for (; y >= 0; y--) {
        a[y] = 0;
    }

}

/// shifts bits. out = c << b
/// \param out
/// \param in
/// \param b
/// \param len  = bitLen
template<typename digit>
void binary_row<digit>::lshRow(digit *out, const digit *in, const size_t b, const size_t len){
    size_t inword_shift = b % BITS_OF_DIGIT;
    size_t extra_words = b / BITS_OF_DIGIT;
    const size_t as = len/BITS_OF_DIGIT;    // so this should be limbs

    if (b == 0)
        return cpyRow(out, in);

    // fast path
    if (inword_shift == 0) {
        return lshdRow(out, in, extra_words);
    }

    out[as] = in[as+extra_words] >> (BITS_OF_DIGIT - inword_shift);
    for (size_t i = as; i >= 2; i--){
        out[i-1] = (in[i - 1 - extra_words] << inword_shift)
                    | (in[i - 2 - extra_words] >> (BITS_OF_DIGIT - inword_shift));
    }

    out[extra_words] = in[0] << inword_shift;

    for(size_t i = 0; i < extra_words; i++) {
        out[i] = 0;
    }

}



/// c <<= b
/// \tparam digit
/// \param c
/// \param b
/// \param len
template<typename digit>
void binary_row<digit>::lshRowInplace(digit *c, const size_t b, const size_t len){
    size_t inword_shift = b % BITS_OF_DIGIT;
    size_t extra_words = b / BITS_OF_DIGIT;
    const size_t as = len/BITS_OF_DIGIT;    // so this should be limbs

    // fast path
    if (inword_shift == 0) {
        return lshdRowInplace(c, extra_words);
    }

    c[as] = c[as+extra_words] >> (BITS_OF_DIGIT - inword_shift);
    for (size_t i = as; i >= 2; i--){
        c[i-1] = (c[i - 1 - extra_words] << inword_shift)
                 | (c[i - 2 - extra_words] >> (BITS_OF_DIGIT - inword_shift));
    }

    //TODO eine resizeable version machen in der die BitMaps dynamisch wachsen
    /*c[as] = c[as+extra_words] >> (BITS_OF_DIGIT - inword_shift);
    for (size_t i = as+extra_words; i >= as + 2; i--){
        c[i-1] = (c[i - 1 - extra_words] << inword_shift)
                   | (c[i - 2 - extra_words] >> (BITS_OF_DIGIT - inword_shift));
    }*/

    c[extra_words] = c[0] << inword_shift;

    for(size_t i = 0; i < extra_words; i++) {
        c[i] = 0;
    }
}


/// rightshift limbs NOT bits
///     out = a << x
/// \param out
/// \param a
/// \param x
template<typename digit>
void binary_row<digit>::rshdRow(digit *out, const digit *a, const size_t x){
    if (unlikely(x > NUMBER_OF_ROW_LIMBS)){
        printf("to many limbs");
        return;
    }

    // fast path
    if (x == 0)
        return cpyRow(out, a);

    int y = NUMBER_OF_ROW_LIMBS;

    LOOP_UNROLL()
    for (; y >= NUMBER_OF_ROW_LIMBS-x; y--) {
        out[y] = 0;
    }

    LOOP_UNROLL()
    for (; y >= 0; y--) {
        out[y] = a[y+x];
    }
}

///
/// \param a
/// \param x
template<typename digit>
void binary_row<digit>::rshdRowInplace(digit *a, const size_t x){
    if (unlikely(x > NUMBER_OF_ROW_LIMBS)){
        printf("to many limbs");
        return;
    }

    LOOP_UNROLL()
    for (int y = NUMBER_OF_ROW_LIMBS; y >= x; y--) {
        a[y-x] = a[y];
    }

    LOOP_UNROLL()
    for (int y = NUMBER_OF_ROW_LIMBS; y >= NUMBER_OF_ROW_LIMBS-x; y--) {
        //a[y] = 0;
    }
}

/// shifts bits. out = c >> b
/// \param out
/// \param in
/// \param b
/// \param len  = bitLen
template<typename digit>
void binary_row<digit>::rshRow(digit *out, const digit *in, const size_t b, const size_t len){
    printf("NOT IMPLEMETED\n");
}


///  c >>= b
/// \tparam digit
/// \param c
/// \param b
/// \param len
template<typename digit>
void binary_row<digit>::rshRowInplace(digit *c, const size_t b, const size_t len){
    printf("Not Implemented\n");
}

///
/// \param out
/// \param in
template<typename digit>
void binary_row<digit>::cpyRow(digit *__restrict__ out, const digit *__restrict__ in){
#ifdef USE_AVX2
    memcpy(out, in, NUMBER_OF_ROW_LIMBS*BYTES_OF_DIGIT);

    /*
    LOOP_UNROLL()
    for(int i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(in + 32*i));
        STOREBYPASS256((__m256i *)(out + i*32), a);
    }
     */
#else
    memcpy(out, in, COLUMNS);
#endif

}


///
/// \param out
/// \param in
/// \param bitLen actually bitlen not number of limbs
template<typename digit>
void binary_row<digit>::cpyRowN(digit *__restrict__ out, const digit *__restrict__ in, const size_t bitLen){
    size_t limbs = bitLen/BITS_OF_DIGIT;
    memcpy(out, in, limbs);

    // Do the rest
    for (size_t i = limbs * BITS_OF_DIGIT; i < bitLen; i++) {
        set_bit(out, i, get_bit(in, i));
    }
}

///
/// \param in1
/// \param in2
/// \return
template<typename digit>
int binary_row<digit>::cmpRow(const digit *__restrict__ in1, const digit *__restrict__ in2){
    return memcmp(in1, in2, BYTES_IN_ROW);
}

/// dont forget to prefetch
/// \tparam digit
/// \param __restrict__
/// \param __restrict__
/// \param bitLen
/// \return
template<typename digit>
int binary_row<digit>::cmpRowN(const digit *__restrict__ in1, const digit *__restrict__ in2, const size_t bitLen){
    // 8 because of bytes
    size_t l = bitLen / 8;
    int r =  memcmp(in1, in2, l);
    //printf("l: %zu, %zu\n", l, bitLen);

    if (r != 0)
        return r;

    // Now we need to check the last bits
    for (size_t i = l * 8; i < bitLen; ++i) {
        //printf("i: %zu\n", i);
        if (get_bit(in1, i) > get_bit(in2, i))
            return MA_GE;
        if (get_bit(in1, i) < get_bit(in2, i))
            return MA_LE;
    }

    return MA_EQ;
}

template<typename digit>
int binary_row<digit>::randomRow(digit *row){
    getrandom(row, BYTES_IN_ROW, 0);
    return MA_OK;
}


///
/// \tparam digit
/// \param row
/// \param bitLen bitLen = bits to set random
/// \return
template<typename digit>
int binary_row<digit>::randomRowN(digit *row, const size_t bitLen){
    size_t bytes = bitLen / 8;
    getrandom(row, bytes, 0);
    for (size_t i = bytes * 8; i < bitLen; ++i) {
        set_bit(row, i, get_random_uint64_t() & (digit)1u);
    }

    return MA_OK;
}

///
/// \param row
/// \param weightMap
/// \param len_weightMap
/// \param len
/// \return
template<typename digit>
int binary_row<digit>::randomRowWithWeightMap(digit *row, const double *weightMap, const size_t len_weightMap, const size_t len){
    for (int i = 0; i < len; ++i) {
        double x = get_random_uint64_t();

        double savedProb = 1.0;
        digit savedVal  = 0;

        for (size_t val = 0; val < len_weightMap; ++val) {
            double prob = weightMap[val];
            if ((prob <= x) || (prob >= savedProb))
                continue;

            savedProb = prob;
            savedVal = val;
        }

        set_bit(row, i, savedVal & 1u);
    }

    return MA_OK;
}

template<typename digit>
int binary_row<digit>::clearRow(digit *row) {
#ifdef USE_AVX2
    // at this point we dont need to distinguish between compressed mem or not.
    // This is the fastest possible
    memset(row, 0, BYTES_IN_ROW);
#else
    LOOP_UNROLL()
    for (int i = 0; i < BITS_IN_ROW; ++i) {
        set_bit(row, i, 0);
    }
#endif
    return MA_OK;
}

///
/// \param row
/// \param bitLen acutal len in bits. not limbs
/// \return
template<typename digit>
int binary_row<digit>::clearRowN(digit *row, const size_t bitLen) {
    size_t bytes = bitLen / 8;
    memset(row, 0, bytes);

    for (size_t i = bytes*8; i < bitLen; i++){
        set_bit(row, i, 0);
    }
    return MA_OK;
}

/// set row [3][2][1][0] = n    // if |digit| = 8
/// \param row
/// \param n
/// \return
template<typename digit>
int binary_row<digit>::setRowUInt32(digit *row, const uint32_t n) {
#if ((BITS_IN_ROW) < 32)
    #error lol whats wrong with you? less then 32 bits in a row?
#endif
    digit mask = (1u << BITS_OF_DIGIT) - 1u;

    LOOP_UNROLL()
    for(uint32_t i = 0; i < 32/BITS_OF_DIGIT; i ++) {
        row[i] = (n >> (i*BITS_OF_DIGIT)) & mask;
    }

    return MA_OK;
}

/// set row [7]...[3][2][1][0] = n    // if |digit| = 8
/// \param row
/// \param n
/// \return
template<typename digit>
int binary_row<digit>::setRowUInt64(digit *row, const uint64_t n) {
#if ((BITS_IN_ROW) < 64)
    #error lol whats wrong with you? less then 64 bits in a row?
#endif
    digit mask = (1u << BITS_OF_DIGIT) - 1;

    LOOP_UNROLL()
    for(uint64_t i = 0; i < 64/BITS_OF_DIGIT; i++) {
        row[i] = (n >> (i*BITS_OF_DIGIT)) & mask;
    }

    return MA_OK;
}

/// interpret the row as a number and add a 'Integer' to it
/// \param row
/// \param n
template<typename digit>
void binary_row<digit>::addRowInt(digit *row, const uint64_t n, const size_t binlen){
    uint64_t *a = (uint64_t *)row;

    uint64_t loAdd = n;
    a[0] += loAdd;


    for (int i = 0; i < binlen/64 - 1; ++i) {
        a[i+1] += 0 + (a[i] < loAdd);
        loAdd = a[i];
    }
}

/// interpret 'in1' and 'in2' as big int. Add them together and write the result into 'out'
/// \param out
/// \param in1
/// \param in2
/// \param binlen
template<typename digit>
void binary_row<digit>::addRow(digit *__restrict__ out, digit *__restrict__ in1, digit *__restrict__ in2, const size_t binlen){
    const size_t y = binlen/BITS_OF_DIGIT;
    size_t x = 0;
    uint16_t t = 0;

    for (x = 0; x < y; x++) {
        t += ((uint16_t) in1[x]) + ((uint16_t) in2[x]);
        out[x] = (digit) t;
        t >>= BITS_OF_DIGIT;
    }
}

#ifdef USE_AVX2

// Source: https://arxiv.org/pdf/1611.07612.pdf
// A C function using AVX2 intrinsics
// implementing a bitwise parallel carry-save adder (CSA).
void CSA(__m256i *h, __m256i *l, __m256i a, __m256i b, __m256i c) {
    __m256i u = _mm256_xor_si256(a, b);
    *h = _mm256_or_si256(_mm256_and_si256(a, b), _mm256_and_si256(u, c));
    *l = _mm256_xor_si256(u, c);
}

static __m256i popcount_pshufb(__m256i v) {

    const __m256i lookup = _mm256_setr_epi8(
            /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
            /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
            /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
            /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4,

            /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
            /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
            /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
            /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4
    );

    const __m256i low_mask = _mm256_set1_epi8(0x0f);

    const __m256i lo  = _mm256_and_si256(v, low_mask);
    const __m256i hi  = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
    const __m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
    const __m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);

    return _mm256_sad_epu8(_mm256_add_epi8(popcnt1, popcnt2), _mm256_setzero_si256());
}

// implementing Harley-Seal’s algorithm. It assumes, for
// simplicity, that the input size in 256-bit vectors is divisible
// by 16. See Fig. 10 for the count function.
///
/// \param data
/// \param size in bytes
/// \return
static uint64_t count_bits_avx2(const __m256i* data, const uint64_t size) {
    __m256i total     = _mm256_setzero_si256();
    __m256i ones      = _mm256_setzero_si256();
    __m256i twos      = _mm256_setzero_si256();
    __m256i fours     = _mm256_setzero_si256();
    __m256i eights    = _mm256_setzero_si256();
    __m256i sixteens  = _mm256_setzero_si256();
    __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

    const uint64_t limit = size - size % 16;
    uint64_t i = 0;

    for(; i < limit; i += 16) {
        CSA(&twosA, &ones, ones, _mm256_lddqu_si256(data + i), _mm256_lddqu_si256(data + i + 1));
        CSA(&twosB, &ones, ones, _mm256_lddqu_si256(data + i + 2), _mm256_lddqu_si256(data + i + 3));
        CSA(&foursA, &twos, twos, twosA, twosB);
        CSA(&twosA, &ones, ones, _mm256_lddqu_si256(data + i + 4), _mm256_lddqu_si256(data + i + 5));
        CSA(&twosB, &ones, ones, _mm256_lddqu_si256(data + i + 6), _mm256_lddqu_si256(data + i + 7));
        CSA(&foursB,& twos, twos, twosA, twosB);
        CSA(&eightsA,&fours, fours, foursA, foursB);
        CSA(&twosA, &ones, ones, _mm256_lddqu_si256(data + i + 8), _mm256_lddqu_si256(data + i + 9));
        CSA(&twosB, &ones, ones, _mm256_lddqu_si256(data + i + 10), _mm256_lddqu_si256(data + i + 11));
        CSA(&foursA, &twos, twos, twosA, twosB);
        CSA(&twosA, &ones, ones, _mm256_lddqu_si256(data + i + 12), _mm256_lddqu_si256(data + i + 13));
        CSA(&twosB, &ones, ones, _mm256_lddqu_si256(data + i + 14), _mm256_lddqu_si256(data + i + 15));
        CSA(&foursB, &twos, twos, twosA, twosB);
        CSA(&eightsB, &fours, fours, foursA, foursB);
        CSA(&sixteens, &eights, eights, eightsA, eightsB);

        total = _mm256_add_epi64(total, popcount_pshufb(sixteens));
    }

    total = _mm256_slli_epi64(total, 4);     // * 16
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcount_pshufb(eights), 3)); // += 8 * ...
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcount_pshufb(fours),  2)); // += 4 * ...
    total = _mm256_add_epi64(total, _mm256_slli_epi64(popcount_pshufb(twos),   1)); // += 2 * ...
    total = _mm256_add_epi64(total, popcount_pshufb(ones));
    for(; i < size; i++)
        total = _mm256_add_epi64(total, popcount_pshufb(_mm256_lddqu_si256(data + i)));


    return (uint64_t)(_mm256_extract_epi64(total, 0))
           + (uint64_t)(_mm256_extract_epi64(total, 1))
           + (uint64_t)(_mm256_extract_epi64(total, 2))
           + (uint64_t)(_mm256_extract_epi64(total, 3));
}

template<typename digit>
uint64_t binary_row<digit>::weightRow(const digit *row){
    // remember, we need to pass the length in bytes
    return count_bits_avx2((const __m256i *)row, BYTES_IN_ROW);
}

template<typename digit>
uint64_t binary_row<digit>::weightRowN(const digit *row, const size_t bitLen){
    return count_bits_avx2((const __m256i *)row, bitLen/BYTES_OF_DIGIT);
}


#endif


/// returns the position of the first bit set.
/// Set LSB
/// \param row
/// \param bitlen
/// \return
template<typename digit>
uint64_t binary_row<digit>::getFirstBitSet(const digit *row, const size_t bitlen){
    for (size_t i = 0; i < bitlen/BITS_OF_DIGIT; i++){
        if (row[i] != 0){
            return __builtin_clzll(row[i]) + i * BITS_OF_DIGIT;
        }
    }

    return bitlen;
}

/// set MSB
/// \param row
/// \param bitlen
/// \return
template<typename digit>
uint64_t binary_row<digit>::getLastBitSet(const digit *row, const size_t bitlen){
    for (size_t i = bitlen/BITS_OF_DIGIT; i > 0; --i){
        if (row[i-1] != 0){
            return (BITS_OF_DIGIT-1) - __builtin_clzll(row[i-1]) + (i - 1) * BITS_OF_DIGIT;
        }
    }
    return 0;

}

/// check if the row == 0**bitlen
/// \param row
/// \param bitlen
/// \return
template<typename digit>
int binary_row<digit>::anyBitSetRow(const digit *row, const size_t bitlen) {
#ifdef USE_AVX2
    __m256i null = _mm256_set_epi64x(0, 0, 0, 0);

    LOOP_UNROLL()
    for(int i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i b = LOADBYPASS256((__m256i *)(row + 32*i));
        __m256i pcmp = _mm256_cmpeq_epi32(b, null);
        unsigned bitmask = _mm256_movemask_epi8(pcmp);
        if (bitmask != 0xffffffffu){
            return MA_YES;
        }
    }
#else
    for (size_t i = 0; i < bitlen/BITS_OF_DIGIT; i++){
        if (row[i] != 0) {
            return MA_YES;
        }
    }
#endif
    return MA_NO;

}

#ifdef __AVX2__
/// returns true if every bit set in 'in2' is also set in 'in1'
/// \param in1
/// \param in2
/// \param bitlen
/// \return
template<typename digit>
int binary_row_avx<digit>::containsAllRow(const digit *in1, const digit *in2, const size_t bitlen){
    // TODO anpassen für bitlen
    LOOP_UNROLL()
    for(int i = 0; i < NUMBER_OF_ROW_LIMBS_256; i++){
        __m256i a = LOADBYPASS256((__m256i *)(in1 + 32*i));
        __m256i b = LOADBYPASS256((__m256i *)(in2 + 32*i));
        __m256i pcmp = _mm256_cmpeq_epi32(a, _mm256_and_si256(a, b));
        unsigned bitmask = _mm256_movemask_epi8(pcmp);
        if (bitmask != 0xffffffffU){
            return MA_NO;
        }
    }

    return MA_YES;
}
#endif
