#ifndef BINARY_ROW_RANDOM_H
#define BINARY_ROW_RANDOM_H

static uint64_t random_x=123456789u, random_y=362436069u, random_z=521288629u;

// Sowas von nicht sicher, Aber egal.
static void random_seed(uint64_t i){
	random_x += i;
	random_y = random_x*4095834;
	random_z = random_x + random_y*98798234;
}

static uint64_t xorshf96() {          //period 2^96-1
	uint64_t t;
	random_x ^= random_x << 16u;
	random_x ^= random_x >> 5u;
	random_x ^= random_x << 1u;

	t = random_x;
	random_x = random_y;
	random_y = random_z;
	random_z = t ^ random_x ^ random_y;

	return random_z;
}

/* n = size of buffer in bytes, */
static int fastrandombytes(void *buf, size_t n){
	uint64_t *a = (uint64_t *)buf;

	const uint32_t rest = n%8;
	const size_t limit = n/8;
	size_t i = 0;

	for (; i < limit; ++i) {
		a[i] = xorshf96();
	}

	// last limb
	uint8_t *b = (uint8_t *)buf;
	b += n - rest;
	uint64_t limb = xorshf96();
	for (size_t j = 0; j < rest; ++j) {
		b[j] = (limb >> (j*8u)) & 0xFFu;
	}

	return 0;
}

/* n = bytes. */
inline static int fastrandombytes_uint64_array(uint64_t *buf, size_t n){
	void *a = (void *)buf;
	fastrandombytes(a, n);
	return 0;
}

static uint64_t fastrandombytes_uint64() {
	constexpr uint32_t UINT64_POOL_SIZE = 512;    // page should be 512 * 8 Byte
	static uint64_t tmp[UINT64_POOL_SIZE];
	static size_t counter = 0;

	if (counter == 0){
		fastrandombytes_uint64_array(tmp, UINT64_POOL_SIZE * 8 );
		counter = UINT64_POOL_SIZE;
	}

	counter -= 1;
	return tmp[counter];
}

#endif //BINARY_ROW_RANDOM_H
