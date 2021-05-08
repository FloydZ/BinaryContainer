IMPORTANT NOTE
----
This is only a pre alpha version.

Usage
---
```C++
#include "binary_container"
using BC = BinaryContainer<1024, uint64_t>;
BC a,b,c;
    
```

Flags
----
The following compiler flags are supported:
```bash
USE_BRANCH_PREDICTION
USE_PREFETCH
USE_LOOP_UNROLL
BINARY_CONTAINER_ALIGNMENT
```

Additionally if you include [M4RI](https://bitbucket.org/malb/m4ri) (alternatively you can define `M4RI_M4RI_H` ) before add this project supporting functions to convert vectors from or to M4RI.

TODO
----

Implement the following functions
```C++
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
```
