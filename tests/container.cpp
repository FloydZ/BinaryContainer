#include <gtest/gtest.h>
#include <cstdint>
#include <bitset>
#include <utility>

// Hack for testing private functions (C++ god)
#define private public

#include "binary_container.h"
#include "test.h"

using ::testing::EmptyTestEventListener;
using ::testing::InitGoogleTest;
using ::testing::Test;
using ::testing::TestEventListeners;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::UnitTest;

using namespace std;

TEST(Internals, get_limbs) {
	BinaryContainer b;
#if N <= 64
	EXPECT_EQ(1, b.limbs());
#elif N <= 128
	EXPECT_EQ(2, b.limbs());
#else
	EXPECT_EQ((N + 64 -1)/64, b.limbs());
#endif
}

TEST(Internals, get_size) {
	BinaryContainer<> b;
	EXPECT_EQ(N, b.size());
}

TEST(Internals, access) {
	BinaryContainer<> b;
	for (int i = 0; i < b.size(); ++i) {
		// this are the explicit cast steps to the final result.
		auto bit = b[i];
		bool bbit = bool(bit);
		EXPECT_EQ(0, bbit);
	}
}

TEST(Internals, access_pass_through) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero();
	b2.random();

	EXPECT_EQ(b1.size(), b2.size());
	for (int i = 0; i < b1.size(); ++i) {
		b2[i] = b1[i];
	}

	for (int i = 0; i < b1.size(); ++i) {
		EXPECT_EQ(b2[i], b1[i]);
	}
}

TEST(Internals, get_bits) {
	using TestBinaryContainer = BinaryContainer<512>;
	// TODO test get_bits
}


/*
TEST(Internals, round_up_to_limb){
	BinaryContainer<64> b;

	EXPECT_EQ(b.round_up_to_limb(0), 1);
	EXPECT_EQ(b.round_up_to_limb(20), 1);
	EXPECT_EQ(b.round_up_to_limb(63), 1);
	EXPECT_EQ(b.round_up_to_limb(64), 2);
}

TEST(Internals, round_down_to_limb){
	BinaryContainer<64> b;

	EXPECT_EQ(b.round_down_to_limb(0), 0);
	EXPECT_EQ(b.round_down_to_limb(20), 0);
	EXPECT_EQ(b.round_down_to_limb(63), 0);
	EXPECT_EQ(b.round_down_to_limb(64), 1);
}
*/
TEST(Internals, compute_limbs){
	BinaryContainer<63> b1;
	EXPECT_EQ(b1.compute_limbs(), 1);

	BinaryContainer<64> b2;
	EXPECT_EQ(b2.compute_limbs(), 1);

	BinaryContainer<65> b3;
	EXPECT_EQ(b3.compute_limbs(), 2);

	BinaryContainer<250> b4;
	EXPECT_EQ(b4.compute_limbs(), 4);
}

TEST(Internals, masks){
	BinaryContainer<> b;

	EXPECT_EQ(b.mask(0), 1);
	EXPECT_EQ(b.mask(1), 2);
	EXPECT_EQ(b.mask(2), 4);
}

TEST(Internals, random_limb_with_limits){
	using TestBinaryContainer = BinaryContainer<64>;

	constexpr uint64_t offset = 20;
	using LimbType = TestBinaryContainer::LimbType;
	for (int i = 0; i < 1; ++i) {
		for (int k_lower = 0; k_lower < 64; ++k_lower) {
			for (int k_upper = k_lower+offset; k_upper < 64; ++k_upper) {
				LimbType a = TestBinaryContainer::random_limb(k_lower, k_upper);
				LimbType lmask = TestBinaryContainer::lower_mask(k_lower%64);
				LimbType umask = TestBinaryContainer::higher_mask(k_upper%64);

				EXPECT_NE(a, 0);
				EXPECT_EQ(a&lmask, 0);
				EXPECT_EQ(a&umask, 0);
			}
		}
	}
}



TEST(Zero, Simple) {
	BinaryContainer<> b;
	std::bitset<N> bb;
	b.zero();
	bb.reset();

	for (int i = 0; i < b.size(); ++i) {
		EXPECT_EQ(0, b[i]);
		EXPECT_EQ(bb[i], b[i]);
	}
}

TEST(Zero, Zero_with_Limits) {
	BinaryContainer<> b;

	const int step = 100;
	for (int k_lower = 1; k_lower < b.size(); k_lower += step) {
		for (int k_upper = k_lower+1; k_upper < b.size(); k_upper += step) {
			b.one();
			b.zero(k_lower, k_upper);

			for (int i = 0; i < k_lower; ++i) {
				EXPECT_EQ(1, b[i]);
			}

			for (int i = k_lower; i < k_upper; ++i) {
				EXPECT_EQ(0, b[i]);
			}

			for (int i = k_upper; i < b.size(); ++i) {
				EXPECT_EQ(1, b[i]);
			}
		}
	}

	for (int k_lower = 1; k_lower < b.size(); k_lower += step) {
		for (int k_upper = k_lower+1; k_upper < b.size(); k_upper += step) {
			b.zero();
			b.zero(k_lower, k_upper);

			for (int i = 0; i < k_lower; ++i) {
				EXPECT_EQ(0, b[i]);
			}

			for (int i = k_lower; i < k_upper; ++i) {
				EXPECT_EQ(0, b[i]);
			}

			for (int i = k_upper; i < b.size(); ++i) {
				EXPECT_EQ(0, b[i]);
			}
		}
	}
}

TEST(One, One) {
	BinaryContainer<> b;
	std::bitset<N> bb;
	b.zero();
	b.one();
	bb.reset();

	for (int i = 0; i < b.size(); ++i) {
		EXPECT_NE(bb[i], b[i]);
		EXPECT_EQ(1, b[i]);
	}
}

TEST(One, One_with_Limits) {
	BinaryContainer<> b;

	for (int k_lower = 0; k_lower < b.size(); ++k_lower) {
		for (int k_upper = k_lower+1; k_upper < b.size(); ++k_upper) {
			b.zero();
			b.one(k_lower, k_upper);

			for (int i = 0; i < k_lower; ++i) {
				EXPECT_EQ(0, b[i]);
			}

			for (int i = k_lower; i < k_upper; ++i) {
				EXPECT_EQ(1, b[i]);
			}

			for (int i = k_upper; i < b.size(); ++i) {
				EXPECT_EQ(0, b[i]);
			}
		}
	}

	for (int k_lower = 0; k_lower < b.size(); ++k_lower) {
		for (int k_upper = k_lower+1; k_upper < b.size(); ++k_upper) {
			b.one();
			b.one(k_lower, k_upper);

			for (int i = 0; i < k_lower; ++i) {
				EXPECT_EQ(1, b[i]);
			}

			for (int i = k_lower; i < k_upper; ++i) {
				EXPECT_EQ(1, b[i]);
			}

			for (int i = k_upper; i < b.size(); ++i) {
				EXPECT_EQ(1, b[i]);
			}
		}
	}
}


TEST(Set, Simple) {
	BinaryContainer<> b;
	std::bitset<N> bb;

	bb.reset();
	b.zero();
	b[0] = true;
	EXPECT_EQ(b[0], 1);

	bb[0] = true;
	for (int i = 0; i < b.size(); ++i) {
		// das sind die expliciten casts.
		auto bit = b[i];
		bool bbit = bit;
		EXPECT_EQ(bb[i], bbit);
	}
}

TEST(Set, Random) {
	for (int i = 0; i < TESTSIZE; ++i) {
		BinaryContainer<> b;
		std::bitset<N> bb;
		bb.reset();
		b.zero();

		auto pos = fastrandombytes_uint64() % N;
		b[pos] = bool(fastrandombytes_uint64() % 2);
		bb[pos] = b[pos];
		for (int j = 0; j < b.size(); ++j) {
			EXPECT_EQ(bb[j], b[j]);
		}
	}
}

TEST(Set, Full_Length_Zero) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.one(); b2.zero();

	BinaryContainer<>::set(b1, b2, 0, BinaryContainer<>::ssize());

	for (int j = 0; j < BinaryContainer<>::ssize(); ++j) {
		EXPECT_EQ(0, b1[j]);
	}
}

TEST(Set, Full_Length_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.zero();

	b2[0] = true;
	BinaryContainer<>::set(b1, b2, 0, N);
	EXPECT_EQ(1, b2[0]);

	for (int j = 1; j < b1.size(); ++j) {
		EXPECT_EQ(0, b1[j]);
	}

	// 2. test.
	b1.zero(); b2.one();
	BinaryContainer<>::set(b1, b2, 0, N);
	for (int j = 0; j < b1.size(); ++j) {
		EXPECT_EQ(true, b1[j]);
		EXPECT_EQ(1, b1[j]);

	}
}

TEST(Set, OffByOne_Lower_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.one(); b2.zero();
	BinaryContainer<>::set(b1, b2, 1, N);
	EXPECT_EQ(1, b1[0]);
	for (int j = 1; j < b1.size(); ++j) {
		EXPECT_EQ(0, b1[j]);
	}

	// 2. test.
	b1.zero(); b2.one();
	BinaryContainer<>::set(b1, b2, 1, N);
	EXPECT_EQ(0, b1[0]);
	EXPECT_EQ(false, b1[0]);
	for (int j = 1; j < b1.size(); ++j) {
		EXPECT_EQ(true, b1[j]);
		EXPECT_EQ(1, b1[j]);
	}
}

TEST(Set, OffByOne_Higher_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.zero();

	b1[N-1] = true;   // this should be ignored.
	BinaryContainer<>::set(b1, b2, 0, N - 1);
	EXPECT_EQ(1, b1[N-1]);

	for (int j = 0; j < b1.size() - 1; ++j) {
		EXPECT_EQ(0, b1[j]);
	}

	// 2. test.
	b1.zero(); b2.one();
	BinaryContainer<>::set(b1, b2, 0, N - 1);
	EXPECT_EQ(0, b1[N-1]);
	EXPECT_EQ(false, b1[N-1]);
	for (int j = 0; j < b1.size() - 1; ++j) {
		EXPECT_EQ(true, b1[j]);
		EXPECT_EQ(1, b1[j]);
	}
}

TEST(Set, Complex_Ones) {
	//  test the following:
	//  [111...111]
	// +[000...000] \forall k_lower k_higher
	// =[0...0 1...1 0...0]
	//    k_lower k_higher
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b1.zero(); b2.one();

			BinaryContainer<>::set(b1, b2, k_lower, k_higher);

			for (int j = 0; j < k_lower; ++j) {
				EXPECT_EQ(0, b1[j]);
			}
			for (int j = k_lower; j < k_higher; ++j) {
				EXPECT_EQ(1, b1[j]);
			}
			for (int j = k_higher; j < b1.size(); ++j) {
				EXPECT_EQ(0, b1[j]);
			}
		}
	}
}

TEST(Set, Complex_Zeros) {
	//  test the following:
	//  [111000...000111]
	// +[111000...000111] \forall k_lower k_higher
	// =[1...1 0...0 1...1]
	//    k_lower k_higher

	BinaryContainer<> b1;
	BinaryContainer<> b2;

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b1.zero(); b2.zero();

			for (int j = 0; j < k_lower; ++j) {
				b1[j] = true;
				b2[j] = true;
			}
			for (int j = k_higher; j < b1.size(); ++j) {
				b1[j] = true;
				b2[j] = true;
			}


			BinaryContainer<>::set(b1, b2, k_lower, k_higher);

			for (int j = 0; j < k_lower; ++j) {
				EXPECT_EQ(1, b1[j]);
			}
			for (int j = k_lower; j < k_higher; ++j) {
				EXPECT_EQ(0, b1[j]);
			}
			for (int j = k_higher; j < b1.size(); ++j) {
				EXPECT_EQ(1, b1[j]);
			}
		}
	}
}

TEST(Static_Add, Probabilistic){
	vector<pair<uint64_t, uint64_t>> boundsSet = {pair(0, 64), pair(0, 10), pair(2, 70), pair(64, 128), pair(0, 65), pair(3, 66)};

	for(auto bounds : boundsSet){
		uint64_t k_lower = bounds.first;
		uint64_t k_upper = bounds.second;

		for(uint64_t i = 0; i < 100; i++){
			uint64_t a = fastrandombytes_uint64();
			uint64_t b = fastrandombytes_uint64();
			uint64_t c = fastrandombytes_uint64();
			uint64_t d = fastrandombytes_uint64();
			uint64_t e = fastrandombytes_uint64();
			uint64_t f = fastrandombytes_uint64();

			BinaryContainer<128> b1;
			BinaryContainer<128> b2;
			BinaryContainer<128> b3;

			b1.data()[0] = a; b1.data()[1] = b;
			b2.data()[0] = c; b2.data()[1] = d;
			b3.data()[0] = e; b3.data()[1] = f;


			BinaryContainer<128>::add(b3, b2, b1, k_lower, k_upper);

			for(uint64_t k = 0; k < k_lower; k++){
				if(k < 64){
					ASSERT_EQ(b3.get_bit_shifted(k), (e>>k) & 1);
				}
				else {
					ASSERT_EQ(b3.get_bit_shifted(k), (f>>k) & 1);
				}
			}
			for(uint64_t k = k_lower; k < k_upper; k++){
				if(k < 64){
					ASSERT_EQ(b3.get_bit_shifted(k), ((a^c) >> k) & 1);
				}
				else {
					ASSERT_EQ(b3.get_bit_shifted(k), ((b^d) >> k) & 1);
				}
			}
			for(uint64_t k = k_upper; k < 128; k++){
				if(k < 64){
					ASSERT_EQ(b3.get_bit_shifted(k), (e>>k) & 1);
				}
				else {
					ASSERT_EQ(b3.get_bit_shifted(k), (f>>(k-64)) & 1);
				}
			}
		}
	}

}


TEST(Add, Probabilistic){
	using BinaryContainerTest = BinaryContainer<128>;
	vector<pair<uint64_t, uint64_t>> boundsSet = {pair(0, 64), pair(0, 10), pair(2, 70), pair(64, 128), pair(0, 65), pair(3, 66)};

	for(auto bounds : boundsSet) {
		uint64_t k_lower = bounds.first;
		uint64_t k_upper = bounds.second;

		for(uint64_t i = 0; i < 100; i++) {
			uint64_t a = fastrandombytes_uint64();
			uint64_t b = fastrandombytes_uint64();
			uint64_t c = fastrandombytes_uint64();
			uint64_t d = fastrandombytes_uint64();

			BinaryContainerTest b1;
			BinaryContainerTest b2;

			b1.data()[0] = a; b1.data()[1] = b;
			b2.data()[0] = c; b2.data()[1] = d;


			b1.add(b2, k_lower, k_upper);

			for(uint64_t k = 0; k < k_lower; k++){
				if(k < 64){
					ASSERT_EQ(b1.get_bit_shifted(k), (a>>k) & 1);
				}
				else {
					ASSERT_EQ(b1.get_bit_shifted(k), (b>>k) & 1);
				}
			}
			for(uint64_t k = k_lower; k < k_upper; k++){
				if(k < 64){
					ASSERT_EQ(b1.get_bit_shifted(k), ((a^c) >> k) & 1);
				}
				else {
					ASSERT_EQ(b1.get_bit_shifted(k), ((b^d) >> k) & 1);
				}
			}
			for(uint64_t k = k_upper; k < 128; k++){
				if(k < 64){
					ASSERT_EQ(b1.get_bit_shifted(k), (a>>k) & 1);
				}
				else {
					ASSERT_EQ(b1.get_bit_shifted(k), (b>>k) & 1);
				}
			}
		}
	}
}

TEST(Add, Norm){
	using BinaryContainerTest = BinaryContainer<128>;

	BinaryContainerTest b1;
	BinaryContainerTest b2;
	BinaryContainerTest b3;

	std::vector<std::pair<uint64_t, uint64_t>> boundsSet = {pair(0, 128)};

	for(auto bounds : boundsSet){
		uint64_t k_lower = bounds.first;
		uint64_t k_upper = bounds.second;
		b1.one();
		b2.zero();


		uint64_t norm = (k_lower+k_upper) / 2;

		bool result = BinaryContainerTest::add(b3, b2, b1, k_lower, k_upper, norm);
		ASSERT_EQ(true, result);
	}
}

TEST(Add, Full_Length_Zero) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	BinaryContainer<>::add(b3, b1, b2, 0, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Add, Full_Length_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	b1[0] = true;
	BinaryContainer<>::add(b3, b1, b2, 0, N);
	EXPECT_EQ(1, b3[0]);

	for (int j = 1; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}

	// 2. test.
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
	}

	BinaryContainer<>::add(b3, b1, b2, 0, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(true, b3[j]);
		EXPECT_EQ(1, b3[j]);

	}

	//3.test
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
		b2[i] = true;
	}

	BinaryContainer<>::add(b3, b1, b2, 0, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(false, b3[j]);
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Add, OffByOne_Lower_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	b1[0] = true;   // this should be ignored.
	BinaryContainer<>::add(b3, b1, b2, 1, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}

	// 2. test.
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
	}

	BinaryContainer<>::add(b3, b1, b2, 1, N);
	EXPECT_EQ(0, b3[0]);
	EXPECT_EQ(false, b3[0]);
	for (int j = 1; j < b3.size(); ++j) {
		EXPECT_EQ(true, b3[j]);
		EXPECT_EQ(1, b3[j]);
	}

	//3.test
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
		b2[i] = true;
	}

	BinaryContainer<>::add(b3, b1, b2, 1, N);
	EXPECT_EQ(0, b3[0]);
	EXPECT_EQ(false, b3[0]);

	for (int j = 1; j < b3.size(); ++j) {
		EXPECT_EQ(false, b3[j]);
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Add, OffByOne_Higher_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	b1[N-1] = true;   // this should be ignored.
	BinaryContainer<>::add(b3, b1, b2, 0, N - 1);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}

	// 2. test.
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
	}

	BinaryContainer<>::add(b3, b1, b2, 0, N - 1);
	EXPECT_EQ(0, b3[N-1]);
	EXPECT_EQ(false, b3[N-1]);
	for (int j = 0; j < b3.size() - 1; ++j) {
		EXPECT_EQ(true, b3[j]);
		EXPECT_EQ(1, b3[j]);
	}

	//3.test
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
		b2[i] = true;
	}

	BinaryContainer<>::add(b3, b1, b2, 0, N - 1);
	EXPECT_EQ(0, b3[N-1]);
	EXPECT_EQ(false, b3[N-1]);

	for (int j = 1; j < b3.size() - 1; ++j) {
		EXPECT_EQ(false, b3[j]);
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Add, Complex_Ones) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b1.zero(); b2.zero(); b3.zero();

			for (int i = 0; i < b1.size(); ++i) {
				b1[i] = true;
			}

			BinaryContainer<>::add(b3, b1, b2, k_lower, k_higher);

			for (int j = 0; j < k_lower; ++j) {
				EXPECT_EQ(0, b3[j]);
			}
			for (int j = k_lower; j < k_higher; ++j) {
				EXPECT_EQ(1, b3[j]);
			}
			for (int j = k_higher; j < b1.size(); ++j) {
				EXPECT_EQ(0, b3[j]);
			}
		}
	}
}

TEST(Sub, Full_Length_Zero) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	BinaryContainer<>::sub(b3, b1, b2, 0, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Sub, Full_Length_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	b1[0] = true;
	BinaryContainer<>::sub(b3, b1, b2, 0, N);
	EXPECT_EQ(1, b3[0]);

	for (int j = 1; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}

	// 2. test.
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
	}

	BinaryContainer<>::sub(b3, b1, b2, 0, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(true, b3[j]);
		EXPECT_EQ(1, b3[j]);

	}

	//3.test
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
		b2[i] = true;
	}

	BinaryContainer<>::sub(b3, b1, b2, 0, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(false, b3[j]);
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Sub, OffByOne_Lower_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	b1[0] = true;   // this should be ignored.
	BinaryContainer<>::sub(b3, b1, b2, 1, N);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}

	// 2. test.
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
	}

	BinaryContainer<>::sub(b3, b1, b2, 1, N);
	EXPECT_EQ(0, b3[0]);
	EXPECT_EQ(false, b3[0]);
	for (int j = 1; j < b3.size(); ++j) {
		EXPECT_EQ(true, b3[j]);
		EXPECT_EQ(1, b3[j]);
	}

	//3.test
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
		b2[i] = true;
	}

	BinaryContainer<>::sub(b3, b1, b2, 1, N);
	EXPECT_EQ(0, b3[0]);
	EXPECT_EQ(false, b3[0]);

	for (int j = 1; j < b3.size(); ++j) {
		EXPECT_EQ(false, b3[j]);
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Sub, OffByOne_Higher_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	b1.zero(); b2.zero(); b3.zero();

	b1[N-1] = true;   // this should be ignored.
	BinaryContainer<>::sub(b3, b1, b2, 0, N - 1);
	for (int j = 0; j < b3.size(); ++j) {
		EXPECT_EQ(0, b3[j]);
	}

	// 2. test.
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
	}

	BinaryContainer<>::sub(b3, b1, b2, 0, N - 1);
	EXPECT_EQ(0, b3[N-1]);
	EXPECT_EQ(false, b3[N-1]);
	for (int j = 0; j < b3.size() - 1; ++j) {
		EXPECT_EQ(true, b3[j]);
		EXPECT_EQ(1, b3[j]);
	}

	//3.test
	b1.zero(); b2.zero(); b3.zero();
	for (int i = 0; i < b1.size(); ++i) {
		b1[i] = true;
		b2[i] = true;
	}

	BinaryContainer<>::sub(b3, b1, b2, 0, N - 1);
	EXPECT_EQ(0, b3[N-1]);
	EXPECT_EQ(false, b3[N-1]);

	for (int j = 1; j < b3.size() - 1; ++j) {
		EXPECT_EQ(false, b3[j]);
		EXPECT_EQ(0, b3[j]);
	}
}

TEST(Sub, Complex_Ones) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;
	BinaryContainer<> b3;

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b1.zero(); b2.zero(); b3.zero();

			for (int i = 0; i < b1.size(); ++i) {
				b1[i] = true;
			}

			BinaryContainer<>::sub(b3, b1, b2, k_lower, k_higher);

			for (int j = 0; j < k_lower; ++j) {
				EXPECT_EQ(0, b3[j]);
			}
			for (int j = k_lower; j < k_higher; ++j) {
				EXPECT_EQ(1, b3[j]);
			}
			for (int j = k_higher; j < b1.size(); ++j) {
				EXPECT_EQ(0, b3[j]);
			}
		}
	}
}


TEST(Cmp, Simple_Everything_False) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.one();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));
		}
	}
}

TEST(Cmp, Simple_Everything_True) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));
		}
	}
}

TEST(Cmp, OffByOne_Lower_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.zero();

	b1[0] = true;
	EXPECT_EQ(1, b1[0]);
	EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, 0, b1.size()));
	for (int j = 1; j < b1.size(); ++j) {
		EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, j, b1.size()));
	}
}

TEST(Cmp, OffByOne_Higher_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.zero();

	b1[N-1] = true;
	EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, 0, N));
	for (int j = 0; j < b1.size() - 1; ++j) {
		EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, j, N - 1));
	}
}

TEST(Cmp, Complex_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b2.zero();
			for (int i = k_lower; i < k_higher; ++i) {
				b2[i] = true;
			}

			if (k_lower > 0)
				EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, 0, k_lower));
			EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));
			EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, k_higher, b1.size()));


			b2.zero();
			EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));

			for (int i = k_higher; i < b2.size(); ++i) {
				b2[i] = true;
			}
			if (k_lower > 0)
				EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, 0, k_lower));
			EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));
			EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, k_higher, b1.size()));
		}
	}
}

TEST(Cmp, Complex_Zero) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero();
	b2.zero();

	for (int i = 0; i < b1.size(); ++i) { b1[i] = true; }

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b2.zero();
			for (int i = k_lower; i < k_higher; ++i) {
				b2[i] = true;
			}
			if (k_lower > 0)
				EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, 0, k_lower));

			EXPECT_EQ(true, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));
			EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, k_higher, b1.size()));


			b2.zero();
			EXPECT_EQ(false, BinaryContainer<>::cmp(b1, b2, k_lower, k_higher));
		}
	}
}

TEST(Cmp, Special_OffByOne_Lower_Zero) {
	constexpr uint64_t size = 64;
	using TestBinaryContainer = BinaryContainer<size>;
	TestBinaryContainer b1;
	TestBinaryContainer b2;

	b2.zero();

	for (int i = 0; i < size-2; ++i) {
		b1.zero();
		b1[i] = true;
		EXPECT_EQ(1, b1[i]);
		EXPECT_EQ(false, TestBinaryContainer::cmp(b1, b2, 0, size));

		if (i > 0)
			EXPECT_EQ(false, TestBinaryContainer::cmp(b1, b2, 0, i+1));

		for (int j = i+1; j < size; ++j) {
			EXPECT_EQ(true, TestBinaryContainer::cmp(b1, b2, j, size));
		}
	}

}

TEST(Cmp, Special_OffByOne_Lower_One) {
	constexpr uint64_t size = 64;
	using TestBinaryContainer = BinaryContainer<size>;
	TestBinaryContainer b1;
	TestBinaryContainer b2;

	b2.one();

	for (int i = 0; i < size-2; ++i) {
		b1.zero();
		b1[i] = true;
		EXPECT_EQ(1, b1[i]);
		EXPECT_EQ(false, TestBinaryContainer::cmp(b1, b2, 0, size));

		if (i > 0)
			EXPECT_EQ(false, TestBinaryContainer::cmp(b1, b2, 0, i+1));

		for (int j = i+1; j < size; ++j) {
			EXPECT_EQ(false, TestBinaryContainer::cmp(b1, b2, j, size));
		}

		b1.one();
		for (int j = 0; j < size; ++j) {
			EXPECT_EQ(true, TestBinaryContainer::cmp(b1, b2, j, size));
		}
	}

}

TEST(IsGreater, Simple_Everything_False) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.one();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			EXPECT_EQ(false, b1.is_greater(b2, k_lower, k_higher));
			EXPECT_EQ(true, b2.is_greater(b1, k_lower, k_higher));
		}
	}
}

TEST(IsGreater, Simple_Everything_True) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.one(); b2.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			// std ::cout << k_lower << " " << k_higher << "\n";
			auto b = b1.is_greater(b2, k_lower, k_higher);
			EXPECT_EQ(true, b);
		}
	}
}

TEST(IsGreater, Complex_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b2.zero();
			for (int i = 0; i < k_higher; ++i) {
				b2[i] = 1;
			}

			EXPECT_EQ(false, b1.is_greater(b2, k_lower, k_higher));
			EXPECT_EQ(false, b1.is_greater(b2, k_higher, b1.size()));
			EXPECT_EQ(true, b2.is_greater(b1, k_lower, k_higher));

			if (k_higher < b1.size() - 1){
				b1[k_higher] = 1;
				EXPECT_EQ(true, b1.is_greater(b2, k_higher, b1.size()));
				b1.zero();
			}


			b2.zero();
			for (int i = k_higher; i < b2.size(); ++i) {
				b2[i] = 1;
			}

			EXPECT_EQ(false, b1.is_greater(b2, k_lower, k_higher));
			EXPECT_EQ(false, b1.is_greater(b2, k_higher, b1.size()));

			EXPECT_EQ(false, b2.is_greater(b1, k_lower, k_higher));
			EXPECT_EQ(true, b2.is_greater(b1, k_higher, b1.size()));

			if (k_higher < b1.size() - 1){
				EXPECT_EQ(true, b2.is_greater(b1, k_lower, k_higher+1));
				EXPECT_EQ(false, b1.is_greater(b2, k_lower, k_higher+1));

			}
		}
	}
}

TEST(IsGreater, Complex) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b2.zero();
			for (int i = 0; i < k_higher; ++i) {
				b2[i] = true;
			}

			EXPECT_EQ(false, b1.is_greater(b2, k_lower, k_higher));
			EXPECT_EQ(false, b1.is_greater(b2, k_higher, b1.size()));
			EXPECT_EQ(true, b2.is_greater(b1, k_lower, k_higher));

			if (k_higher < b1.size() - 1){
				b1[k_higher] = 1;
				EXPECT_EQ(true, b1.is_greater(b2, k_higher, b1.size()));
				b1.zero();
			}
		}
	}
}


TEST(IsLower, Simple_Everything_False) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero(); b2.one();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			EXPECT_EQ(false, b2.is_lower(b1, k_lower, k_higher));
			EXPECT_EQ(true, b1.is_lower(b2, k_lower, k_higher));
		}
	}
}

TEST(IsLower, Simple_Everything_True) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.one(); b2.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			EXPECT_EQ(false, b1.is_lower(b2, k_lower, k_higher));
			EXPECT_EQ(true, b2.is_lower(b1, k_lower, k_higher));

		}
	}
}

TEST(IsLower, Complex_One) {
	BinaryContainer<> b1;
	BinaryContainer<> b2;

	b1.zero();

	for (int k_lower  = 0; k_lower < b1.size(); ++k_lower) {
		for (int k_higher = k_lower + 1; k_higher < b1.size(); ++k_higher) {
			b2.zero();
			for (int i = 0; i < k_higher; ++i) {
				b2[i] = 1;
			}

			EXPECT_EQ(true, b1.is_lower(b2, k_lower, k_higher));
			EXPECT_EQ(false, b2.is_lower(b1, k_lower, k_higher));
			EXPECT_EQ(false, b1.is_lower(b2, k_higher, b1.size()));

			if (k_higher < b1.size() - 1){
				b1[k_higher] = 1;
				EXPECT_EQ(false, b1.is_lower(b2, k_higher, b1.size()));
				EXPECT_EQ(true, b2.is_lower(b1, k_higher, b1.size()));

				b1.zero();
			}


			b2.zero();
			for (int i = k_higher; i < b2.size(); ++i) {
				b2[i] = 1;
			}

			EXPECT_EQ(false, b1.is_lower(b2, k_lower, k_higher));
			EXPECT_EQ(true, b1.is_lower(b2, k_higher, b1.size()));
			EXPECT_EQ(false, b2.is_lower(b1, k_higher, b1.size()));

			if (k_higher < b1.size() - 1){
				EXPECT_EQ(false, b2.is_lower(b1, k_lower, k_higher+1));
				EXPECT_EQ(true, b1.is_lower(b2, k_lower, k_higher+1));
			}
		}
	}
}


TEST(weight, Simple_Everything_True) {
	BinaryContainer<> b1;
	b1.zero();

	for (int k_lower  = 1; k_lower < b1.size(); ++k_lower) {
		b1[k_lower-1] = true;
		EXPECT_EQ(k_lower, b1.weight());
		EXPECT_EQ(k_lower, b1.weight(0, k_lower));
		EXPECT_EQ(k_lower, b1.weight(0, b1.size()));
		if (k_lower + 1 < b1.size())
			EXPECT_EQ(0, b1.weight(k_lower +1, b1.size()));
	}
}

int main(int argc, char **argv) {
	InitGoogleTest(&argc, argv);
	random_seed(time(NULL));
	return RUN_ALL_TESTS();
}
