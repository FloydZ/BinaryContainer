#include <gtest/gtest.h>
#include <cstdint>
#include <bitset>
#include <utility>

#include "random.h"
#include "binary_container.h"
#include "example.h"

using ::testing::EmptyTestEventListener;
using ::testing::InitGoogleTest;
using ::testing::Test;
using ::testing::TestEventListeners;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::UnitTest;

TEST(Internals, Alignment) {

}

int main(int argc, char **argv) {
	InitGoogleTest(&argc, argv);
	random_seed(time(NULL));
	return RUN_ALL_TESTS();
}
