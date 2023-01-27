//
// Created by user on 08/09/22.
//

#include "lie_fixture.h"

#include <sstream>

TEST_F(LieFixture, testDegreeLetter) {
    auto k = key(1);

    EXPECT_EQ(1, basis->degree(k));
}

TEST_F(LieFixture, testDegreeKey) {
    auto k = key(1, 2);
    EXPECT_EQ(2, basis->degree(k));
}

TEST_F(LieFixture, testDegreeCompoundKey) {
    auto k = key(key(1), key(1, 2));
    EXPECT_EQ(3, basis->degree(k));
}


TEST_F(LieFixture, testParentsLetter) {
    auto k = key(1);

    ASSERT_EQ(k, key_type(1, 0));

    auto parents = basis->parents(k);

    EXPECT_EQ(parents.first, key_type(0, 0));
    EXPECT_EQ(parents.second, k);

}

TEST_F(LieFixture, testParentsKey) {
    auto k = key(1, 2);

    ASSERT_EQ(k, key_type(2, 0));

    auto parents = basis->parents(k);

    EXPECT_EQ(key(1), parents.first);
    EXPECT_EQ(key(2), parents.second);
}

TEST_F(LieFixture, testParentsCompoundKey) {
    auto k = key(key(1, 2), key(3, 4));

    auto parents = basis->parents(k);

    EXPECT_EQ(key(1, 2), parents.first);
    EXPECT_EQ(key(3, 4), parents.second);
}



TEST_F(LieFixture, testPrintLetter) {
    auto k = key(1);

    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("1", ss.str());
}

TEST_F(LieFixture, testPrintBracketLetters) {
    auto k = key(1, 2);

    ASSERT_EQ(k, key_type(2, 0));
    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("[1,2]", ss.str());
}

TEST_F(LieFixture, testPrintBracketKeys) {
    auto k = key(key(1, 2), key(3, 4));

    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("[[1,2],[3,4]]", ss.str());
}


TEST_F(LieFixture, testStartOfDegree) {
    EXPECT_EQ(0, basis->start_of_degree(0));
    EXPECT_EQ(0, basis->start_of_degree(1));
    EXPECT_EQ(5, basis->start_of_degree(2));
    EXPECT_EQ(5 + 10, basis->start_of_degree(3));
    EXPECT_EQ(5 + 10 + 40, basis->start_of_degree(4));
    EXPECT_EQ(5 + 10 + 40 + 150, basis->start_of_degree(5));
    EXPECT_EQ(5 + 10 + 40 + 150 + 624, basis->start_of_degree(6));
}

TEST_F(LieFixture, testStartOfDegreeSize) {
    for (deg_t i=0; i<=5; ++i) {
        EXPECT_EQ(basis->size(i), basis->start_of_degree(i+1));
    }
}

TEST_F(LieFixture, testSizeOfDegreeVsSize) {
    for (deg_t i=0; i<5; ++i) {
        EXPECT_EQ(basis->size(i+1) - basis->size(i), basis->size_of_degree(i+1));
    }
}


TEST_F(LieFixture, testIndexToKeyRoundtrip) {

    for (lal::dimn_t i=0; i<basis->size(-1); ++i) {
        auto k = basis->index_to_key(i);
        auto index = basis->key_to_index(k);

        EXPECT_EQ(i, index);
    }


}
