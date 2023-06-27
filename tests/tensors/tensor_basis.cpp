//
// Created by user on 08/09/22.
//

#include "tensor_fixture.h"
#include <sstream>

TEST_F(TensorFixture, testDegreeEmpty) {
    auto k = key();
    EXPECT_EQ(0, basis->degree(k));
}

TEST_F(TensorFixture, testDegreeLetter) {
    EXPECT_EQ(1, basis->degree(key(1)));
}

TEST_F(TensorFixture, testDegreeWord) {
    EXPECT_EQ(3, basis->degree(key(1, 2, 3)));
}

TEST_F(TensorFixture, testParentsEmpty) {
    auto k = key();
    EXPECT_EQ(key(), basis->lparent(k));
    EXPECT_EQ(key(), basis->rparent(k));
}

TEST_F(TensorFixture, testParentsLetter) {
    auto k = key(1);
    EXPECT_EQ(key(1), basis->lparent(k));
    EXPECT_EQ(key(), basis->rparent(k));
}

TEST_F(TensorFixture, testParentsWord2) {
    auto k = key(1, 2);
    EXPECT_EQ(key(1), basis->lparent(k));
    EXPECT_EQ(key(2), basis->rparent(k));
}

TEST_F(TensorFixture, testParentsWord3) {
    auto k = key(1, 2, 3);
    EXPECT_EQ(key(1), basis->lparent(k));
    EXPECT_EQ(key(2, 3), basis->rparent(k));
}

TEST_F(TensorFixture, testPrintEmptyKey) {
    auto k = key();

    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("", ss.str());
}

TEST_F(TensorFixture, testPrintLetter) {
    auto k = key(1);

    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("1", ss.str());
}

TEST_F(TensorFixture, testPrintLength2) {
    auto k = key(1, 2);

    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("1,2", ss.str());
}

TEST_F(TensorFixture, testPrintLongKey) {
    auto k = key(1, 2, 3, 4, 5);

    std::stringstream ss;
    basis->print_key(ss, k);

    EXPECT_EQ("1,2,3,4,5", ss.str());
}


TEST_F(TensorFixture, testIndexToKeyRoundtrip) {

    for (lal::dimn_t i=0; i<basis->size(-1); ++i) {
        auto k = basis->index_to_key(i);
        EXPECT_EQ(i, basis->key_to_index(k));
    }

}


TEST_F(TensorFixture, testStartOfDegree) {
    auto size = 0;
    for (auto i=0; i<=depth; ++i) {
        EXPECT_EQ(size, basis->start_of_degree(i));
        size *= width;
        size += 1;
    }
}

TEST_F(TensorFixture, testSize) {
    auto sz = 0;
    for (auto i=0; i<=depth; ++i) {
        sz *= width;
        sz += 1;
        EXPECT_EQ(sz, basis->size(i));
    }
}

TEST_F(TensorFixture, testSizeOfDegreeVsSize) {
    for (auto i=0; i<=depth; ++i) {
        EXPECT_EQ(basis->size(i) - basis->start_of_degree(i),
                basis->size_of_degree(i));
    }
}

TEST_F(TensorFixture, testLetterToKey) {
    for (lal::let_t l=1; l<=static_cast<lal::let_t>(width); ++l) {
        EXPECT_EQ(key(l), basis->key_of_letter(l));
    }
}
