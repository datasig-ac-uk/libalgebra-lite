//
// Created by user on 07/02/23.
//


#include "tensor_fixture.h"

namespace {

struct SparseTensorTests : TensorFixture {

    using tensor_type = lal::free_tensor<lal::polynomial_ring, lal::sparse_vector, lal::dtl::standard_storage>;


    tensor_type generic_tensor(char letter) const
    { return TensorFixture::template generic_tensor<tensor_type>(letter); }

};


}


TEST_F(SparseTensorTests, testUminus) {
    auto lhs = this->generic_tensor('x');

    auto result = -lhs;

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = -scalar_type(this->key_to_monomial('x', k), 1);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(result, expected);
}

TEST_F(SparseTensorTests, testAddition) {
    auto lhs = this->generic_tensor('x');
    auto rhs = this->generic_tensor('y');

    auto result = lhs + rhs;

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) + scalar_type(this->key_to_monomial('y', k), 1);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(result, expected);
}

TEST_F(SparseTensorTests, testSubtraction) {
    auto lhs = this->generic_tensor('x');
    auto rhs = this->generic_tensor('y');

    auto result = lhs - rhs;

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) - scalar_type(this->key_to_monomial('y', k), 1);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(result, expected);
}

TEST_F(SparseTensorTests, testScalarMultiplyLeft) {
    auto lhs = this->generic_tensor('x');
    scalar_type scalar(lal::monomial(typename lal::monomial::letter_type('z', 0)), 1);

    auto result = scalar*lhs;

    tensor_type expected(basis);
    this->for_each_key([this, &expected, &scalar](key_type k) {
        auto tmp = scalar*scalar_type(this->key_to_monomial('x', k), 1);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(result, expected);
}

TEST_F(SparseTensorTests, testScalarMultiplyRight) {
    auto lhs = this->generic_tensor('x');
    scalar_type scalar(lal::monomial(typename lal::monomial::letter_type('z', 0)), 1);

    auto result = lhs*scalar;

    tensor_type expected(basis);
    this->for_each_key([this, &expected, &scalar](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1)*scalar;
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(result, expected);
}

TEST_F(SparseTensorTests, testRationalDivide) {
    auto lhs = this->generic_tensor('x');

    auto result = lhs / rational_type(2);

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) / rational_type(2);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(result, expected);
}


TEST_F(SparseTensorTests, testInplaceAddition) {
    auto lhs = this->generic_tensor('x');
    auto rhs = this->generic_tensor('y');

    lhs += rhs;

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) + scalar_type(this->key_to_monomial('y', k), 1);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(lhs, expected);
}

TEST_F(SparseTensorTests, testInplaceSubtraction) {
    auto lhs = this->generic_tensor('x');
    auto rhs = this->generic_tensor('y');

    lhs -= rhs;

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) - scalar_type(this->key_to_monomial('y', k), 1);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(lhs, expected);
}

TEST_F(SparseTensorTests, testInplaceScalarMultiplyRight) {
    auto lhs = this->generic_tensor('x');
    scalar_type scalar(lal::monomial(typename lal::monomial::letter_type('z', 0)), 1);

    lhs *= scalar;

    tensor_type expected(basis);
    this->for_each_key([this, &expected, &scalar](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) * scalar;
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(lhs, expected);
}

TEST_F(SparseTensorTests, testInplaceRationalDivide) {
    auto lhs = this->generic_tensor('x');

    lhs /= rational_type(2);

    tensor_type expected(basis);
    this->for_each_key([this, &expected](key_type k) {
        auto tmp = scalar_type(this->key_to_monomial('x', k), 1) / rational_type(2);
        expected.add_scal_prod(k, tmp);
    });

    ASSERT_EQ(lhs, expected) << lhs  << "\n\n" << expected;
}


TEST_F(SparseTensorTests, testAddScalProdKeyScalar) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});
    scalar_type polykey(this->key_to_monomial('z', key), 1);

    lhs.add_scal_prod(key, polykey);

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect += polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}

TEST_F(SparseTensorTests, testAddScalProdKeyInto) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});

    lhs.add_scal_prod(key, 1);
    scalar_type polykey(1);

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect += polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}


TEST_F(SparseTensorTests, testSubScalProdKeyScalar) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});
    scalar_type polykey(this->key_to_monomial('z', key), 1);

    lhs.sub_scal_prod(key, polykey);

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect -= polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}

TEST_F(SparseTensorTests, testSubScalProdKeyInto) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});

    lhs.sub_scal_prod(key, 1);
    scalar_type polykey(1);

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect -= polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}



TEST_F(SparseTensorTests, testAddScalDivKeyScalar) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});
    scalar_type polykey(1);
    polykey /= 2;

    lhs.add_scal_div(key, rational_type(2));

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect += polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}

TEST_F(SparseTensorTests, testAddScalDivKeyInto) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});
    scalar_type polykey(1);
    polykey /= 2;

    lhs.add_scal_div(key, 2);

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect += polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}


TEST_F(SparseTensorTests, testSubScalDivKeyScalar) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});
    scalar_type polykey(1);
    polykey /= 2;


    lhs.sub_scal_div(key, rational_type(2));

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect -= polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}

TEST_F(SparseTensorTests, testSubScalDivKeyInto) {
    auto lhs = this->generic_tensor('x');

    auto key = this->key({1, 2, 3, 4, 5});
    scalar_type polykey(1);
    polykey /= 2;

    lhs.sub_scal_div(key, 2);

    this->for_each_key([this, &key, &polykey, &lhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        if (k == key) {
            expect -= polykey;
        }
        EXPECT_EQ(expect, lhs[k]) << expect << "\n\n" << lhs[k];
    });

}


TEST_F(SparseTensorTests, testSubScalProdVectorScalar) {
    auto lhs = this->generic_tensor('x');
    const auto rhs = this->generic_tensor('y');
    scalar_type poly(lal::monomial(typename lal::monomial::letter_type('z', 0)), 1);

    lhs.sub_scal_prod(rhs, poly);

    this->for_each_key([this, &lhs, &rhs, &poly](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        expect -= (rhs[k]*poly);
        EXPECT_EQ(expect, lhs[k]);
    });
}


TEST_F(SparseTensorTests, testAddScalDivVectorScalar) {
    auto lhs = this->generic_tensor('x');
    const auto rhs = this->generic_tensor('y');

    lhs.add_scal_div(rhs, 2);

    this->for_each_key([this, &lhs, &rhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        expect += (rhs[k]/2);
        EXPECT_EQ(expect, lhs[k]);
    });
}

TEST_F(SparseTensorTests, testSubScalDivVectorScalar) {
    auto lhs = this->generic_tensor('x');
    const auto rhs = this->generic_tensor('y');

    lhs.sub_scal_div(rhs, 2);

    this->for_each_key([this, &lhs, &rhs](key_type k) {
        scalar_type expect(this->key_to_monomial('x', k), 1);
        expect -= (rhs[k]/2);
        EXPECT_EQ(expect, lhs[k]);
    });
}
