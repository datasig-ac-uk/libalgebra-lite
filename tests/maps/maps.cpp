//
// Created by user on 09/09/22.
//


#include "maps_fixture.h"

TEST_F(MapsFixture, testExpandLetter) {
    const auto& result = maps.expand(lkey(1));
    ASSERT_EQ(1, result.size());
    EXPECT_EQ(tkey(1), result[0].first);
    EXPECT_EQ(1, result[0].second);
}

TEST_F(MapsFixture, testExpandWord2) {
    const auto& result = maps.expand(lkey(1, 2));

    ASSERT_EQ(2, result.size());
    EXPECT_EQ(tkey(1, 2), result[0].first);
    EXPECT_EQ(1, result[0].second);

    EXPECT_EQ(tkey(2, 1), result[1].first);
    EXPECT_EQ(-1, result[1].second);
}

TEST_F(MapsFixture, testExpandWord3) {
    const auto& result = maps.expand(lkey(1, lkey(1, 2)));

    // [1, [1, 2]] = 1 x (1 x 2 - 2 x 1) - (1 x 2 - 2 x 1) x 1
    //             = 1 x 1 x 2 - 1 x 2 x 1 - 1 x 2 x 1 + 2 x 1 x 1
    ASSERT_EQ(3, result.size());

    EXPECT_EQ(tkey(1, 1, 2), result[0].first);
    EXPECT_EQ(1, result[0].second);

    EXPECT_EQ(tkey(1, 2, 1), result[1].first);
    EXPECT_EQ(-2, result[1].second);

    EXPECT_EQ(tkey(2, 1, 1), result[2].first);
    EXPECT_EQ(1, result[2].second);
}



TEST_F(MapsFixture, testRBracketingLetter) {
    const auto& result = maps.rbracketing(tkey(1));

    ASSERT_EQ(1, result.size());

    EXPECT_EQ(lkey(1), result[0].first);
    EXPECT_EQ(1, result[0].second);
}

TEST_F(MapsFixture, testRBracketingWord2Ordered) {
    const auto& result = maps.rbracketing(tkey(1, 2));

    ASSERT_EQ(1, result.size());

    EXPECT_EQ(lkey(1, 2), result[0].first);
    EXPECT_EQ(1, result[0].second);
}

TEST_F(MapsFixture, testRBracketingWord2UnOrdered) {
    const auto& result = maps.rbracketing(tkey(2, 1));

    ASSERT_EQ(1, result.size());

    EXPECT_EQ(lkey(1, 2), result[0].first);
    EXPECT_EQ(-1, result[0].second);
}

TEST_F(MapsFixture, testRBracketingWord2SameLetter) {
    const auto& result = maps.rbracketing(tkey(1, 1));

    ASSERT_EQ(0, result.size());
}

TEST_F(MapsFixture, testRBracketingWord3Ordered) {
    const auto& result = maps.rbracketing(tkey(1, 1, 2));

    // [1, [1, 2]]
    ASSERT_EQ(1, result.size());

    EXPECT_EQ(lkey(1, lkey(1, 2)), result[0].first);
    EXPECT_EQ(1, result[0].second);
}

TEST_F(MapsFixture, testRBracketingWord3UnOrdered) {
    const auto& result = maps.rbracketing(tkey(1, 2, 1));

    // [1, [2, 1]] = -[1, [1, 2]]
    ASSERT_EQ(1, result.size());

    EXPECT_EQ(lkey(1, lkey(1, 2)), result[0].first);
    EXPECT_EQ(-1, result[0].second);
}

TEST_F(MapsFixture, testRBracketingWord3RepeatedLetter) {
    const auto& result = maps.rbracketing(tkey(2, 1, 1));

    ASSERT_EQ(0, result.size());
}


TEST_F(MapsFixture, testLieToTensor) {

    const auto test_lie = generic_lie('x');

    const auto test_tensor = maps.lie_to_tensor(test_lie);
    const auto result_lie = maps.tensor_to_lie(test_tensor);

    ASSERT_EQ(result_lie, test_lie);
}

// This test doesn't pass because test_tensor is not log of a group-like element.
//TEST_F(MapsFixture, testTensorToLie) {
//
//    const auto test_tensor = generic_dtensor('x');
//    const auto test_lie = maps.tensor_to_lie(test_tensor);
//    const auto result_tensor = maps.lie_to_tensor(test_lie);
//
//    ASSERT_EQ(result_tensor, test_tensor) << (result_tensor - test_tensor);
//}
