#include "vec2.hh"
#include <gtest/gtest.h>

TEST(Vec2Test, AssignOpsVec) {
    Vec2 a = {24, 195};
    Vec2 b = {6, 5};
    a /= b;
    ASSERT_EQ(a, Vec2(4, 39));
    a -= b;
    ASSERT_EQ(a, Vec2(-2, 34));
    a += b;
    ASSERT_EQ(a, Vec2(4, 39));
    a *= b;
    ASSERT_EQ(a, Vec2(24, 195));
}

TEST(Vec2Test, AssignOpsFloat) {
    Vec2 a = {85, 289};
    float b = 17;
    a /= b;
    ASSERT_EQ(a, Vec2(5, 17));
    a -= b;
    ASSERT_EQ(a, Vec2(-12, 0));
    a += b;
    ASSERT_EQ(a, Vec2(5, 17));
    a *= b;
    ASSERT_EQ(a, Vec2(85, 289));
}

TEST(Vec2Test, ScalarProduct) {
    Vec2 a = {12, -65};
    Vec2 b = {13, 2};
    ASSERT_EQ(a.dot(b), 26);
    ASSERT_EQ(b.dot(a), 26);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
