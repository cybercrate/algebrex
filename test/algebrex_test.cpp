#include <gtest/gtest.h>

#include "algebrex.h"

TEST(algebrex_test, simple_expr_1) {
    const auto result = expr_solver{}.solve("2 + 2 * 2");
    const auto expected = "6";

    EXPECT_EQ(result, expected);
}

TEST(algebrex_test, simple_expr_2) {
    const auto result = expr_solver{}.solve("(2 + 2) * 2");
    const auto expected = "8";

    EXPECT_EQ(result, expected);
}

TEST(algebrex_test, complex_expr) {
    const auto result = expr_solver{}.solve("1 + ((2 - 3 * 4) / pi) ^ 6");
    const auto expected = "1041.161744";

    EXPECT_EQ(result, expected);
}

TEST(algebrex_test, error) {
    EXPECT_THROW(expr_solver{}.solve("2 +"), std::exception);
    EXPECT_THROW(expr_solver{}.solve("(2 +)"), std::exception);
    EXPECT_THROW(expr_solver{}.solve("2..2"), std::exception);
}