#include <algo_dm/alias.hpp>
#include <algo_dm/inference.hpp>

#include <gtest/gtest.h>

#include <vector>

TEST(ProductTest, Example3P1)
{
  const auto x = ad::Variable("x", 2);
  const auto y = ad::Variable("y", 2);
  const auto z = ad::Variable("z", 2);

  const auto f1_table =
      ad::FactorTable(ad::assign({x, y}), std::vector({0.3, 0.4, 0.2, 0.1}));
  const auto f1 = ad::Factor({x, y}, f1_table);

  const auto f2_table =
      ad::FactorTable(ad::assign({y, z}), std::vector({0.2, 0.0, 0.3, 0.5}));
  const auto f2 = ad::Factor({y, z}, f2_table);

  auto expected_table = ad::FactorTable(
      ad::assign({x, y, z}),
      std::vector({0.06, 0.00, 0.12, 0.20, 0.04, 0.00, 0.03, 0.05}));
  auto expected_product = ad::Factor({x, y, z}, expected_table);

  auto product = f1 * f2;

  EXPECT_EQ(product, expected_product);
  EXPECT_EQ(product, f2 * f1);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}