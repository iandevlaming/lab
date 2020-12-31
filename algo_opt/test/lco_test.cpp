
#include <algo_opt/alias.hpp>
#include <algo_opt/lco.hpp>

#include <gtest/gtest.h>

#include <algorithm>

namespace lco = ao::lco;
namespace simplex = lco::simplex;

// tests taken from Example 11.7 - Algorithms for Optimization
struct LcoSimplexTest : public ::testing::Test
{
public:
  LcoSimplexTest() : B({2, 3})
  {
    // clang-format off
    ao::Matrixd<M, N> A;
    A << 1,  1, 1, 0,
        -4, 2, 0, 1;
    // clang-format on

    ao::Vectord<M> b;
    b << 9.0, 2.0;

    ao::Vectord<N> c;
    c << 3.0, -1.0, 0.0, 0.0;

    LP = lco::LinearProgram<N, M>({A, b, c});
  }

  static constexpr int N = 4;
  static constexpr int M = 2;

  std::array<int, M> B;

  lco::LinearProgram<N, M> LP;
};

TEST_F(LcoSimplexTest, get_vertex_simple)
{
  ao::Vectord<N> x = simplex::get_vertex<N, M>(B, LP);
  ao::Vectord<N> expected_x;
  expected_x << 0.0, 0.0, LP.b(0), LP.b(1);

  EXPECT_EQ(x, expected_x) << x;
}

TEST_F(LcoSimplexTest, edge_transition_simple)
{
  size_t q = 1;
  auto [p, xq] = simplex::edge_transition<N, M>(B, LP, q);

  size_t expected_p = 1;
  EXPECT_EQ(p, expected_p);

  auto expected_xq = 1.0;
  EXPECT_DOUBLE_EQ(xq, expected_xq);
}

TEST_F(LcoSimplexTest, step_lp_simple)
{
  // first iteration
  auto [B1, f1] = simplex::step_lp<N, M>(B, LP);

  std::ranges::sort(B1);
  auto expected_B1 = std::array<int, 2>({1, 2});

  EXPECT_EQ(B1, expected_B1);
  EXPECT_FALSE(f1);

  // second iteration
  auto [B2, f2] = simplex::step_lp<N, M>(B1, LP);
  EXPECT_EQ(B2, B1);
  EXPECT_TRUE(f2);
}

TEST_F(LcoSimplexTest, minimize_simple)
{
  B = simplex::minimize<N, M>(B, LP);
  std::ranges::sort(B);

  auto expected_B = std::array<int, M>({1, 2});
  EXPECT_EQ(B, expected_B);
}

TEST_F(LcoSimplexTest, minimize_simplex_simple)
{
  ao::Vectord<N> x = minimize_simplex<N, M>(LP);

  ao::Vectord<N> expected_x;
  expected_x << 0.0, 1.0, 8.0, 0.0;
  EXPECT_EQ(x, expected_x) << x;
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}