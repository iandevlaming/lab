
#include <algo_opt/bracketing.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>

namespace ao = algo_opt;

class BracketingTest : public ::testing::Test
{
public:
  BracketingTest()
      : f_([](double x) { return -cos(M_PI * x / 10.0); }),
        df_([](double x) { return M_PI / 10.0 * sin(M_PI * x / 10.0); }),
        x_min_(0.0)
  {
  }

protected:
  std::function<double(double)> f_;
  std::function<double(double)> df_;
  double x_min_;
};

TEST_F(BracketingTest, bracketing_minimum_simple)
{
  const auto x0 = -1.0;
  const auto s = 0.01;
  const auto k = 2.0;
  auto [a, b] = ao::bracket_minimum(f_, x0, s, k);

  EXPECT_LT(a, x_min_);
  EXPECT_GT(b, x_min_);
}

TEST_F(BracketingTest, fibonacci_search_simple)
{
  const auto x0 = -1.0;
  const auto s = 0.01;
  const auto k = 2.0;
  auto bracket = ao::bracket_minimum(f_, x0, s, k);

  auto [c1, d1] = ao::fibonacci_search(f_, bracket, 2);
  auto [c2, d2] = ao::fibonacci_search(f_, bracket, 3);

  EXPECT_LT(c1, x_min_);
  EXPECT_LT(c2, x_min_);
  EXPECT_GT(d1, x_min_);
  EXPECT_GT(d2, x_min_);

  EXPECT_LT(d2 - c2, d1 - c1);
}

TEST_F(BracketingTest, quadratic_fit_search_sorted_output)
{
  auto bracket = std::array<double, 3>({1.0, 0.3, -1.0});
  auto n = 3;
  auto fit_bracket = ao::quadratic_fit_search(f_, bracket, n);

  EXPECT_TRUE(std::ranges::is_sorted(fit_bracket))
      << fit_bracket[0] << " " << fit_bracket[1] << " " << fit_bracket[2];
}

TEST_F(BracketingTest, quadratic_fit_search_simple)
{
  auto bracket = std::array<double, 3>({1.0, 0.3, -1.0});
  auto n = 3;

  auto fit_bracket = ao::quadratic_fit_search(f_, bracket, n);
  EXPECT_LT(fabs(fit_bracket[1] - x_min_), fabs(bracket[1] - x_min_));
}

TEST_F(BracketingTest, quadratic_fit_search_tighten)
{
  auto bracket = std::array<double, 3>({1.0, 0.3, -1.0});
  auto small_n = 3;
  auto big_n = small_n * 2;

  auto less_fit_bracket = ao::quadratic_fit_search(f_, bracket, small_n);
  auto more_fit_bracket = ao::quadratic_fit_search(f_, bracket, big_n);

  EXPECT_LT(more_fit_bracket[2] - more_fit_bracket[0],
            less_fit_bracket[2] - less_fit_bracket[0]);
}

TEST_F(BracketingTest, bisection_simple)
{
  auto bracket = std::array<double, 2>({-2.0, 1.0});
  auto eps = 0.001;

  auto [a, b] = ao::bisection(df_, bracket, eps);

  EXPECT_LE(a, x_min_);
  EXPECT_GE(b, x_min_);
  EXPECT_LT(b - a, eps);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}