#pragma once

#include <algo_opt/alias.hpp>
#include <algo_opt/bracketing.hpp>

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <optional>
#include <type_traits>

namespace algo_opt
{
binary_f_t makeRosenbrock(double a = 1, double b = 100)
{
  // f = (a - x)^2 + b * (y - x^2)^2
  return [a, b](const Eigen::Vector2d &pt) {
    const auto &x = pt(0);
    const auto &y = pt(1);
    return std::pow((a - x), 2) + b * std::pow((y - x * x), 2);
  };
}

// auto diff (or even getting funclib working) would be a better move _but_ here
// we are
binary_grad_t makeGradRosenbrock(double a = 1, double b = 100)
{
  // df/dx = -2 * (a - x) - 4 * x * b * (y - x^2)
  // df/dy = 2 * b * (y - x^2)
  return [a, b](const Eigen::Vector2d &pt) -> Eigen::Vector2d {
    const auto &x = pt(0);
    const auto &y = pt(1);
    const auto dfdx = -2.0 * (a - x) - 4.0 * x * b * (y - x * x);
    const auto dfdy = 2.0 * b * (y - x * x);
    Eigen::Vector2d grad(dfdx, dfdy);
    return grad;
  };
}

binary_level_t makeLevelRosenbrock(double a = 1, double b = 100)
{
  return [a, b](double x, double f) -> std::optional<std::array<double, 2>> {
    // y = sqrt((f - (a - x)^2) / b) + x^2
    auto g = (f - pow((a - x), 2)) / b;
    if (g < 0.0)
      return {};
    auto h = std::sqrt(g);
    auto x2 = x * x;
    return std::array<double, 2>({-h + x2, h + x2});
  };
}

template <typename F,
          typename G,
          int N,
          typename = enable_if_multivar_f<F, N>::type,
          typename = enable_if_multivar_grad<G, N>::type>
Vectord<N> line_search(F f,
                       G df,
                       const Vectord<N> &x,
                       const Vectord<N> &d,
                       double eps = 1e-4)
{
  auto f_alpha = [&](double alpha) { return f(x + alpha * d); };
  auto df_alpha = [&](double alpha) { return df(x + alpha * d).dot(d); };
  auto bracket = bracket_minimum(f_alpha);
  auto alpha = brent_min(df_alpha, bracket, eps);
  return x + alpha * d;
}

template <typename F,
          typename G,
          int N,
          typename = enable_if_multivar_f<F, N>::type,
          typename = enable_if_multivar_grad<G, N>::type>
Eigen::Matrix<double, N, 1>
backtracking_line_search(F f,
                         G df,
                         const Eigen::Matrix<double, N, 1> &x,
                         const Eigen::Matrix<double, N, 1> &d,
                         double alpha,
                         double p = 0.5,
                         double beta = 1e-4)
{
  const auto y = f(x);
  const Eigen::Matrix<double, N, 1> dydx = df(x);
  while (f(x + alpha * d) > y + beta * alpha * dydx.dot(d))
    alpha *= p;
  return x + alpha * d;
}

template <typename F,
          typename G,
          int N,
          typename = std::enable_if_t<is_multivar_f<F, N>::value &&
                                      is_multivar_grad<G, N>::value>>
Eigen::Matrix<double, N, 1>
strong_backtracking_line_search(F f,
                                G df,
                                const Eigen::Matrix<double, N, 1> &x,
                                const Eigen::Matrix<double, N, 1> &d,
                                double alpha = 1.0,
                                double beta = 1e-4,
                                double sigma = 0.1)
{
  const auto y0 = f(x);
  auto y_prev = f(x + alpha * d) - 1.0; // set lower than initial y
  const auto g0 = df(x).dot(d);
  auto alpha_prev = 0.0;
  auto alpha_lo = 0.0;
  auto alpha_hi = 0.0;

  // bracketing phase - get a region which is guaranteed to contain the strong
  // wolfe condition
  while (1)
  {
    const auto y = f(x + alpha * d);
    // leave bracketing if violated first wolfe condition or overshot local min
    if ((y > y0 + beta * alpha * g0) || y >= y_prev)
    {
      alpha_lo = alpha_prev;
      alpha_hi = alpha;
      break;
    }

    const auto g = df(x + alpha * d).dot(d);
    // done if satisfy strong wolfe condition
    if (abs(g) <= -sigma * g0)
      return x + alpha * d;

    // leave bracketing if overshot local min
    if (g >= 0.0)
    {
      alpha_lo = alpha;
      alpha_hi = alpha_prev;
      break;
    }

    y_prev = y;
    alpha_prev = alpha;
    alpha *= 2.0; // magic number, alpha increase factor
  }

  // zoom phase - bisect alpha region until strong wolfe condition is met
  const auto y_lo = f(x + alpha_lo * d);
  while (1)
  {
    alpha = (alpha_lo + alpha_hi) / 2.0;
    const auto y = f(x + alpha * d);
    if ((y > y0 + beta * alpha * g0) || y >= y_lo)
      alpha_hi = alpha;
    else
    {
      const auto g = df(x + alpha * d).dot(d);
      if (abs(g) <= -sigma * g0 || (alpha_hi - alpha_lo < 1e-4))
        return x + alpha * d;
      if (g * (alpha_hi - alpha_lo) >= 0.0)
        alpha_hi = alpha_lo;
      alpha_lo = alpha;
    }
  }
}

// the trust region version requires a minimization solver. we could use
// quadratic fit solver, but would prefer another method - but I need to provide
// a derivative function approximation. Might want to implementa a taylor
// expansion algo
} // namespace algo_opt