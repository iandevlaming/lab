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
template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N> central_diff(F f, const Vectord<N> &x, double h = 1e-4)
{
  Vectord<N> df = Vectord<N>::Zero();
  for (unsigned int i = 0; i < N; ++i)
  {
    Vectord<N> xa = x;
    xa(i) += h / 2.0;
    Vectord<N> xb = x;
    xb(i) -= h / 2.0;
    df(i) = (f(xa) - f(xb)) / h;
  }
  return df;
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

template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N>
line_search(F f, const Vectord<N> &x, const Vectord<N> &d, double eps = 1e-4)
{
  auto f_alpha = [&](double alpha) { return f(x + alpha * d); };
  auto df_alpha = [&](double alpha) {
    Vectord<N> xn = x + alpha * d;
    Vectord<N> df = central_diff(f, xn, eps);
    return df.dot(d);
  };
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