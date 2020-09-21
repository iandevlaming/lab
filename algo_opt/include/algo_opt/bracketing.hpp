#pragma once

#include <algo_opt/alias.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <math.h>
#include <string>
#include <type_traits>
#include <vector>

namespace algo_opt
{
template <typename F, typename = enable_if_unary_f<F>::type>
bracket_t<2> bracket_minimum(F f,
                             double x = 0.0,
                             double s = 0.01,
                             double k = 2.0,
                             std::vector<bracket_t<2>> *log = nullptr)
{
  auto a = Point2d({x, f(x)});
  x = a.x + s;
  auto b = Point2d({x, f(x)});

  if (b.y > a.y)
  {
    std::swap(a, b);
    s *= -1.0;
  }

  auto c = Point2d();
  while (1)
  {
    x = b.x + s;
    c = Point2d({x, f(x)});

    if (log)
      log->push_back({a.x, c.x});

    if (c.y > b.y)
      break;

    a = b;
    b = c;
    s *= k;
  }

  auto bracket = bracket_t<2>({a.x, c.x});
  std::ranges::sort(bracket);
  return bracket;
}

template <typename F, typename = enable_if_unary_f<F>::type>
bracket_t<2>
fibonacci_search(F f, bracket_t<2> bracket, unsigned int n, double eps = 0.01)
{
  eps = abs(eps);
  auto &[a, b] = bracket;

  // from Binet's formula Fn = (phi^n - (1 - phi)^n)/sqrt(5), with phi as the
  // golden ratio
  // ratio of adjacent elements Fi of fibonacci sequence Rn = Fn/Fn-1 = phi * (1
  // - s^(n + 1)) / (1 - s^n)
  static constexpr auto s = (1.0 - sqrt(5.0)) / (1.0 + sqrt(5.0));
  static constexpr auto phi = (1.0 + sqrt(5.0)) / 2.0; // 1 / Rn

  auto rho = 1.0 / (phi * (1.0 - pow(s, n + 1)) / (1.0 - pow(s, n)));
  auto d = rho * b + (1.0 - rho) * a;
  auto yd = f(d);
  auto c = 0.0;
  for (unsigned int i = 1; i < n; ++i)
  {
    if (i == n - 1)
      c = eps * a + (1.0 - eps) * d;
    else
      c = rho * a + (1.0 - rho) * b;
    auto yc = f(c);
    if (yc < yd)
    {
      b = d;
      d = c;
      yd = yc;
    }
    else
    {
      a = b;
      b = c;
    }
    rho = 1.0 / (phi * (1.0 - pow(s, n - i + 1)) / (1.0 - pow(s, n - i)));
  }

  std::ranges::sort(bracket);
  return bracket;
}

template <typename F, typename = enable_if_unary_f<F>::type>
bracket_t<2> golden_section_search(F f, bracket_t<2> bracket, unsigned int n)
{
  auto &[a, b] = bracket;

  static constexpr auto phi = (1.0 + sqrt(5.0)) / 2.0; // 1 / Rn
  auto rho = phi - 1.0;
  auto d = rho * b + (1.0 - rho) * a;
  auto yd = f(d);

  auto c = 0.0;
  for (unsigned int i = 1; i < n; ++i)
  {
    c = rho * a + (1.0 - rho) * b;
    auto yc = f(c);
    if (yc < yd)
    {
      b = d;
      d = c;
      yd = yc;
    }
    else
    {
      a = b;
      b = c;
    }
  }

  std::ranges::sort(bracket);
  return bracket;
}

struct QuadraticFitSearchLog
{
  std::array<double, 3> coefficients;
  bracket_t<3> bracket;
};

template <typename F, typename = enable_if_unary_f<F>::type>
bracket_t<3>
quadratic_fit_search(F f,
                     bracket_t<3> bracket,
                     unsigned int n,
                     std::vector<QuadraticFitSearchLog> *log = nullptr)
{
  std::ranges::sort(bracket);
  assert(std::ranges::adjacent_find(bracket) == bracket.end());
  auto &[a, b, c] = bracket;

  auto ya = f(a);
  auto yb = f(b);
  auto yc = f(c);

  auto fit = [](const auto &x1,
                const auto &x2,
                const auto &x3,
                const auto &y1,
                const auto &y2,
                const auto &y3) {
    // analytically get coefficients (don't have eigen yet)
    auto c2 = (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)) /
              ((x1 - x2) * (x1 - x3) * (x2 - x3));
    auto c1 = (y2 - y1) / (x2 - x1) - c2 * (x1 + x2);
    auto c0 = y1 - c2 * x1 * x1 - c1 * x1;
    return bracket_t<3>({c0, c1, c2});
  };

  for (unsigned int i = 0; i < n; ++i)
  {
    auto coefficients = fit(a, b, c, ya, yb, yc);

    if (log)
      log->push_back(QuadraticFitSearchLog({coefficients, bracket}));

    auto x = -coefficients[1] / (2.0 * coefficients[2]);
    auto yx = f(x);

    if (x > b)
    {
      if (yx > yb)
      {
        c = x;
        yc = yx;
      }
      else
      {
        a = b;
        b = x;
        ya = yb;
        yb = yx;
      }
    }
    else if (x < b)
    {
      if (yx > yb)
      {
        a = x;
        ya = yx;
      }
      else
      {
        c = b;
        b = x;
        yc = yb;
        yb = yx;
      }
    }
    else
      break; // unclear how to continue, so bail
  }

  if (log)
  {
    auto coefficients = fit(a, b, c, ya, yb, yc);
    log->push_back(QuadraticFitSearchLog({coefficients, bracket}));
  }

  return bracket;
}

template <typename F, typename = enable_if_unary_f<F>::type>
bracket_t<2> bisection(F df,
                       bracket_t<2> bracket,
                       double eps,
                       std::vector<bracket_t<2>> *log = nullptr)
{
  eps = abs(eps);
  std::ranges::sort(bracket);
  auto &[a, b] = bracket;
  assert(a != b);

  auto ya = df(a);
  auto yb = df(b);

  if (ya == 0.0)
    b = a;
  if (yb == 0.0)
    a = b;

  while ((b - a) > eps)
  {
    if (log)
      log->push_back({a, b});

    auto x = (a + b) / 2.0;
    auto y = df(x);

    if (y == 0.0) // something something don't use float equality checks
      a = b = x;
    else if ((y * ya) > 0.0) // y and ya have same sign
      a = x;
    else
      b = x;
  }

  if (log)
    log->push_back({a, b});

  return bracket;
}

// obviously, in production-grade code, you'd need additional termination
// guarantees
template <typename F, typename = enable_if_unary_f<F>::type>
double bisection_min(F df, bracket_t<2> bracket, double eps)
{
  eps = abs(eps);
  bracket = bisection(df, bracket, eps);
  return abs(df(bracket[0])) < abs(df(bracket[1])) ? bracket[0] : bracket[1];
}

template <typename F, typename = enable_if_unary_f<F>::type>
double secant_min(F df, bracket_t<2> bracket, double eps)
{
  eps = abs(eps);
  while (abs(bracket[0] - bracket[1]) > eps)
  {
    auto [xn_2, xn_1] = bracket;
    auto xn = xn_1 - df(xn_1) * (xn_1 - xn_2) / (df(xn_1) - df(xn_2));
    bracket = bracket_t<2>({xn_1, xn});
  }

  return abs(df(bracket[0])) < abs(df(bracket[1])) ? bracket[0] : bracket[1];
}

unary_f_t lagrange_poly_fit(const std::vector<Point2d> &points)
{
  assert(points.size() >= 2);
  const auto k = points.size();

  auto x = std::vector<double>(k);
  auto y = std::vector<double>(k);
  for (unsigned int i = 0; i < k; ++i)
  {
    x[i] = points[i].x;
    y[i] = points[i].y;
  }

  auto l = std::vector<unary_f_t>(k);
  for (unsigned int j = 0; j < k; ++j)
  {
    auto lj = std::vector<unary_f_t>();
    const auto &xj = x[j];
    for (unsigned int m = 0; m < k; ++m)
    {
      if (m == j)
        continue;
      const auto &xm = x[m];
      lj.push_back([xj, xm](double xi) { return (xi - xm) / (xj - xm); });
    }
    l[j] = [lj](double xi) {
      auto yi = 1.0;
      std::for_each(lj.begin(), lj.end(), [&](auto &lji) { yi *= lji(xi); });
      return yi;
    };
  }

  auto L = [l, y, k](double xi) {
    auto yi = 0.0;
    for (unsigned int j = 0; j < k; ++j)
      yi += y[j] * l[j](xi);
    return yi;
  };

  return L;
}

unary_f_t inverse_lagrange_poly_fit(const std::vector<Point2d> &points)
{
  auto inv_points = std::vector<Point2d>();
  auto invPoint = [](auto p) { return Point2d({p.y, p.x}); };
  std::transform(
      points.cbegin(), points.cend(), std::back_inserter(inv_points), invPoint);
  return lagrange_poly_fit(inv_points);
}

/**
 * @see https://en.wikipedia.org/wiki/Brent%27s_method
 */
template <typename F, typename = enable_if_unary_f<F>::type>
double brent_min(F df, const bracket_t<2> &bracket, double eps)
{
  eps = abs(eps);
  auto a = Point2d({bracket[0], df(bracket[0])});
  auto b = Point2d({bracket[1], df(bracket[1])});
  assert(a.y * b.y < 0.0); // function must be bracketed

  if (abs(a.y) < abs(b.y))
    std::swap(a, b); // b should be our best guess so far

  auto c = a;
  auto d = c;
  auto s = d;
  auto mflag = true;

  // bisection conditions (god have mercy on me for all this awful code)
  auto c1 = [&]() {
    return s.x >= std::min((3.0 * a.x + b.x) / 4.0, b.x) &&
           s.x <= std::max((3.0 * a.x + b.x) / 4.0, b.x);
  };
  auto c2 = [&]() { return mflag && abs(s.x - b.x) >= abs(b.x - c.x) / 2.0; };
  auto c3 = [&]() { return !mflag && abs(s.x - b.x) >= abs(c.x - c.x) / 2.0; };
  auto c4 = [&]() { return mflag && abs(b.x - c.x) < eps; };
  auto c5 = [&]() { return !mflag && abs(c.x - d.x) < eps; };

  while (abs(b.x - a.x) > eps) // original algo also checks f(b) or f(s) == 0
  {
    if (a.y != c.y && b.y != c.y)
      s.x =
          inverse_lagrange_poly_fit({a, b, c})(0.0); // inverse quadratic interp
    else
      s.x = b.x - b.y * (b.x - a.x) / (b.y - a.y); // secant method

    if (mflag = c1() || c2() || c3() || c4() || c5())
      s.x = (a.x + b.x) / 2.0; // bisection

    s.y = df(s.x);
    d = c;
    c = b;
    if (a.y * s.y < 0.0)
      b = s;
    else
      a = s;

    if (abs(a.y) < abs(b.y))
      std::swap(a, b);
  }

  return abs(b.y) < abs(s.y) ? b.x : s.x;
}

} // namespace algo_opt