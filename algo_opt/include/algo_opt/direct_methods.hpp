#pragma once

#include <algo_opt/alias.hpp>
#include <algo_opt/local_descent.hpp>

#include <Eigen/QR>

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <vector>

namespace algo_opt
{
namespace descent
{
template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N> cyclic_coordinate_descent(F f,
                                     const Vectord<N> x,
                                     double eps = 1e-4,
                                     std::vector<Vectord<N>> *log = nullptr)
{
  if (log)
    log->push_back(x);

  auto delta = eps + 1.0;
  Vectord<N> xn = x;

  while (delta > eps)
  {
    Vectord<N> xi = xn;
    for (unsigned int i = 0; i < N; ++i)
    {
      Vectord<N> d = Vectord<N>::Zero();
      d(i) = 1.0;
      xn = line_search(f, xn, d);

      if (log)
        log->push_back(xn);
    }
    Vectord<N> x_diff = xn - xi;
    delta = x_diff.norm();
  }
  return xn;
}

template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N>
accelerated_cyclic_coordinate_descent(F f,
                                      const Vectord<N> x,
                                      double eps = 1e-4,
                                      std::vector<Vectord<N>> *log = nullptr)
{
  if (log)
    log->push_back(x);

  auto delta = eps + 1.0;
  Vectord<N> xn = x;

  while (delta > eps)
  {
    Vectord<N> xi = xn;
    for (unsigned int i = 0; i < N; ++i)
    {
      Vectord<N> d = Vectord<N>::Zero();
      d(i) = 1.0;
      xn = line_search(f, xn, d);

      if (log)
        log->push_back(xn);
    }

    Vectord<N> x_diff = xn - xi;
    xn = line_search(f, xn, x_diff);

    if (log)
      log->push_back(xn);

    x_diff = xn - xi;
    delta = x_diff.norm();
  }
  return xn;
}

template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N> powell(F f,
                  const Vectord<N> &x,
                  double eps = 1e-4,
                  std::vector<Vectord<N>> *log = nullptr)
{
  auto delta = eps + 1.0;
  Vectord<N> xn = x;

  if (log)
    log->push_back(xn);

  std::vector<Vectord<N>> u(N, Vectord<N>::Zero());
  for (unsigned int i = 0; i < N; ++i)
    u[i](i) = 1.0;

  while (delta > eps)
  {
    Vectord<N> xi = xn;
    for (unsigned int i = 0; i < N; ++i)
    {
      xi = line_search(f, xi, u[i]);

      if (log)
        log->push_back(xi);
    }

    u.erase(u.cbegin());
    Vectord<N> di = xi - xn;
    di.normalize();
    u.push_back(di);

    xi = line_search(f, xi, di);
    Vectord<N> x_diff = xi - xn;
    delta = x_diff.norm();
    xn = xi;

    if (log)
      log->push_back(xi);
  }

  return xn;
}

template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N> hooke_jeeves(F f,
                        const Vectord<N> &x,
                        double alpha,
                        double eps = 1e-4,
                        double gamma = 0.5,
                        std::vector<std::vector<Vectord<N>>> *log = nullptr)
{
  Vectord<N> xn = x;
  auto yn = f(xn);
  while (alpha > eps)
  {
    auto improved = false;
    Vectord<N> x_opt = xn;
    auto y_opt = yn;

    if (log)
      log->push_back(std::vector<Vectord<N>>({x_opt}));

    for (unsigned int i = 0; i < N; ++i)
    {
      for (auto sgn : {-1.0, 1.0})
      {
        Vectord<N> d = Vectord<N>::Zero();
        d(i) = sgn;
        Vectord<N> xi = xn + alpha * d;
        auto yi = f(xi);
        if (yi < y_opt)
        {
          x_opt = xi;
          y_opt = yi;
          improved = true;
        }

        if (log)
          log->back().push_back(xi);
      }
    }
    xn = x_opt;
    yn = y_opt;

    if (!improved)
    {
      alpha *= gamma;

      if (log)
        log->pop_back();
    }
  }

  return xn;
}

template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N>
generalized_pattern_search(F f,
                           const Vectord<N> &x,
                           double alpha,
                           const std::vector<Vectord<N>> &D,
                           double eps = 1e-4,
                           double gamma = 0.5,
                           std::vector<std::vector<Vectord<N>>> *log = nullptr)
{
  Vectord<N> xn = x;
  auto yn = f(xn);
  auto d = D;
  while (alpha > eps)
  {
    auto improved = false;

    if (log)
      log->push_back(std::vector<Vectord<N>>({xn}));

    for (unsigned int i = 0; i < d.size(); ++i)
    {
      Vectord<N> di = d[i];
      Vectord<N> xi = xn + alpha * di;
      auto yi = f(xi);

      if (log)
        log->back().push_back(xi);

      if (yi < yn)
      {
        xn = xi;
        yn = yi;
        improved = true;

        d.erase(d.cbegin() + i);
        d.insert(d.cbegin(), di);
        break;
      }
    }
    if (!improved)
    {
      alpha *= gamma;

      if (log)
        log->pop_back();
    }
  }
  return xn;
}

template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N>
nelder_mead(F f,
            const std::array<Vectord<N>, N + 1> &S,
            double epsilon,
            double alpha = 1.0,
            double beta = 2.0,
            double gamma = 0.5,
            std::vector<std::array<Vectord<N>, N + 1>> *log = nullptr)
{
  auto s = S;
  auto delta = std::numeric_limits<double>::infinity();
  auto y_arr = std::array<double, N + 1>();
  std::ranges::transform(s, y_arr.begin(), f);

  while (delta > epsilon)
  {
    std::ranges::sort(y_arr);
    std::ranges::sort(s, {}, f);

    if (log)
      log->push_back(s);

    // could use const refs here effectively, but don't want to get confused
    // get lowest
    Vectord<N> xl = s.front();
    auto yl = y_arr.front();
    // get highest
    Vectord<N> xh = s.back();
    auto yh = y_arr.back();
    // get second highest
    auto ys = y_arr[N - 2];

    // I should really just be taking advantage of Eigen here but oh well
    // centroid excluding highest point
    Vectord<N> init = Vectord<N>::Zero();
    Vectord<N> xbar = std::accumulate(s.cbegin(), s.cend() - 1, init) / N;
    // reflection point about centroid
    Vectord<N> xr = xbar + alpha * (xbar - xh);
    auto yr = f(xr);

    // expand if reflection descends
    if (yr < yl)
    {
      Vectord<N> xe = xbar + beta * (xr - xbar);
      auto ye = f(xe);

      s.back() = ye < yr ? xe : xr;
      y_arr.back() = ye < yr ? ye : yr;
    }
    else if (yr > ys)
    {
      // shrink highest point if possible
      if (yr <= yh)
      {
        s.back() = xh = xr;
        y_arr.back() = yh = yr;
      }
      // contract
      Vectord<N> xc = xbar + gamma * (xh - xbar);
      auto yc = f(xc);
      if (yc > yh)
      {
        std::for_each(s.begin() + 1, s.end(), [&](Vectord<N> &si) {
          si = (si + xl) / 2.0;
        });
        std::transform(s.cbegin() + 1, s.cend(), y_arr.begin() + 1, f);
      }
      else
      {
        s.back() = xc;
        y_arr.back() = yc;
      }
    }
    else
    {
      s.back() = xr;
      y_arr.back() = yr;
    }

    // delta is the standard deviation of y
    // this is going to suck to compute here
    auto y_mean =
        std::accumulate(y_arr.cbegin(), y_arr.cend(), 0.0) / y_arr.size();
    auto y_dev = y_arr;
    std::ranges::for_each(
        y_dev, [&](auto &yi) { yi = (yi - y_mean) * (yi - y_mean); });
    auto y_std = sqrt(std::accumulate(y_dev.cbegin(), y_dev.cend(), 0.0) / N);

    delta = y_std;
  }

  if (log)
    log->push_back(s);

  return *std::ranges::min_element(s, {}, f);
}
} // namespace descent
} // namespace algo_opt