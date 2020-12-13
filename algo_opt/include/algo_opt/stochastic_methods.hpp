#pragma once

#include <algo_opt/alias.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <numeric>
#include <random>

namespace algo_opt
{
template <int N>
Matrixd<N> rand_lower_basis(double alpha)
{
  // caveat - for the dist to work, alpha needs to be no greater than 4
  if (alpha <= 0.0 || alpha > 4.0)
    throw(std::domain_error("step size must be in the range (0.0, 4.0]"));

  auto delta = static_cast<int>(round(1.0 / sqrt(alpha)));
  Matrixd<N> L = Matrixd<N>::Zero();

  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  auto gen = std::default_random_engine(seed);
  auto sgn_dist = std::bernoulli_distribution(0.5);
  auto disc_dist = std::uniform_int_distribution(std::min(-delta + 1, -0),
                                                 std::max(delta - 1, 0));
  auto get_diagonal = [&]() {
    return static_cast<float>(delta * (sgn_dist(gen) ? 1.0 : -1.0));
  };
  auto get_off_diagonal = [&]() { return static_cast<float>(disc_dist(gen)); };

  for (unsigned int i = 0; i < N; ++i)
  {
    L(i, i) = get_diagonal();
    for (unsigned int j = 0; j < i; ++j)
      L(i, j) = get_off_diagonal();
  }

  return L;
}

template <int N>
std::array<Vectord<N>, N + 1> rand_pos_spanning_set(double alpha)
{
  Matrixd<N> L = rand_lower_basis<N>(alpha);

  auto idx = std::array<unsigned int, N>();
  std::iota(idx.begin(), idx.end(), static_cast<unsigned int>(0));

  auto perm_rows = idx;
  std::random_shuffle(perm_rows.begin(), perm_rows.end());

  auto perm_cols = idx;
  std::random_shuffle(perm_cols.begin(), perm_cols.end());

  Matrixd<N> L_temp = Matrixd<N>::Zero();
  for (const auto i : idx)
    L_temp.row(i) = L.row(perm_rows[i]);

  for (const auto i : idx)
    L.col(i) = L_temp.col(perm_cols[i]);

  Vectord<N> L_extra = Vectord<N>::Zero();
  auto S = std::array<Vectord<N>, N + 1>();
  for (const auto &i : idx)
  {
    L_extra(i) = -1.0 * L.row(i).sum();
    S[i] = L.col(i);
  }
  S.back() = L_extra;

  return S;
}

namespace descent
{
template <typename F, int N, typename = enable_if_multivar_f<F, N>::type>
Vectord<N> mesh_adaptive_direct_search(F f,
                                       const Vectord<N> &x,
                                       double eps = 1e-4,
                                       std::vector<Vectord<N>> *log = nullptr)
{
  auto alpha = 1.0;
  Vectord<N> xi = x;
  auto yi = f(xi);

  while (alpha > eps)
  {
    if (log)
      log->push_back(xi);

    auto improved = false;
    auto S = rand_pos_spanning_set<N>(alpha);
    for (unsigned int i = 0; i < N + 1; ++i)
    {
      // originally, no norm, but I wanna
      Vectord<N> d = S[i].normalized();
      Vectord<N> xn = xi + alpha * d;
      auto yn = f(xn);
      if (improved = (yn < yi))
      {
        xi = xn;
        yi = yn;
        xn = xi + 3 * alpha * d;
        yn = f(xn);
        if (yn < yi)
        {
          xi = xn;
          yi = yn;
        }
        break;
      }
    }
    if (improved)
      alpha = std::min(4.0 * alpha, 1.0);
    else
    {
      alpha /= 4.0;
      if (log)
        log->pop_back();
    }
  }

  if (log)
    log->push_back(xi);

  return xi;
}
} // namespace descent
} // namespace algo_opt