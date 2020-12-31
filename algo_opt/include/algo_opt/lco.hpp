#pragma once

#include <algo_opt/alias.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <exception>
#include <limits>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

namespace algo_opt
{
namespace lco
{
/**
 * @brief data structure representing a linear program
 *
 * @param A MxN constraint matrix for inequality constrain Ax <= b
 * @param b Mx1 constraint matrix for inequality constrain Ax <= b
 * @param c Nx1 objective matrix for min c'x
 *
 * @tparam N int representing problem dimensionality, or number of design
 *         variables
 * @tparam M int representing number of constraints
 */
template <int N, int M>
struct LinearProgram
{
  Matrixd<M, N> A;
  Vectord<M> b;
  Vectord<N> c;
};

// helper functions for implementing the simplex algorithm for solving linear
// constrained optimization problems
namespace simplex
{
/**
 * @brief get x associated with a vertex partition
 *
 * @param[in] B   initial partition of x as array of indices
 * @param[in] LP  Linear Program being solved - must be in equality form
 *
 * @tparam N  int representing problem dimensionality / number of design problem
 * @tparam M  int representing number of constraints - must satisfy M <= N
 *
 * @returns x associated with partition B for LP
 */
template <int N, int M, typename = std::enable_if_t<N >= M>>
Vectord<N> get_vertex(std::array<int, M> B, const LinearProgram<N, M> &LP)
{
  std::ranges::sort(B);

  Matrixd<M> AB = Matrixd<M>::Zero();
  for (size_t i = 0; i < M; ++i)
    AB.col(i) = LP.A.col(B[i]);

  Vectord<M> xB = AB.colPivHouseholderQr().solve(LP.b);

  Vectord<N> x = Vectord<N>::Zero();
  for (size_t i = 0; i < M; ++i)
    x(B[i]) = xB(i);

  return x;
}

/**
 * @brief locate vertex associated with an entering index
 *
 * @param[in] B   initial partition of x as array of indices
 * @param[in] LP  Linear Program being solved - must be in equality form
 * @param[in] q   entering index
 *
 * @tparam N  int representing problem dimensionality / number of design problem
 * @tparam M  int representing number of constraints - must satisfy M <= N
 *
 * @returns tuple of (exiting index, xq)
 */
template <int N, int M, typename = std::enable_if_t<N >= M>>
std::tuple<size_t, double>
edge_transition(std::array<int, M> B, const LinearProgram<N, M> &LP, size_t q)
{
  std::ranges::sort(B);

  auto V = std::vector<int>();
  for (size_t i = 0; i < N; i++)
    if (std::ranges::find(B, i) == B.end())
      V.push_back(i);

  Matrixd<M> AB = Matrixd<M>::Zero();
  for (size_t i = 0; i < M; ++i)
    AB.col(i) = LP.A.col(B[i]);

  Eigen::ColPivHouseholderQR<Matrixd<M>> AB_QR = AB.colPivHouseholderQr();
  Vectord<M> d = AB_QR.solve(LP.A.col(V[q]));
  Vectord<M> xB = AB_QR.solve(LP.b);

  size_t p = 0;
  auto xq = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < M; ++i)
    if (d(i) > 0.0)
    {
      auto v = xB(i) / d(i);
      if (v < xq)
      {
        p = i;
        xq = v;
      }
    }

  return {p, xq};
}

/**
 * @brief attempt to step from one vertex in a LCO problem to a more optimal one
 *
 * @param[in] B   initial partition of x as array of indices
 * @param[in] LP  Linear Program being solved - must be in equality form
 *
 * @tparam N  int representing problem dimensionality / number of design problem
 * @tparam M  int representing number of constraints - must satisfy M <= N
 *
 * @returns next most optimal partition with flag indicating if this partition
 *          matches the input as (B, flag)
 */
template <int N, int M, typename = std::enable_if_t<N >= M>>
std::tuple<std::array<int, M>, bool> step_lp(std::array<int, M> B,
                                             const LinearProgram<N, M> &LP)
{
  std::ranges::sort(B);

  auto V = std::vector<int>();
  for (size_t i = 0; i < N; i++)
    if (std::ranges::find(B, i) == B.end())
      V.push_back(i);

  Matrixd<M> AB = Matrixd<M>::Zero();
  Vectord<M> cB = Vectord<M>::Zero();
  for (size_t i = 0; i < M; ++i)
  {
    AB.col(i) = LP.A.col(B[i]);
    cB(i) = LP.c(B[i]);
  }

  Matrixd<M, N - M> AV = Matrixd<M, N - M>::Zero();
  Vectord<N - M> cV = Vectord<N - M>::Zero();
  for (size_t i = 0; i < N - M; ++i)
  {
    AV.col(i) = LP.A.col(V[i]);
    cV(i) = LP.c(V[i]);
  }

  Vectord<M> lambda = AB.transpose().colPivHouseholderQr().solve(cB);
  Vectord<N - M> muV = cV - AV.transpose() * lambda;

  size_t q = 0;
  size_t p = 0;
  auto xq = std::numeric_limits<double>::infinity();
  auto delta = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < N - M; ++i)
  {
    if (muV(i) < 0.0)
    {
      auto [pi, xi] = edge_transition<N, M>(B, LP, i);
      if (muV(i) * xi < delta)
      {
        q = i;
        p = pi;
        xq = xi;
        delta = muV(i) * xi;
      }
    }
  }

  if (q == 0)
    return {B, true}; // input was the optimal point

  if (xq == std::numeric_limits<double>::infinity())
    throw(std::runtime_error("solution is unbounded"));

  B[p] = V[q];
  return {B, false}; // input was not optimal point, found new vertex
}

/**
 * @brief solve a linear constrained optimization problem using the simplex
 *        method
 *
 * @param[in] B   initial partition of x as array of indices
 * @param[in] LP  Linear Program to solve - must be in equality form
 *
 * @tparam N  int representing problem dimensionality / number of design problem
 * @tparam M  int representing number of constraints - must satisfy M <= N
 *
 * @returns partition which produces x*
 */
template <int N, int M, typename = std::enable_if_t<N >= M>>
std::array<int, M> minimize(std::array<int, M> B, const LinearProgram<N, M> &LP)
{
  auto done = false;

  while (!done)
    std::tie(B, done) = step_lp<N, M>(B, LP);

  return B;
}
} // namespace simplex

/**
 * @brief solve a linear constrained optimization problem using the simplex
 *        method
 *
 * @param[in] LP  Linear Program to solve - must be in equality form
 *
 * @tparam N  int representing problem dimensionality / number of design problem
 * @tparam M  int representing number of constraints - must satisfy M <= N
 *
 * @returns x* (x value which optimizes the linear program)
 */
template <int N, int M, typename = std::enable_if_t<N >= M>>
Vectord<N> minimize_simplex(const LinearProgram<N, M> &LP)
{
  Matrixd<M> Z = Matrixd<M>::Zero();
  for (size_t i = 0; i < M; ++i)
  {
    if (LP.b(i) < 0.0)
      Z(i, i) = -1.0;
    else
      Z(i, i) = 1.0;
  }

  Matrixd<M, N + M> AZ;
  AZ << LP.A, Z;

  Vectord<N + M> cz = Vectord<N + M>::Zero();
  cz.tail(M) = Vectord<M>::Ones();

  auto LP_init = LinearProgram<N + M, M>({AZ, LP.b, cz});
  auto B = std::array<int, M>();
  std::iota(B.begin(), B.end(), N);
  B = simplex::minimize<N + M, M>(B, LP_init);

  // ordinarily, we would want to agument the LP with the z values from the
  // auxillary problem in order to allow partition indices into z values, but it
  // is unclear how to select an additional M indices for the augmented
  // partition when this is the case
  for (const auto &bi : B)
    if (bi >= N)
      throw(std::runtime_error("LP is infeasible"));

  B = simplex::minimize<N, M>(B, LP);
  return simplex::get_vertex<N, M>(B, LP);
}
} // namespace lco
} // namespace algo_opt