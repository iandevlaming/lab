#pragma once

#include <algo_opt/alias.hpp>

#include <Eigen/QR>

#include <algorithm>

namespace algo_opt
{
namespace descent
{
template <typename G,
          typename H,
          int N,
          typename = enable_if_multivar_grad<G, N>::type,
          typename = enable_if_multivar_hess<H, N>::type>
Vectord<N>
newton_step(G g, H h, const Vectord<N> x, double eps, unsigned int k_max)
{
  auto k = static_cast<unsigned int>(1);
  Vectord<N> delta = Vectord<N>::Constant(eps);
  Vectord<N> xi = x;

  while (delta.norm() > eps && k <= k_max)
  {
    Vectord<N> gi = g(x);
    Matrixd<N> hi = h(x);
    delta = hi.colPivHouseholderQr().solve(gi);
    xi -= delta;
    ++k;
  }
  return xi;
}

} // namespace descent
} // namespace algo_opt