#pragma once

#include <algo_opt/bracketing.hpp>
#include <algo_opt/local_descent.hpp>

#include <algorithm>
#include <vector>

namespace algo_opt
{
namespace descent
{
struct GradientDescentParams
{
  double alpha{1.0};
};

struct GradientDescentOptParams
{
  double eps{1e-4};
};

template <int N>
struct ConjugateGradientDescentParams
{
  ConjugateGradientDescentParams() = default;
  template <typename G, typename = enable_if_multivar_grad<G, N>::type>
  ConjugateGradientDescentParams(G grad,
                                 const Vectord<N> &x,
                                 double epsilon = 1e-4)
      : eps(epsilon)
  {
    g = grad(x);
    d = -g;
  }

  Vectord<N> g;
  Vectord<N> d;
  double eps{1e-4};
};

template <int N>
struct MomentumParams
{
  MomentumParams() { v = Vectord<N>::Zero(); }
  MomentumParams(double alpha_in, double beta_in, const Vectord<N> &v_in)
      : alpha(alpha_in), beta(beta_in), v(v_in)
  {
  }
  double alpha{0.0};
  double beta{0.0};
  Vectord<N> v;
};

template <int N>
struct NesterovMomentumParams
{
  NesterovMomentumParams() { v = Vectord<N>::Zero(); }
  NesterovMomentumParams(double alpha_in,
                         double beta_in,
                         const Vectord<N> &v_in)
      : alpha(alpha_in), beta(beta_in), v(v_in)
  {
  }
  double alpha{0.0};
  double beta{0.0};
  Vectord<N> v;
};

template <int N>
struct AdagradParams
{
  AdagradParams() { s = Vectord<N>::Zero(); }
  AdagradParams(double alpha_in, double eps_in, const Vectord<N> &s_in)
      : alpha(alpha_in), eps(eps_in), s(s_in)
  {
  }
  double alpha{0.01};
  double eps{1e-4};
  Vectord<N> s;
};

template <int N>
struct RMSPropParams
{
  RMSPropParams() { s = Vectord<N>::Zero(); }
  RMSPropParams(double alpha_in,
                double gamma_in,
                double eps_in,
                const Vectord<N> &s_in)
      : alpha(alpha_in), gamma(gamma_in), eps(eps_in), s(s_in)
  {
  }
  double alpha{0.01};
  double gamma{0.9};
  double eps{1e-4};
  Vectord<N> s;
};

template <int N>
struct AdadeltaParams
{
  AdadeltaParams()
  {
    s = Vectord<N>::Zero();
    u = Vectord<N>::Zero();
  }
  AdadeltaParams(double gamma_s_in,
                 double gamma_x_in,
                 double eps_in,
                 const Vectord<N> &s_in,
                 const Vectord<N> &u_in)
      : gamma_s(gamma_s_in), gamma_x(gamma_x_in), eps(eps_in), s(s_in), u(u_in)
  {
  }
  double gamma_s{0.9};
  double gamma_x{0.9};
  double eps{1e-4};
  Vectord<N> s;
  Vectord<N> u;
};

template <int N>
struct AdamParams
{
  AdamParams()
  {
    s = Vectord<N>::Zero();
    v = Vectord<N>::Zero();
  }
  AdamParams(double alpha_in,
             double gamma_s_in,
             double gamma_v_in,
             double eps_in,
             unsigned int k_in,
             const Vectord<N> &s_in,
             const Vectord<N> &v_in)
      : alpha(alpha_in),
        gamma_s(gamma_s_in),
        gamma_v(gamma_v_in),
        eps(eps_in),
        k(k_in),
        s(s_in),
        v(v_in)
  {
  }
  double alpha{0.001};
  double gamma_s{0.999};
  double gamma_v{0.9};
  double eps{1e-8};
  unsigned int k{0};
  Vectord<N> s;
  Vectord<N> v;
};

template <int N>
struct HyperGradientDescentParams
{
  HyperGradientDescentParams() { g = Vectord<N>::Zero(); }
  HyperGradientDescentParams(double alpha_in,
                             double mu_in,
                             const Vectord<N> &g_in)
      : alpha(alpha_in), mu(mu_in), g(g_in)
  {
  }
  double alpha{0.01};
  double mu{1e-6};
  Vectord<N> g;
};

template <int N>
struct DFPParams
{
  DFPParams() { q = Matrixd<N>::Identity(); }
  DFPParams(const Matrixd<N> &q_in) : q(q_in) {}
  Matrixd<N> q;
};

template <int N>
struct BFGSParams
{
  BFGSParams() { q = Matrixd<N>::Identity(); }
  BFGSParams(const Matrixd<N> &q_in) : q(q_in) {}
  Matrixd<N> q;
};

template <int N>
struct LimitedMemoryBFGSParams
{
  unsigned int m_max{0};
  std::vector<Vectord<N>> delta_s;
  std::vector<Vectord<N>> gamma_s;
  std::vector<Vectord<N>> q_s;
};

template <typename S, typename P, typename = enable_if_unary_f<S>::type>
struct NoisyDescentParams
{
  NoisyDescentParams()
  {
    sigma = [](const auto &ki) { return 1.0 / ki; };
  }
  NoisyDescentParams(const P &params_in, const S &sigma_in, unsigned int k_in)
      : descent_params(params_in), sigma(sigma_in), k(k_in)
  {
  }
  P descent_params;
  S sigma;
  unsigned int k = 1;
};

// we can use dispatch if we match function interfaces, but the interface for
// each function is not the same... the goto solution here, since each step
// requires some internal state, is probably an abstract class and a factory. I
// like the functional approach of being able to inspect input and predict
// output but that means we have to track a state var anyway... some of these
// need to modify the param struct... well, let's try using it and see how
// frustrating it is
template <typename F,
          typename G,
          int N,
          typename = enable_if_multivar_f<F, N>::type,
          typename = enable_if_multivar_grad<G, N>::type>
struct Descent
{
  static std::tuple<Vectord<N>, GradientDescentParams>
  step(const GradientDescentParams &params, F, G g, const Vectord<N> &x)
  {
    return {x - (params.alpha * g(x)), params};
  }

  static std::tuple<Vectord<N>, GradientDescentOptParams>
  step(const GradientDescentOptParams &params, F f, G g, const Vectord<N> &x)
  {
    Vectord<N> d = -1.0 * g(x).normalized();
    Vectord<N> xi = line_search(f, g, x, d, params.eps);
    return {xi, params};
  }

  static std::tuple<Vectord<N>, ConjugateGradientDescentParams<N>>
  step(const ConjugateGradientDescentParams<N> &params,
       F f,
       G g,
       const Vectord<N> &x)
  {
    const Vectord<N> &g_prev = params.g;
    const Vectord<N> &d_prev = params.d;

    Vectord<N> g_next = g(x);
    auto beta = (g_next.dot(g_next - g_prev)) /
                (g_prev.dot(g_prev)); // polak-ribiere update
    beta = std::max(0.0, beta);

    Vectord<N> d_next = -g_next + beta * d_prev;
    Vectord<N> x_next = line_search(f, g, x, d_next, params.eps);

    auto new_params = params;
    new_params.g = g_next;
    new_params.d = d_next;

    return {x_next, new_params};
  }

  static std::tuple<Vectord<N>, MomentumParams<N>>
  step(const MomentumParams<N> &params, F, G g, const Vectord<N> &x)
  {
    const auto &alpha = params.alpha;
    const auto &beta = params.beta;
    const Vectord<N> &v_prev = params.v;

    // could do better I think with optimal alpha, but that's not the algo
    Vectord<N> gi = g(x);
    Vectord<N> vi = beta * v_prev - alpha * gi;

    auto new_params = MomentumParams<N>(alpha, beta, vi);
    return {x + vi, new_params};
  }

  static std::tuple<Vectord<N>, NesterovMomentumParams<N>>
  step(const NesterovMomentumParams<N> &params, F, G g, const Vectord<N> &x)
  {
    const auto &alpha = params.alpha;
    const auto &beta = params.beta;
    const Vectord<N> &v_prev = params.v;

    // could do better I think with optimal alpha, but that's not the algo
    Vectord<N> vi = beta * v_prev - alpha * g(x + beta * v_prev);

    auto new_params = NesterovMomentumParams<N>(alpha, beta, vi);
    return {x + vi, new_params};
  }

  static std::tuple<Vectord<N>, AdagradParams<N>>
  step(const AdagradParams<N> &params, F, G g, const Vectord<N> &x)
  {
    const auto &alpha = params.alpha;
    const auto &eps = params.eps;
    auto s = params.s;

    Vectord<N> gi = g(x);
    s += gi.cwiseProduct(gi);
    Vectord<N> xi = x;
    for (unsigned int i = 0; i < N; ++i)
      xi(i) -= alpha * gi(i) / (eps + sqrt(s(i)));

    auto new_params = AdagradParams<N>(alpha, eps, s);
    return {xi, new_params};
  }

  static std::tuple<Vectord<N>, RMSPropParams<N>>
  step(const RMSPropParams<N> &params, F, G g, const Vectord<N> &x)
  {
    const auto &alpha = params.alpha;
    const auto &gamma = params.gamma;
    const auto &eps = params.eps;
    Vectord<N> s = params.s;

    Vectord<N> gi = g(x);
    s = gamma * s + (1.0 - gamma) * gi.cwiseProduct(gi);
    Vectord<N> xi = x;
    for (unsigned int i = 0; i < N; ++i)
      xi(i) -= alpha * gi(i) / (eps + sqrt(s(i)));

    auto new_params = RMSPropParams<N>(alpha, gamma, eps, s);
    return {xi, new_params};
  }

  static std::tuple<Vectord<N>, AdadeltaParams<N>>
  step(const AdadeltaParams<N> &params, F, G g, const Vectord<N> &x)
  {
    const auto &gamma_s = params.gamma_s;
    const auto &gamma_x = params.gamma_x;
    const auto &eps = params.eps;
    Vectord<N> s = params.s;
    Vectord<N> u = params.u;

    Vectord<N> gi = g(x);
    s = gamma_s * s + (1.0 - gamma_s) * gi.cwiseProduct(gi);
    Vectord<N> delta_x = Vectord<N>::Zero();
    for (unsigned int i = 0; i < N; ++i)
      delta_x(i) -= (sqrt(u(i)) + eps) * gi(i) / (sqrt(s(i)) + eps);
    u = gamma_x * u + (1.0 - gamma_x) * delta_x.cwiseProduct(delta_x);

    auto new_params = AdadeltaParams<N>(gamma_s, gamma_x, eps, s, u);
    return {x + delta_x, new_params};
  }

  static std::tuple<Vectord<N>, AdamParams<N>>
  step(const AdamParams<N> &params, F, G g, const Vectord<N> &x)
  {
    const auto &alpha = params.alpha;
    const auto &gamma_s = params.gamma_s;
    const auto &gamma_v = params.gamma_v;
    const auto &eps = params.eps;
    auto k = params.k + 1;
    Vectord<N> s = params.s;
    Vectord<N> v = params.v;

    Vectord<N> gi = g(x);
    s = gamma_s * s + (1.0 - gamma_s) * gi.cwiseProduct(gi);
    v = gamma_v * v + (1.0 - gamma_v) * gi;
    Vectord<N> s_hat = s * (1.0 / (1.0 - pow(gamma_s, k)));
    Vectord<N> v_hat = v * (1.0 / (1.0 - pow(gamma_v, k)));
    Vectord<N> xi = x;
    for (unsigned int i = 0; i < N; ++i)
      xi(i) -= alpha * v_hat(i) / (sqrt(s_hat(i)) + eps);

    auto new_params = AdamParams<N>(alpha, gamma_s, gamma_v, eps, k, s, v);
    return {xi, new_params};
  }

  static std::tuple<Vectord<N>, HyperGradientDescentParams<N>>
  step(const HyperGradientDescentParams<N> &params, F, G g, const Vectord<N> &x)
  {
    auto alpha = params.alpha;
    const auto &mu = params.mu;
    const Vectord<N> &g_prev = params.g;

    Vectord<N> gi = g(x);
    alpha += mu * gi.dot(g_prev);

    auto new_params = HyperGradientDescentParams<N>(alpha, mu, gi);
    return {x - alpha * gi, new_params};
  }

  static std::tuple<Vectord<N>, DFPParams<N>>
  step(const DFPParams<N> &params, F f, G g, const Vectord<N> &x)
  {
    Matrixd<N> q = params.q;

    Vectord<N> gi = g(x);
    Vectord<N> d = -q * gi;
    Vectord<N> xn = line_search(f, g, x, d);
    Vectord<N> gn = g(xn);
    Vectord<N> delta = xn - x;
    Vectord<N> gamma = gn - gi;
    q -= (q * gamma * gamma.transpose() * q) / (gamma.transpose() * q * gamma);
    q += (delta * delta.transpose()) / (delta.transpose() * gamma);

    auto new_params = DFPParams(q);
    return {xn, new_params};
  }

  static std::tuple<Vectord<N>, BFGSParams<N>>
  step(const BFGSParams<N> &params, F f, G g, const Vectord<N> &x)
  {
    Matrixd<N> q = params.q;

    Vectord<N> gi = g(x);
    Vectord<N> d = -q * gi;
    Vectord<N> xn = line_search(f, g, x, d);
    Vectord<N> gn = g(xn);
    Vectord<N> delta = xn - x;
    Vectord<N> gamma = gn - gi;
    q -= (delta * gamma.transpose() * q + q * gamma * delta.transpose()) /
         (delta.transpose() * gamma);
    q += (1.0 +
          (gamma.transpose() * q * gamma)(0) / (delta.transpose() * gamma)) *
         (delta * delta.transpose()) / (delta.transpose() * gamma);

    auto new_params = BFGSParams(q);
    return {xn, new_params};
  }

  static std::tuple<Vectord<N>, LimitedMemoryBFGSParams<N>>
  step(const LimitedMemoryBFGSParams<N> &params, F f, G g, const Vectord<N> &x)
  {
    std::vector<Vectord<N>> delta_s = params.delta_s;
    std::vector<Vectord<N>> gamma_s = params.gamma_s;
    std::vector<Vectord<N>> q_s = params.q_s;
    auto m_max = params.m_max;

    Vectord<N> d = Vectord<N>::Zero();
    Vectord<N> xn = Vectord<N>::Zero();
    Vectord<N> gi = g(x);
    auto m = static_cast<int>(delta_s.size());

    if (m > 0)
    {
      Vectord<N> q = gi;
      for (auto i = m - 1; i >= 0; --i)
      {
        q_s[i] = q;
        q -= (delta_s[i].dot(q)) / (gamma_s[i].dot(delta_s[i])) * gamma_s[i];
      }

      Vectord<N> z =
          (gamma_s.back().cwiseProduct(delta_s.back()).cwiseProduct(q)) /
          (gamma_s.back().dot(gamma_s.back()));
      for (auto i = 0; i < m; ++i)
        z += delta_s[i] * (delta_s[i].dot(q_s[i]) - gamma_s[i].dot(z)) /
             (gamma_s[i].dot(delta_s[i]));
      d = -1.0 * z;
    }
    else
    {
      d = -1.0 * gi;
    }

    xn = line_search(f, g, x, d);
    Vectord<N> gn = g(xn);

    delta_s.push_back(xn - x);
    gamma_s.push_back(gn - gi);
    q_s.push_back(Vectord<N>::Zero());

    auto n_remove =
        std::max(static_cast<int>(delta_s.size()) - static_cast<int>(m_max), 0);
    delta_s.erase(delta_s.begin(), delta_s.begin() + n_remove);
    gamma_s.erase(gamma_s.begin(), gamma_s.begin() + n_remove);
    q_s.erase(q_s.begin(), q_s.begin() + n_remove);

    auto new_params = LimitedMemoryBFGSParams<N>(m_max, delta_s, gamma_s, q_s);
    return {xn, new_params};
  }

  template <typename S, typename P>
  static std::tuple<Vectord<N>, NoisyDescentParams<S, P>>
  step(const NoisyDescentParams<S, P> &params, F f, G g, const Vectord<N> &x)
  {
    const auto &descent_params = params.descent_params;
    const auto &sigma = params.sigma;
    const auto &k = params.k;

    Vectord<N> xn;
    P new_descent_params;
    std::tie(xn, new_descent_params) =
        Descent<F, G, N>::step(descent_params, f, g, x);
    xn += sigma(k) * Vectord<N>::Random();

    auto new_params = NoisyDescentParams(new_descent_params, sigma, k + 1);
    return {xn, new_params};
  }
};

template <typename P, typename F, typename G, int N>
std::tuple<Vectord<N>, P> step(const P &params, F f, G g, const Vectord<N> &x)
{
  return Descent<F, G, N>::step(params, f, g, x);
}

} // namespace descent
} // namespace algo_opt