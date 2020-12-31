#pragma once

#include <Eigen/Core>

#include <type_traits>

namespace algo_opt
{
// type alias templates
template <int N>
using Vectord = Eigen::Matrix<double, N, 1>;
template <int N, int M = N>
using Matrixd = Eigen::Matrix<double, N, M>;

// function alias templates
template <int N>
using multivar_f_t = std::function<double(Vectord<N>)>;
template <int N>
using multivar_grad_t = std::function<Vectord<N>(Vectord<N>)>;
template <int N>
using multivar_hess_t = std::function<Matrixd<N>(Vectord<N>)>;

// is_invocable_r alias templates
template <typename F, int N>
using is_multivar_f = std::is_invocable_r<double, F, Vectord<N>>;
template <typename G, int N>
using is_multivar_grad = std::is_invocable_r<Vectord<N>, G, Vectord<N>>;
template <typename H, int N>
using is_multivar_hess = std::is_invocable_r<Matrixd<N>, H, Vectord<N>>;
template <typename F>
using is_unary_f = std::is_invocable_r<double, F, double>;

// enable_if alias_templates
template <typename F, int N>
using enable_if_multivar_f = std::enable_if<is_multivar_f<F, N>::value>;
template <typename G, int N>
using enable_if_multivar_grad = std::enable_if<is_multivar_grad<G, N>::value>;
template <typename H, int N>
using enable_if_multivar_hess = std::enable_if<is_multivar_hess<H, N>::value>;

template <typename F>
using enable_if_unary_f = std::enable_if<is_unary_f<F>::value>;

// invocable alias specializations
using unary_f_t = std::function<double(double)>;

using binary_f_t = multivar_f_t<2>;
using binary_grad_t = multivar_grad_t<2>;
using binary_hess_t = multivar_hess_t<2>;
using binary_level_t =
    std::function<std::optional<std::array<double, 2>>(double, double)>;

template <int N>
using bracket_t = std::array<double, N>;

struct Point2d
{
  double x;
  double y;
};

} // namespace algo_opt

namespace ao = algo_opt;