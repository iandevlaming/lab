#pragma once

#include <Eigen/Core>

#include <type_traits>

namespace algo_opt
{
// type alias templates
template <int N>
using Vectord = Eigen::Matrix<double, N, 1>;

// function alias templates
template <int N>
using multivar_f_t = std::function<double(Vectord<N>)>;
template <int N>
using multivar_grad_t = std::function<Vectord<N>(Vectord<N>)>;

// is_invocable_r alias templates
template <typename F, int N>
using is_multivar_f = std::is_invocable_r<double, F, Vectord<N>>;
template <typename G, int N>
using is_multivar_grad = std::is_invocable_r<Vectord<N>, G, Vectord<N>>;

// enable_if alias_templates
template <typename F, int N>
using enable_if_multivar_f = std::enable_if<is_multivar_f<F, N>::value>;
template <typename G, int N>
using enable_if_multivar_grad = std::enable_if<is_multivar_grad<G, N>::value>;

// invocable alias specializations
using binary_f_t = multivar_f_t<2>;
using binary_grad_t = multivar_grad_t<2>;
using binary_level_t =
    std::function<std::optional<std::array<double, 2>>(double, double)>;
} // namespace algo_opt