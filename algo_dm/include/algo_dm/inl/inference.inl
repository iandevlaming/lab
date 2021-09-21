#pragma once

#include <numeric>
#include <stdexcept>

namespace algo_dm
{
template <typename ItrT>
Factor product(const ItrT& begin, const ItrT& end)
{
  if (begin == end)
    throw std::invalid_argument(
        "Cannot compute the product of an empty container");

  auto product = [](const auto& lhs, const auto& rhs) { return lhs * rhs; };
  return std::accumulate(std::next(begin), end, *begin, product);
}

// template <template <typename...> typename C, typename... Ts>
// Factor product(const C<Factor, Ts...>& factors)
// {
//   if (factors.cbegin() == factors.cend())
//     throw std::invalid_argument(
//         "Cannot compute the product of an empty container");

//   auto product = [](const auto& lhs, const auto& rhs) { return lhs * rhs; };
//   return std::accumulate(
//       factors.cbegin() + 1, factors.cend(), factors.front(), product);
// }

template <template <typename...> typename C, typename... Ts>
Factor product(const C<Factor, Ts...>& factors)
{
  return product(factors.cbegin(), factors.cend());
}
} // namespace algo_dm