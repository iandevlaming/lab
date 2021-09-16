#pragma once

#include <algo_dm/types.hpp>

#include <ranges>

namespace algo_dm
{
template <template <typename...> typename C, typename... Ts>
C<Ts...> setDiff(const C<Ts...>& super_set, const C<Ts...>& disjoint_set)
{
  auto is_not_in = [](const auto& set_1) {
    return [&set_1](const auto& elem_2) {
      return std::ranges::find(set_1, elem_2) == set_1.cend();
    };
  };

  auto sub_set_view = std::views::filter(super_set, is_not_in(disjoint_set));
  auto sub_set = C<Ts...>(sub_set_view.begin(), sub_set_view.end());

  return sub_set;
}

Factor operator*(const Factor& lhs, const Factor& rhs);

Factor marginalize(const Factor& factor, const Variable::Name name);
} // namespace algo_dm