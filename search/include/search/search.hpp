#pragma once

#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace search
{
template <typename ReturnT, typename FunctionT, typename... Args>
using has_signature =
    std::enable_if_t<std::is_invocable_r_v<ReturnT, FunctionT, Args...>>;

namespace search_impl
{
template <typename X>
void visit(std::unordered_map<X, std::pair<X, int>> &p, X parent, X child)
{
  auto count = p.contains(parent) ? p[parent].second : 0;
  ++count;
  p[child] = std::pair(parent, count);
}

template <typename X>
bool isVisited(std::unordered_map<X, std::pair<X, int>> const &p, X x)
{
  return p.contains(x);
}
} // namespace search_impl

/**
 * @brief
 *
 * @param[in]
 *
 * @tparam
 *
 * @returns
 */
template <typename Q,
          typename X,
          typename PushF,
          typename PopF,
          typename IsEmptyF,
          typename IsGoalF,
          typename NeighborsF,
          typename P = std::unordered_map<X, std::pair<X, int>>,
          typename VisitF = decltype(search_impl::visit<X>),
          typename IsVisitedF = decltype(search_impl::isVisited<X>),
          typename = has_signature<void, PushF, Q &, X>,
          typename = has_signature<X, PopF, Q &>,
          typename = has_signature<bool, IsEmptyF, Q const &>,
          typename = has_signature<bool, IsGoalF, X>,
          typename = has_signature<std::vector<X>, NeighborsF, X>,
          typename = has_signature<void, VisitF, P &, X, X>,
          typename = has_signature<bool, IsVisitedF, P const &, X>>
std::pair<X, P> search(Q q,
                       P p,
                       X start,
                       PushF push,
                       PopF pop,
                       IsEmptyF is_empty,
                       IsGoalF is_goal,
                       NeighborsF neighbors,
                       VisitF visit,
                       IsVisitedF is_visited)
{
  push(q, start);
  visit(p, start, start);

  while (!is_empty(q))
  {
    X x = pop(q);
    if (is_goal(x))
      return {x, p};
    for (auto const &xn : neighbors(x))
    {
      if (!is_visited(p, xn))
      {
        push(q, xn);
        visit(p, x, xn);
      }
      // else
      // {
      //   resolve(q, v, xn);
      // }
    }
  }

  return {start, p};
}

template <typename X>
std::vector<X> trace(std::unordered_map<X, std::pair<X, int>> const &p, X goal)
{
  auto n_x = p.at(goal).second;
  auto path = std::vector<X>(n_x);
  path.back() = goal;
  for (int i = n_x - 2; i >= 0; --i)
    path[i] = p.at(path[i + 1]).first;

  return path;
}
} // namespace search