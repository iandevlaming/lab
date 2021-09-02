#pragma once

#include <search/search.hpp>

#include <queue>
#include <type_traits>
#include <utility>

namespace search
{
namespace bfs_impl
{
template <typename X>
void push(std::queue<X> &q, X x)
{
  q.push(x);
}

template <typename X>
X pop(std::queue<X> &q)
{
  X x = q.front();
  q.pop();
  return x;
}

template <typename X>
bool isEmpty(std::queue<X> const &q)
{
  return q.empty();
}
} // namespace bfs_impl

/**
 * @brief
 *
 * @param[in]
 *
 * @tparam
 *
 * @returns
 */
template <typename X,
          typename IsGoalF,
          typename NeighborsF,
          typename Q = std::queue<X>,
          typename P = std::unordered_map<X, std::pair<X, int>>,
          typename PushF = decltype(bfs_impl::push<X>),
          typename PopF = decltype(bfs_impl::pop<X>),
          typename IsEmptyF = decltype(bfs_impl::isEmpty<X>),
          typename VisitF = decltype(search_impl::visit<X>),
          typename IsVisitedF = decltype(search_impl::isVisited<X>),
          typename = has_signature<bool, IsGoalF, X>,
          typename = has_signature<std::vector<X>, NeighborsF, X>,
          typename = has_signature<void, PushF, Q &, X>,
          typename = has_signature<X, PopF, Q &>,
          typename = has_signature<bool, IsEmptyF, Q const &>,
          typename = has_signature<void, VisitF, P &, X, X>,
          typename = has_signature<bool, IsVisitedF, P const &, X>>
std::vector<X>
breadth_first_search(X start,
                     IsGoalF is_goal,
                     NeighborsF neighbors,
                     Q q = Q(),
                     P p = P(),
                     PushF push = bfs_impl::push<X>,
                     PopF pop = bfs_impl::pop<X>,
                     IsEmptyF is_empty = bfs_impl::isEmpty<X>,
                     VisitF visit = search_impl::visit<X>,
                     IsVisitedF is_visited = search_impl::isVisited<X>)
{
  auto [goal, p_res] = search(
      q, p, start, push, pop, is_empty, is_goal, neighbors, visit, is_visited);

  if (!is_goal(goal))
    return {};

  return trace(p_res, goal);
}
} // namespace search