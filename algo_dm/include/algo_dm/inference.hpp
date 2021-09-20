#pragma once

#include <algo_dm/types.hpp>

#include <ranges>
#include <unordered_map>

namespace algo_dm
{
template <template <typename...> typename C, typename... Ts>
Factor product(const C<Factor, Ts...>& factors);

Factor marginalize(const Factor& factor, const Variable::Name name);

Factor marginalize(const Factor& factor,
                   const std::vector<Variable::Name> names);

Factor condition(const Factor& factor,
                 const Variable::Name name,
                 const Assignment::Value& value);

Factor condition(
    const Factor& factor,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence);

std::vector<Factor> condition(
    const std::vector<Factor>& factors,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence);

Factor
infer(const BayesianNetwork& bn,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence);

Factor
infer(const std::vector<Factor>& factors,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence);

Factor
infer(const BayesianNetwork& bn,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence,
      const std::vector<Variable::Name>& ordering);

Factor
infer(const std::vector<Factor>& factors,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence,
      const std::vector<Variable::Name>& ordering);
} // namespace algo_dm

#include <algo_dm/inl/inference.inl>