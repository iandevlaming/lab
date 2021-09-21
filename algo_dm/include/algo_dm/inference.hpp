#pragma once

#include <algo_dm/types.hpp>

#include <ranges>
#include <unordered_map>

namespace algo_dm
{
template <typename ItrT>
Factor product(const ItrT& begin, const ItrT& end);

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

Factor condition(const Factor& factor, const Assignment& assignment);

Assignment sample(const Factor& factor);

class SampleGenerator
{
public:
  SampleGenerator(const BayesianNetwork& bn);
  Assignment sample() const;

private:
  const BayesianNetwork& bn_;
  std::vector<Variable::Name> order_;
};

Assignment sample(const BayesianNetwork& bn);

bool isConsistent(
    const Assignment& assignment,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence);

Factor computeBlanket(const BayesianNetwork& bn,
                      const Assignment& assignment,
                      const Variable::Name& variable);

struct ExactInference
{
};

struct VariableElimination
{
  std::vector<Variable::Name> ordering;
};

struct DirectSampling
{
  int num_samples;
};

struct LikelihoodWeightedSampling
{
  int num_samples;
};

template <typename MethodT>
Factor
infer(const MethodT& method,
      const BayesianNetwork& bn,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence);
} // namespace algo_dm

#include <algo_dm/inl/inference.inl>