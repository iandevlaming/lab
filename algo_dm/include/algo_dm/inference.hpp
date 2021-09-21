#pragma once

#include <algo_dm/types.hpp>

#include <ranges>

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

Factor condition(const Factor& factor, const Assignment& assignment);

std::vector<Factor> condition(const std::vector<Factor>& factors,
                              const Assignment& evidence);

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

bool isConsistent(const Assignment& assignment, const Assignment& evidence);

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

struct GibbsSampling
{
  int num_samples;
  int num_burnin;
  int num_skip;
  std::vector<Variable::Name> ordering;
};

template <typename MethodT>
Factor infer(const MethodT& method,
             const BayesianNetwork& bn,
             const std::vector<Variable::Name>& query,
             const Assignment& evidence);

Assignment updateGibbsSample(const Assignment& assignment,
                             const BayesianNetwork& bn,
                             const Assignment& evidence,
                             const std::vector<Variable::Name>& ordering);

Assignment sampleGibbs(const Assignment& assignment,
                       const BayesianNetwork& bn,
                       const Assignment& evidence,
                       const std::vector<Variable::Name>& ordering,
                       int num_samples);
} // namespace algo_dm

#include <algo_dm/inl/inference.inl>