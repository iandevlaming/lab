#include <algo_dm/inference.hpp>

#include <algorithm>
#include <cstdlib>
#include <numeric>

namespace algo_dm
{
Factor marginalize(const Factor& factor, const Variable::Name name)
{
  auto marginalized_table = FactorTable();
  const auto& table = factor.getFactorTable();
  for (const auto& assignment : table.getAssignments())
  {
    auto pruned_assignment = assignment;
    pruned_assignment.erase(name);

    auto new_p = table.get(assignment);
    if (marginalized_table.contains(pruned_assignment))
      new_p += marginalized_table.get(pruned_assignment);

    marginalized_table.set(pruned_assignment, new_p);
  }

  auto reduced_vars = erase(factor.getVariables(), name);
  return Factor(reduced_vars, marginalized_table);
}

Factor marginalize(const Factor& factor,
                   const std::vector<Variable::Name> names)
{
  if (names.empty())
    return factor;

  auto apply_margin = [](const auto& f, const auto& name) {
    return marginalize(f, name);
  };
  return std::accumulate(names.cbegin(), names.cend(), factor, apply_margin);
}

Factor condition(const Factor& factor,
                 const Variable::Name name,
                 const Assignment::Value& value)
{
  if (!isInScope(factor, name))
    return factor;

  auto conditioned_table = FactorTable();
  const auto& table = factor.getFactorTable();
  for (const auto& assignment : table.getAssignments())
  {
    if (assignment.get(name) == value)
    {
      auto conditioned_assignment = assignment;
      conditioned_assignment.erase(name);
      conditioned_table.set(conditioned_assignment, table.get(assignment));
    }
  }

  auto reduced_vars = erase(factor.getVariables(), name);
  return Factor(reduced_vars, conditioned_table);
}

Factor
condition(const Factor& factor,
          const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  auto condition_fold = [](const auto& f, const auto& e) {
    return condition(f, e.first, e.second);
  };
  return std::accumulate(
      evidence.cbegin(), evidence.cend(), factor, condition_fold);
}

Factor condition(const Factor& factor, const Assignment& assignment)
{
  const auto& vars = assignment.getVariableNames();
  auto condition_fold = [&assignment](const auto& f, const auto& v) {
    return condition(f, v, assignment.get(v));
  };
  return std::accumulate(vars.cbegin(), vars.cend(), factor, condition_fold);
}

std::vector<Factor>
condition(const std::vector<Factor>& factors,
          const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  auto condition_on = [](const auto& e) {
    return [&e](const auto& factor) { return condition(factor, e); };
  };
  auto condition_view =
      factors | std::views::transform(condition_on((evidence)));
  return std::vector(condition_view.begin(), condition_view.end());
}

Assignment sample(const Factor& factor)
{
  const auto& table = factor.getFactorTable();
  const auto& assignments = table.getAssignments();

  auto get_prob = [&table](const auto& a) { return table.get(a); };
  auto prob_view = assignments | std::views::transform(get_prob);

  const auto w = std::accumulate(prob_view.begin(), prob_view.end(), 0.0);
  auto total_prob = 0.0;
  auto p = static_cast<double>(std::rand()) / RAND_MAX;

  for (const auto& assignment : assignments)
  {
    total_prob += table.get(assignment) / w;
    if (total_prob >= p)
      return assignment;
  }
  return assignments.back(); // should never happen
}

SampleGenerator::SampleGenerator(const BayesianNetwork& bn)
    : bn_(bn), order_(topoSort(bn.getGraph()))
{
}

Assignment SampleGenerator::sample() const
{
  auto assignment = Assignment();
  for (const auto& node : order_)
  {
    const auto& factor = bn_.getFactor(node);
    auto conditioned_factor = condition(factor, assignment);
    auto random_assignment = algo_dm::sample(conditioned_factor);
    auto value = random_assignment.get(node);
    assignment.set(node, value);
  }
  return assignment;
}

Assignment sample(const BayesianNetwork& bn)
{
  auto generator = SampleGenerator(bn);
  return generator.sample();
}

bool isConsistent(
    const Assignment& assignment,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  auto vars = evidence | std::views::keys;

  auto get_e = [& e = evidence](const auto& var) { return e.at(var); };
  auto e_vals = vars | std::views::transform(get_e);

  auto get_a = [& a = assignment](const auto& var) {
    return a.contains(var) ? a.get(var) : -1;
  };
  auto a_vals = vars | std::views::transform(get_a);

  return std::ranges::equal(e_vals, a_vals);
}

template <>
Factor infer<ExactInference>(
    const ExactInference&,
    const BayesianNetwork& bn,
    const std::vector<Variable::Name>& query,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  const auto& factors = bn.getFactors();

  auto inference = product(factors);
  inference = condition(inference, evidence);

  auto non_query = setDiff(getNames(inference.getVariables()), query);
  inference = marginalize(inference, non_query);

  inference.normalize();
  return inference;
}

template <>
Factor infer<VariableElimination>(
    const VariableElimination& method,
    const BayesianNetwork& bn,
    const std::vector<Variable::Name>& query,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  const auto& factors = bn.getFactors();

  auto reduced_factors = condition(factors, evidence);

  auto all_names = setUnion(method.ordering, getNames(reduced_factors));

  auto depends_on = [](const auto& n) {
    return [&n](const auto& f) { return isInScope(f, n); };
  };
  for (const auto& name : all_names)
  {
    if (!contains(query, name))
    {
      auto dependent_factors_view =
          reduced_factors | std::views::filter(depends_on(name));
      auto dependent_factors = std::vector(dependent_factors_view.begin(),
                                           dependent_factors_view.end());

      if (!dependent_factors.empty())
      {
        auto reduced_factor = product(dependent_factors);
        reduced_factor = marginalize(reduced_factor, name);

        auto remove_itrs =
            std::ranges::remove_if(reduced_factors, depends_on(name));
        reduced_factors.erase(remove_itrs.begin(), remove_itrs.end());
        reduced_factors.push_back(reduced_factor);
      }
    }
  }

  auto reduced_factor = product(reduced_factors);
  reduced_factor.normalize();
  return reduced_factor;
}

template <>
Factor infer<DirectSampling>(
    const DirectSampling& method,
    const BayesianNetwork& bn,
    const std::vector<Variable::Name>& query,
    const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  auto table = FactorTable();
  auto generator = SampleGenerator(bn);

  for (auto i = 0; i < method.num_samples; ++i)
  {
    auto assignment = generator.sample();
    if (isConsistent(assignment, evidence))
    {
      auto sub_assignment = select(assignment, query);
      auto p = 1.0;
      if (table.contains(sub_assignment))
        p += table.get(sub_assignment);
      table.set(sub_assignment, p);
    }
  }
  table.normalize();

  auto vars = select(bn.getVariables(), query);
  return Factor(vars, table);
}

} // namespace algo_dm