#include <algo_dm/inference.hpp>

#include <algorithm>
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

    auto existing_p = marginalized_table.get(pruned_assignment);
    auto new_p = existing_p.has_value() ? existing_p.value() : 0.0;
    new_p += table.get(assignment).value();

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
    if (assignment.get(name).value() == value)
    {
      auto conditioned_assignment = assignment;
      conditioned_assignment.erase(name);
      conditioned_table.set(conditioned_assignment,
                            table.get(assignment).value());
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

Factor
infer(const BayesianNetwork& bn,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  return infer(bn.getFactors(), query, evidence);
}

Factor
infer(const std::vector<Factor>& factors,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence)
{
  auto inference = product(factors);
  inference = condition(inference, evidence);

  auto non_query = setDiff(getNames(inference.getVariables()), query);
  inference = marginalize(inference, non_query);

  inference.normalize();
  return inference;
}

Factor
infer(const std::vector<Factor>& factors,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence,
      const std::vector<Variable::Name>& ordering)
{
  auto reduced_factors = condition(factors, evidence);

  auto all_names = setUnion(ordering, getNames(reduced_factors));

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

Factor
infer(const BayesianNetwork& bn,
      const std::vector<Variable::Name>& query,
      const std::unordered_map<Variable::Name, Assignment::Value>& evidence,
      const std::vector<Variable::Name>& ordering)
{
  return infer(bn.getFactors(), query, evidence, ordering);
}
} // namespace algo_dm