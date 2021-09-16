#include <algo_dm/inference.hpp>

#include <algorithm>

namespace algo_dm
{
Factor operator*(const Factor& lhs, const Factor& rhs)
{
  const auto& lhs_vars = lhs.getVariables();
  const auto& rhs_vars = rhs.getVariables();

  const auto& lhs_table = lhs.getFactorTable();
  const auto& rhs_table = rhs.getFactorTable();

  const auto rhs_names = getNames(rhs_vars);

  const auto rhs_only_vars = setDiff(rhs_vars, lhs_vars);
  const auto rhs_only_assignments = assign(rhs_only_vars);

  auto merged_table = FactorTable();
  for (const auto& lhs_assignment : lhs_table.getAssignments())
  {
    for (const auto& rhs_only_assignment : rhs_only_assignments)
    {
      auto merged_assignment = rhs_only_assignment;
      for (const auto& lhs_var_name : lhs_assignment.getVariableNames())
        merged_assignment.set(lhs_var_name,
                              lhs_assignment.get(lhs_var_name).value());

      const auto rhs_assignment = select(merged_assignment, rhs_names);
      const auto merged_probability = lhs_table.get(lhs_assignment).value() *
                                      rhs_table.get(rhs_assignment).value();

      merged_table.set(merged_assignment, merged_probability);
    }
  }

  auto merged_vars = lhs_vars;
  std::ranges::copy(rhs_only_vars, std::back_inserter(merged_vars));

  return Factor(merged_vars, merged_table);
}

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
} // namespace algo_dm