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
} // namespace algo_dm