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

Factor condition(const Factor& factor, const Assignment& assignment)
{
  const auto& vars = assignment.getVariableNames();
  auto condition_fold = [&assignment](const auto& f, const auto& v) {
    return condition(f, v, assignment.get(v));
  };
  return std::accumulate(vars.cbegin(), vars.cend(), factor, condition_fold);
}

std::vector<Factor> condition(const std::vector<Factor>& factors,
                              const Assignment& evidence)
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

bool isConsistent(const Assignment& assignment, const Assignment& evidence)
{
  auto vars = evidence.getVariableNames();

  auto get_from = [](const auto& c) {
    return [&c](const auto& v) { return c.contains(v) ? c.get(v) : -1; };
  };
  auto e_vals = vars | std::views::transform(get_from(evidence));
  auto a_vals = vars | std::views::transform(get_from(assignment));

  return std::ranges::equal(e_vals, a_vals);
}

Factor computeBlanket(const BayesianNetwork& bn,
                      const Assignment& assignment,
                      const Variable::Name& variable)
{
  auto sub_assignment = assignment;
  sub_assignment.erase(variable);

  auto depends_on = [](const auto& v) {
    return [&v](const auto& f) { return isInScope(f, v); };
  };
  auto condition_on = [](const auto& a) {
    return [&a](const auto& f) { return condition(f, a); };
  };

  const auto& factors = bn.getFactors();
  auto condition_view = factors | std::views::filter(depends_on(variable)) |
                        std::views::transform(condition_on(sub_assignment));

  auto blanket = product(condition_view.begin(), condition_view.end());
  blanket.normalize();
  return blanket;
}

template <>
Factor infer<ExactInference>(const ExactInference&,
                             const BayesianNetwork& bn,
                             const std::vector<Variable::Name>& query,
                             const Assignment& evidence)
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
Factor infer<VariableElimination>(const VariableElimination& method,
                                  const BayesianNetwork& bn,
                                  const std::vector<Variable::Name>& query,
                                  const Assignment& evidence)
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
Factor infer<DirectSampling>(const DirectSampling& method,
                             const BayesianNetwork& bn,
                             const std::vector<Variable::Name>& query,
                             const Assignment& evidence)
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

template <>
Factor
infer<LikelihoodWeightedSampling>(const LikelihoodWeightedSampling& method,
                                  const BayesianNetwork& bn,
                                  const std::vector<Variable::Name>& query,
                                  const Assignment& evidence)
{
  const auto& factors = bn.getFactors();
  const auto ordering = topoSort(bn.getGraph());

  auto table = FactorTable();
  for (auto i = 0; i < method.num_samples; ++i)
  {
    auto assignment = Assignment();
    auto weight = 1.0;

    for (const auto& var : ordering)
    {
      const auto& factor = bn.getFactor(var);
      if (evidence.contains(var))
      {
        assignment.set(var, evidence.get(var));
        const auto var_names = getNames(factor.getVariables());
        const auto sub_assignment = select(assignment, var_names);
        const auto& f_table = factor.getFactorTable();
        weight *= f_table.get(sub_assignment);
      }
      else
      {
        auto conditioned_factor = condition(factor, assignment);
        auto sample_assignment = sample(conditioned_factor);
        assignment.set(var, sample_assignment.get(var));
      }
    }

    auto entry_assignment = select(assignment, query);
    if (table.contains(entry_assignment))
      weight += table.get(entry_assignment);
    table.set(entry_assignment, weight);
  }

  table.normalize();
  auto vars = select(bn.getVariables(), query);
  return Factor(vars, table);
}

Assignment updateGibbsSample(const Assignment& assignment,
                             const BayesianNetwork& bn,
                             const Assignment& evidence,
                             const std::vector<Variable::Name>& ordering)
{
  auto updated_assignment = assignment;
  for (const auto& var : ordering)
  {
    if (!evidence.contains(var))
    {
      const auto blanket = computeBlanket(bn, updated_assignment, var);
      const auto random_assignment = sample(blanket);
      updated_assignment.set(var, random_assignment.get(var));
    }
  }

  return updated_assignment;
}

Assignment sampleGibbs(const Assignment& assignment,
                       const BayesianNetwork& bn,
                       const Assignment& evidence,
                       const std::vector<Variable::Name>& ordering,
                       int num_samples)
{
  auto updated_assignment = assignment;
  for (auto i = 0; i < num_samples; ++i)
    updated_assignment =
        updateGibbsSample(updated_assignment, bn, evidence, ordering);
  return updated_assignment;
}

template <>
Factor infer<GibbsSampling>(const GibbsSampling& method,
                            const BayesianNetwork& bn,
                            const std::vector<Variable::Name>& query,
                            const Assignment& evidence)
{
  auto table = FactorTable();
  auto ordering = setUnion(method.ordering, getNames(bn.getVariables()));

  auto random_sample = sample(bn);
  for (const auto& e : evidence)
    random_sample.set(e.first, e.second);

  random_sample =
      sampleGibbs(random_sample, bn, evidence, ordering, method.num_burnin);
  for (auto i = 0; i < method.num_samples; ++i)
  {
    random_sample =
        sampleGibbs(random_sample, bn, evidence, ordering, method.num_skip);
    const auto& random_assignment = select(random_sample, query);

    auto n_observations = 1.0;
    if (table.contains(random_assignment))
      n_observations += table.get(random_assignment);
    table.set(random_assignment, n_observations);
  }
  table.normalize();

  const auto vars = select(bn.getVariables(), query);
  return Factor(vars, table);
}
} // namespace algo_dm