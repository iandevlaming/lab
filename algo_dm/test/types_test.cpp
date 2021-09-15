#include <algo_dm/alias.hpp>
#include <algo_dm/types.hpp>

#include <gtest/gtest.h>

#include <map>
#include <vector>

TEST(GetKeysTest, Simple)
{
  auto expected_keys = std::vector<std::string>();
  expected_keys.emplace_back("key_1");
  expected_keys.emplace_back("key_2");
  expected_keys.emplace_back("key_3");

  auto map = std::map<std::string, int>();
  for (const auto& key : expected_keys)
    map[key] = 1;

  auto keys = ad::getKeys(map);

  EXPECT_EQ(keys, expected_keys);
}

TEST(VariableTest, NegativeNumVals)
{
  // copy arguments
  auto name = ad::Variable::Name("name");
  auto num_vals = -1;
  EXPECT_THROW(ad::Variable(name, num_vals), std::domain_error);

  // move arguments
  EXPECT_THROW(ad::Variable("name", -1), std::domain_error);
}

TEST(VariableTest, Equality)
{
  auto name = ad::Variable::Name("name");
  auto num_vals = 1;
  auto var_copy = ad::Variable(name, num_vals);
  auto var_move = ad::Variable(ad::Variable::Name("name"), 1);

  EXPECT_EQ(var_copy, var_move);
}

TEST(VariableTest, GetNumVals)
{
  auto name = ad::Variable::Name("name");
  auto num_vals = 1;
  auto var = ad::Variable(name, num_vals);
  EXPECT_EQ(var.getNumVals(), num_vals);
}

TEST(VariableTest, GetName)
{
  auto name = ad::Variable::Name("name");
  auto num_vals = 1;
  auto var = ad::Variable(name, num_vals);
  EXPECT_EQ(var.getName(), name);
}

TEST(GetNamesTest, SimpleTest)
{
  auto expected_names = std::vector<ad::Variable::Name>();
  expected_names.emplace_back("name_1");
  expected_names.emplace_back("name_2");
  expected_names.emplace_back("name_3");

  auto variables = std::vector<ad::Variable>();
  for (const auto& name : expected_names)
    variables.emplace_back(name, 1);

  auto names = ad::getNames(variables);

  EXPECT_EQ(names, expected_names);
}

TEST(AssignmentTest, ConstructionTest)
{
  auto map = std::map<ad::Variable::Name, int>();
  map["var_1"] = 1;
  map["var_2"] = 2;
  map["var_3"] = 3;

  auto assignment = ad::Assignment(map);
  EXPECT_EQ(assignment.getKeys(), ad::getKeys(map));

  for (const auto& item : map)
    EXPECT_EQ(assignment.get(item.first).value(), item.second);
}

TEST(AssignmentTest, AssignmentTest)
{
  auto map = std::map<ad::Variable::Name, int>();
  map["var_1"] = 1;
  map["var_2"] = 2;
  map["var_3"] = 3;

  auto assignment = ad::Assignment();
  for (const auto& item : map)
    assignment.set(item.first, item.second);

  EXPECT_EQ(assignment.getKeys(), ad::getKeys(map));

  for (const auto& item : map)
    EXPECT_EQ(assignment.get(item.first).value(), item.second);
}

TEST(AssignmentTest, EqualityTest)
{
  auto map = std::map<ad::Variable::Name, int>();
  map["var_1"] = 1;
  map["var_2"] = 2;
  map["var_3"] = 3;

  auto assignment_copy = ad::Assignment(map);
  auto assignment_assign = ad::Assignment();
  for (const auto& item : map)
    assignment_assign.set(item.first, item.second);

  EXPECT_EQ(assignment_copy, assignment_assign);
}

TEST(FactorTableTest, BadConstructionTest)
{
  auto assignment_1 = ad::Assignment();
  assignment_1.set("var_1", 1);

  auto assignment_2 = ad::Assignment();
  assignment_2.set("var_2", 1);

  auto table =
      std::unordered_map<ad::Assignment, double, ad::Assignment::Hash>();
  table[assignment_1] = 0.2;
  table[assignment_2] = 0.8;

  auto bad_constructor = [&table]() { return ad::FactorTable(table); };
  EXPECT_THROW(bad_constructor(), std::invalid_argument);
}

TEST(FactorTableTest, ConstructionTest)
{
  auto assignment_1 = ad::Assignment();
  assignment_1.set("var_1", 1);
  assignment_1.set("var_2", 0);

  auto assignment_2 = ad::Assignment();
  assignment_2.set("var_1", 0);
  assignment_2.set("var_2", 1);

  auto table =
      std::unordered_map<ad::Assignment, double, ad::Assignment::Hash>();
  table[assignment_1] = 0.2;
  table[assignment_2] = 0.8;

  auto factor_table = ad::FactorTable(table);
  auto factor_table_keys = factor_table.getKeys();
  EXPECT_EQ(factor_table_keys, ad::getKeys(table));
  for (const auto& key : factor_table_keys)
    EXPECT_EQ(factor_table.get(key).value(), table[key]);
}

TEST(FactorTableTest, EmptyConstructionTest)
{
  auto empty_tables = std::vector<ad::FactorTable>();
  empty_tables.emplace_back(
      std::unordered_map<ad::Assignment, double, ad::Assignment::Hash>());
  empty_tables.emplace_back();

  for (const auto& table : empty_tables)
  {
    EXPECT_TRUE(table.getKeys().empty());
    EXPECT_TRUE(table.getAssignmentKeys().empty());
  }

  auto is_matching = [&t = empty_tables.front()](const auto& table) {
    return table == t;
  };
  EXPECT_TRUE(std::ranges::all_of(empty_tables, is_matching));
}

TEST(FactorTableTest, BadGetTest)
{
  auto assignment = ad::Assignment();
  assignment.set("var", 1);

  auto table = ad::FactorTable();
  EXPECT_FALSE(table.get(assignment).has_value());
}

TEST(FactorTableTest, BadSetTest)
{
  auto assignment_1 = ad::Assignment();
  assignment_1.set("var_1", 1);

  auto assignment_2 = ad::Assignment();
  assignment_2.set("var_2", 0);

  auto table = ad::FactorTable();
  table.set(assignment_1, 0.0);

  EXPECT_THROW(table.set(assignment_2, 0.0), std::invalid_argument);
}

TEST(FactorTableTest, SetTest)
{
  auto var_name = std::string("var");

  auto assignment_1 = ad::Assignment();
  assignment_1.set(var_name, 1);
  auto prob_1 = 0.2;

  auto assignment_2 = ad::Assignment();
  assignment_2.set(var_name, 0);
  auto prob_2 = 0.8;

  auto table = ad::FactorTable();
  table.set(assignment_1, prob_1);
  table.set(assignment_2, prob_2);

  auto assignments = std::vector<ad::Assignment>({assignment_1, assignment_2});
  EXPECT_EQ(table.getAssignmentKeys(),
            std::vector<ad::Assignment::Key>({var_name}));
  EXPECT_EQ(table.getKeys().size(), 2);

  for (const auto& assignment : table.getKeys())
    EXPECT_NE(std::ranges::find(assignments, assignment), assignments.cend());

  EXPECT_EQ(table.get(assignment_1).value(), prob_1);
  EXPECT_EQ(table.get(assignment_2).value(), prob_2);
}

TEST(FactorTableTest, LateInitializationTest)
{
  auto test_tables = std::vector<ad::FactorTable>();
  test_tables.emplace_back(
      std::unordered_map<ad::Assignment, double, ad::Assignment::Hash>());
  test_tables.emplace_back();

  auto var_name = std::string("var");
  auto assignment = ad::Assignment();
  assignment.set(var_name, 1);
  auto probability = 1.0;
  for (auto& table : test_tables)
  {
    table.set(assignment, probability);
    EXPECT_EQ(table.get(assignment).value(), probability);
    EXPECT_EQ(table.getKeys(), std::vector<ad::Assignment>({assignment}));
    EXPECT_EQ(table.getAssignmentKeys(),
              std::vector<ad::Assignment::Key>({var_name}));
  }

  auto is_matching = [&t = test_tables.front()](const auto& table) {
    return table == t;
  };
  EXPECT_TRUE(std::ranges::all_of(test_tables, is_matching));
}

TEST(FactorTableTest, NormalizeTest)
{
  auto var_name = std::string("var");

  auto assignment_1 = ad::Assignment();
  assignment_1.set(var_name, 0);
  auto assignment_2 = ad::Assignment();
  assignment_2.set(var_name, 1);

  auto prob = 0.1;
  auto table = ad::FactorTable();
  table.set(assignment_1, prob);
  table.set(assignment_2, prob);

  table.normalize();
  EXPECT_DOUBLE_EQ(table.get(assignment_1).value(), 0.5);
  EXPECT_DOUBLE_EQ(table.get(assignment_2).value(), 0.5);
}

TEST(FactorTest, ConstructionTest)
{
  auto var = ad::Variable("var", 2);
  auto vars = std::vector<ad::Variable>({var});

  auto assignment_1 = ad::Assignment();
  assignment_1.set(var.getName(), 0);
  auto assignment_2 = ad::Assignment();
  assignment_2.set(var.getName(), 1);

  auto prob = 0.1;
  auto table = ad::FactorTable();
  table.set(assignment_1, prob);
  table.set(assignment_2, prob);

  auto factor = ad::Factor(vars, table);

  EXPECT_EQ(factor.getVariables(), vars);
  EXPECT_EQ(factor.getFactorTable(), table);
}

TEST(FactorTest, NormalizeTest)
{
  auto var = ad::Variable("var", 2);

  auto assignment_1 = ad::Assignment();
  assignment_1.set(var.getName(), 0);
  auto assignment_2 = ad::Assignment();
  assignment_2.set(var.getName(), 1);

  auto prob = 0.1;
  auto table = ad::FactorTable();
  table.set(assignment_1, prob);
  table.set(assignment_2, prob);

  auto factor = ad::Factor({var}, table);
  factor.normalize();

  EXPECT_DOUBLE_EQ(factor.getFactorTable().get(assignment_1).value(), 0.5);
  EXPECT_DOUBLE_EQ(factor.getFactorTable().get(assignment_2).value(), 0.5);
}

TEST(SelectTest, SimpleTest)
{
  auto vars = std::vector<ad::Variable::Name>({"var_1", "var_2", "var_3"});
  auto assignment = ad::Assignment();
  for (auto i = 0; i < static_cast<int>(vars.size()); ++i)
    assignment.set(vars[i], i);

  auto sub_vars = std::vector<ad::Variable::Name>(vars.begin(), vars.end() - 1);
  auto sub_assignment = ad::select(assignment, sub_vars);

  EXPECT_EQ(sub_assignment.getKeys(), sub_vars);
  for (auto i = 0; i < static_cast<int>(sub_vars.size()); ++i)
    EXPECT_DOUBLE_EQ(sub_assignment.get(sub_vars[i]).value(),
                     assignment.get(sub_vars[i]).value());
}

TEST(SelectTest, BadVariableNameTest)
{
  auto assignment = ad::Assignment();
  assignment.set("var_1", 1);

  EXPECT_THROW(ad::select(assignment, {std::string("var_2")}),
               std::invalid_argument);
}

TEST(CartesianProductTest, CartesianProductSingle)
{
  auto num_vals = 4;
  auto variable = ad::Variable("variable", num_vals);
  auto variable_vec = std::vector<ad::Variable>({variable});

  auto expected_product = std::vector<std::vector<int>>();
  for (auto i = 0; i < num_vals; ++i)
    expected_product.push_back({i});

  auto product = ad::cartesianProduct(variable_vec);

  EXPECT_EQ(product, expected_product);
}

TEST(CartesianProductTest, CartesianProductDouble)
{
  auto variable_2 = ad::Variable("variable_2", 2);
  auto variable_3 = ad::Variable("variable_3", 3);
  auto variable_vec = std::vector<ad::Variable>({variable_2, variable_3});

  auto expected_product = std::vector<std::vector<int>>(
      {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}});
  auto product = ad::cartesianProduct(variable_vec);

  EXPECT_EQ(product, expected_product);
}

TEST(AssignTest, SimpleTest)
{
  auto variable_1 = ad::Variable("variable_1", 2);
  auto variable_2 = ad::Variable("variable_2", 2);
  auto variables = std::vector<ad::Variable>({variable_1, variable_2});

  auto products = ad::cartesianProduct(variables);
  auto expected_assignments = std::vector<ad::Assignment>();
  for (const auto& product : products)
  {
    auto assignment = ad::Assignment();
    assignment.set(variable_1.getName(), product.front());
    assignment.set(variable_2.getName(), product.back());
    expected_assignments.push_back(assignment);
  }

  auto assignments = ad::assign(variables);

  EXPECT_EQ(assignments.size(), expected_assignments.size());
  for (const auto& assignment : expected_assignments)
    EXPECT_NE(std::ranges::find(assignments, assignment), assignments.cend());
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}