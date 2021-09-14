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
  auto bad_copy_constructor = []() {
    auto name = ad::Variable::Name("name");
    auto num_vals = -1;
    return ad::Variable(name, num_vals);
  };
  EXPECT_THROW(bad_copy_constructor(), std::domain_error);

  auto bad_move_constructor = []() { return ad::Variable("name", -1); };
  EXPECT_THROW(bad_move_constructor(), std::domain_error);
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
    EXPECT_EQ(assignment[item.first], item.second);
}

TEST(AssignmentTest, AssignmentTest)
{
  auto map = std::map<ad::Variable::Name, int>();
  map["var_1"] = 1;
  map["var_2"] = 2;
  map["var_3"] = 3;

  auto assignment = ad::Assignment();
  for (const auto& item : map)
    assignment[item.first] = item.second;

  EXPECT_EQ(assignment.getKeys(), ad::getKeys(map));

  for (const auto& item : map)
  {
    EXPECT_EQ(assignment[item.first], item.second);
    EXPECT_EQ(assignment.at(item.first), item.second);
  }
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
    assignment_assign[item.first] = item.second;

  EXPECT_EQ(assignment_copy, assignment_assign);
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
  auto variable_2 = ad::Variable("variable", 2);
  auto variable_3 = ad::Variable("variable", 3);
  auto variable_vec = std::vector<ad::Variable>({variable_2, variable_3});

  auto expected_product = std::vector<std::vector<int>>(
      {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}});
  auto product = ad::cartesianProduct(variable_vec);

  EXPECT_EQ(product, expected_product);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}