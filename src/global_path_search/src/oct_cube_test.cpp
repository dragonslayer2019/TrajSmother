// created by yuqing.wu on 20/06/24

#include "../include/oct_cube.h"

#include <gtest/gtest.h>

#include <iostream>
#include <memory>
#include <vector>
#include <queue>

namespace global_path_search {

class OctCubeTest : public ::testing::Test {
 public:
  OctCubeTest() {
    auto FILE = freopen("oct_cube_test_stdout", "w", stdout);
  }
  ~OctCubeTest() {}

 private:
  void TestBody() override { std::cout << "Test body!\n"; }

 public:
  bool Run() {
    std::cout << "Run() fundtion. \n";
    TestBody();
    return true;
  }

 public:
    std::shared_ptr<OctCube> oct_cube_;
    size_t GetCubeLevel() { return oct_cube_->cube_level; }
    CubeType GetCubeType() { return oct_cube_->cube_type; }
};

/*
TEST_F(OctCubeTest, first_gtest) {
    EXPECT_NEAR(1.0, 1.0, 1e-6);
}

TEST_F(OctCubeTest, second_gtest) {
    int x = 1;
    EXPECT_EQ(x, 1);
}
*/

TEST_F(OctCubeTest, build_cube_test) {

    std::shared_ptr<OctCubeTest> p_test = std::make_shared<OctCubeTest>();
    printf("check test\n");
    OcTree::OctNode* oct_node = new OcTree::OctNode();
    p_test->oct_cube_ = std::make_shared<OctCube>(oct_node, 0, 0, 0);
    printf("oct_cube build\n");

    EXPECT_EQ(p_test->GetCubeLevel(), 0);
    EXPECT_EQ(p_test->GetCubeType(), CubeType::CUBE_TWO);
}

} // namespace global_path_search
