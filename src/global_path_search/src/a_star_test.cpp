// created by yuqing.wu on 27/06/24

#include "a_star.h"
#include "vertice.h"

#include <gtest/gtest.h>

#include <iostream>
#include <memory>
#include <vector>
#include <queue>

namespace global_path_search {

class AstarTest : public::testing::Test {
 public:
    AstarTest() {}
    ~AstarTest() {}

 private:
    void TestBody() override { std::cout << "Test body!\n"; }

 public:
    bool Run() {
        std::cout << "Run() fundtion. \n";
        TestBody();
        return true;
    }

 public:
    std::shared_ptr<Astar> a_star_;

};

TEST_F(AstarTest, astar_essential) {
    auto FILE = freopen("astar_stdout", "w", stdout);
    std::shared_ptr<AstarTest> p_test = std::make_shared<AstarTest>();
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    A->MainUpdate(A->GetRoot(), {}, {{4.26, 1.34, 3.12}});
    // auto searcher = std::make_shared<Searcher>(A);
    std::cout << "make test" << std::endl;
    // p_test->a_star_ = std::make_shared<Astar>(nullptr, nullptr);
    std::cout << "check place" << std::endl;
    
    // EXPECT_EQ(0, 0);
}

} // namespace global_path_search
