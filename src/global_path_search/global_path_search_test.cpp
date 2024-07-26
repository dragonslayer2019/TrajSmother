// created by yuqing.wu on 21/06/24
#include "global_path_search.h"
#include "zyb_octree.h"
#include "oct_cube.h"
#include "a_star.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <chrono>

#include <iostream>
#include <memory>
#include <vector>
#include <queue>

namespace global_path_search {

bool test_flag = false;

class GlobalPathSearchTest : public ::testing::Test {
 public:
    GlobalPathSearchTest() {
    }
    ~GlobalPathSearchTest() {}

 private:
    void TestBody() override { std::cout << "Test body!\n"; }

 public:
    bool Run() {
        std::cout << "Run() fundtion. \n";
        TestBody();
        return true;
    }

 public:
    std::shared_ptr<Searcher> searcher_;

 public:
    void AddOctNode(OcTree::OctNode* oct_node, size_t cur_level = 0, size_t rela_to_node = 0, size_t cur_index = 0) {
        searcher_->AddOctNode(oct_node, cur_level, rela_to_node, cur_index);
        return;
    }

    NodeId GetIndex(Point3D pt, int level) {
        std::cout << "test get index" << std::endl;
        return searcher_->GetIndex(pt, level);
    }

    void HalfSolve() {
        std::cout << "wh ??" << std::endl;
        searcher_->HalfSolve(Point3D(1., 1., 1.), Point3D(127., 127., 127.));
        return;
    }

    Point3D BivariateGaussian(double x, double y, double sigma_x, double sigma_y, double mu_x, double mu_y, double norm = 62.) {
        // double norm = 62.;
        double gauss = norm * std::exp(-(std::pow(x - mu_x, 2.) / (2. * std::pow(sigma_x, 2.)) + std::pow(y - mu_y, 2.) / (2. * std::pow(sigma_y, 2.))));
        return Point3D(x, y, std::max(gauss - 12., 0.));
    }

    std::vector<Point3D> DrawMountain(double x, double y, double len, double dense, double x_ra = 1., double y_ra = 1., double norm = 62.) {
        std::vector<Point3D> res;
        double eps = 1e-5;
        double x_min = std::max(x - len / 2., eps);
        double x_max = x + len / 2.;
        double y_min = std::max(y - len / 2., eps);
        double y_max = y + len / 2.;
        for (double cur_x = x_min; cur_x < x_max + eps; cur_x += dense) {
            for (double cur_y = y_min; cur_y < y_max + eps; cur_y += dense) {
                res.push_back(BivariateGaussian(cur_x, cur_y, 32., 32., cur_x + x_ra * (cur_x - x), cur_y + y_ra * (cur_y - y), norm));
            }
        }
        for (auto pt: res) {
            std::cout << "mountain pt: " << pt.ToString() << std::endl;
        }
        return res;
    }

    std::vector<Point3D> DrawObsCube(Point3D centre_pt, double len_x, double len_y, double len_z, double dense) {
        double eps = 1e-5;
        std::vector<Point3D> res;
        double x_min = centre_pt.x - len_x / 2.;
        double x_max = centre_pt.x + len_x / 2.;
        double y_min = centre_pt.y - len_y / 2.;
        double y_max = centre_pt.y + len_y / 2.;
        double z_min = centre_pt.z - len_z / 2.;
        double z_max = centre_pt.z + len_z / 2.;
        for (double x = x_min; x < x_max + eps; x += dense) {
            for (double y = y_min; y < y_max + eps; y += dense) {
                res.push_back(Point3D(x, y, z_min));
                res.push_back(Point3D(x, y, z_max));
            }
        }
        for (double x = x_min; x < x_max + eps; x += dense) {
            for (double z = z_min; z < z_max + eps; z += dense) {
                res.push_back(Point3D(x, y_min, z));
                res.push_back(Point3D(x, y_max, z));
            }
        }
        for (double y = y_min; y < y_max + eps; y += dense) {
            for (double z = z_min; z < z_max + eps; z += dense) {
                res.push_back(Point3D(x_min, y, z));
                res.push_back(Point3D(x_max, y, z));
            }
        }
        return res;
    }

    std::vector<Point3D> DrawWall(Point3D centre_pt, double len_x, double len_y, double len_z, size_t dir, double dense) {
        double eps = 1e-5;
        std::vector<Point3D> res;
        double x_min = centre_pt.x - len_x / 2.;
        double x_max = centre_pt.x + len_x / 2.;
        double y_min = centre_pt.y - len_y / 2.;
        double y_max = centre_pt.y + len_y / 2.;
        double z_min = centre_pt.z - len_z / 2.;
        double z_max = centre_pt.z + len_z / 2.;
        if(dir == 0) {
            for (double y = y_min; y < y_max + eps; y += dense) {
                for (double z = z_min; z < z_max + eps; z += dense) {
                    res.push_back(Point3D(x_min, y, z));
                    res.push_back(Point3D(x_max, y, z));
                }
            }
            return res;
        } else if (dir == 1) {
            for (double x = x_min; x < x_max + eps; x += dense) {
                for (double z = z_min; z < z_max + eps; z += dense) {
                    res.push_back(Point3D(x, y_min, z));
                    res.push_back(Point3D(x, y_max, z));
                }
            }
            return res;
        } else {
            for (double x = x_min; x < x_max + eps; x += dense) {
                for (double y = y_min; y < y_max + eps; y += dense) {
                    res.push_back(Point3D(x, y, z_min));
                    res.push_back(Point3D(x, y, z_max));
                }
            }
            return res;
        }
        return res;
    }

/*
    void CubesBuild(std::shared_ptr<OcTree> octree) {
        std::cout << "cubes build" << std::endl;
        std::cout << "res_path size: " << searcher_->res_path_.size() << std::endl;
        searcher_->PrintCubePoolSize();
        auto all_leaf_nodes = octree->GetAllLeafNodes();
        std::cout << "all leaf nodes get" << std::endl;
        searcher_->PrintCubePoolSize();
        for (auto oct_node_pair : all_leaf_nodes) {
            std::cout << "get pair" << std::endl;
            auto oct_node = oct_node_pair.first;
            bool is_leaf = oct_node_pair.second;
            if (!is_leaf) {
                std::cout << "not leaf" << std::endl;
                searcher_->id2cube_[oct_node->index] = nullptr;
            }
            searcher_->PrintCubePoolSize();
            searcher_->AddOctNode(oct_node, oct_node->level, 0, oct_node->index);
        }
    }
*/

};

TEST_F(GlobalPathSearchTest, first_gtest) {
    auto FILE = freopen("stdout1", "w", stdout);
    std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
    // OcTree A(10);
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    p_test->searcher_ = std::make_shared<Searcher> (A);
    // printf("GlobalPathSearchTest Begin!\n");
    std::cout << "GlobalPathSearchTest Begin!" << std::endl;
    // OctCube::max_dis = 0.;

    // EXPECT_EQ(0, 0);
    Point3D ptA(2., 2., 2.);
    EXPECT_EQ(p_test->GetIndex(ptA, 3), 3) << std::endl;
    EXPECT_EQ(p_test->GetIndex(ptA, 5), 5) << std::endl;
    std::cout << p_test->GetIndex(ptA, 3) << std::endl;
    std::cout << p_test->GetIndex(ptA, 5) << std::endl;
}

TEST_F(GlobalPathSearchTest, second_gtest) {
    auto FILE = freopen("stdout2", "w", stdout);
    std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
    // OcTree A(10);
    // A.MainUpdate(A.GetRoot(), {}, {{4.26, 1.34, 3.12}});
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    p_test->searcher_ = std::make_shared<Searcher> (A);
    // p_test->searcher_->PrintCubePoolSize();
    printf("GlobalPathSearchTest Begin!\n");
    std::cout << "GlobalPathSearchTest Begin!" << std::endl;

    // A->PrintTree(A->GetRoot());
    std::cout << "first print tree" << std::endl;
    A->MainUpdate(A->GetRoot(), {}, {{4.26, 1.34, 3.12}});
    std::cout << "finish main update" << std::endl;
    // A->PrintTree(A->GetRoot());
    auto all_nodes = A->GetAllLeafNodes();
    std::cout << "all_nodes size: " << all_nodes.size() << std::endl;
    size_t leaf_cnt = 0;
    for(auto node_pair: all_nodes) {
        if(node_pair.second == true) {
            leaf_cnt ++;
        }
    }
    std::cout << "leaf nodes cnt: " << leaf_cnt << std::endl;
    p_test->searcher_->PrintCubePoolSize();
    // CubesBuild(A);
    // std::cout << "check place" << std::endl;
    p_test->searcher_->HalfSolve(Point3D(1., 1., 3.), Point3D(35., 35., 1.8));

    EXPECT_EQ(0, 0);
}

TEST_F(GlobalPathSearchTest, cut_cube_test) {
    auto FILE = freopen("stdout3", "w", stdout);
    std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    p_test->searcher_ = std::make_shared<Searcher> (A);
    A->MainUpdate(A->GetRoot(), {}, {{2.528153, 1.591540, 1.772057}});
    std::vector<size_t> level_leaf_node_cnt(13);
    std::vector<size_t> level_node_cnt(13);
    auto all_nodes = A->GetAllLeafNodes();
    for(auto node_pair: all_nodes) {
        auto oct_node = node_pair.first;
        if (oct_node->level == 1) {
            std::cout << "LEVEL 1 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
            std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
        }
        if (oct_node->level == 2) {
            std::cout << "LEVEL 2 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
            std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
        }
        if(node_pair.second == true) {
            level_leaf_node_cnt[node_pair.first->level] ++;
        } else {
            if (oct_node->level == 1) {
                OctCube test_cube(oct_node, 2, 8, 2);
                std::cout << "test_cube size: x = " << test_cube.Xlen() << " , y = " << test_cube.Ylen() << " , z = " << test_cube.Zlen() << std::endl;
                std::cout << "central point: " << test_cube.CentralPoint().ToString() << std::endl;
                OctCube test_cube2(oct_node, 3, 64, 3);
                std::cout << "test_cube2 size: x = " << test_cube2.Xlen() << " , y = " << test_cube2.Ylen() << " , z = " << test_cube2.Zlen() << std::endl;
                std::cout << "central point2: " << test_cube2.CentralPoint().ToString() << std::endl;
                OctCube test_cube3(oct_node, 4099, 96, 3);
                std::cout << "test_cube3 size: x = " << test_cube3.Xlen() << " , y = " << test_cube3.Ylen() << " , z = " << test_cube3.Zlen() << std::endl;
                std::cout << "central point3: " << test_cube3.CentralPoint().ToString() << std::endl;
            }
        }
        level_node_cnt[node_pair.first->level] ++;
    }
    for (size_t i = 0; i < 12; i++) {
        std::cout << "LEVEL : " << i << " leaf cnt = " << level_leaf_node_cnt[i] << std::endl;
        std::cout << "LEVEL : " << i << " cnt = " << level_node_cnt[i] << std::endl;
    }
}

TEST_F(GlobalPathSearchTest, final_gtest) {
    auto FILE = freopen("stdout4", "w", stdout);
    std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    p_test->searcher_ = std::make_shared<Searcher> (A);
    // A->MainUpdate(A->GetRoot(), {}, {{5.528153, 11.591540, 1.772057}});
    auto all_nodes = A->GetAllLeafNodes();
    for(auto node_pair: all_nodes) {
        auto oct_node = node_pair.first;
        if (oct_node->level == 1) {
            std::cout << "LEVEL 1 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
        }
        if (oct_node->level == 2) {
            std::cout << "LEVEL 2 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
        }
    }
    std::cout << "check index1: " << A->GetIndex(Point3D(38.4, 89.6, 38.4), 1) << std::endl;
    std::cout << "check index2: " << A->GetIndex(Point3U(896, 384, 128), 2) << std::endl;
    std::cout << "check index3: " << A->GetIndex(Point3U(640, 128, 384), 2) << std::endl;
    std::cout << "check index4: " << A->GetIndex(Point3D(38.4, 89.6, 38.4), 2) << std::endl;
    p_test->searcher_->PrintCubePoolSize();
    // CubesBuild(A);
    // std::cout << "check place" << std::endl;
    p_test->searcher_->HalfSolve(Point3D(10., 1., 12.), Point3D(50., 61., 32.));
    p_test->searcher_->AstarTest();
}

TEST_F(GlobalPathSearchTest, addition_gtest1) {
    auto FILE = freopen("stdout5", "w", stdout);
    std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    p_test->searcher_ = std::make_shared<Searcher> (A);
    A->MainUpdate(A->GetRoot(), {}, {{16.528153, 10.591540, 1.772057}});
    A->MainUpdate(A->GetRoot(), {}, {{80.638753, 50.984510, 13.205734}});
    A->CheckValidAll(1);
    A->MainUpdate(A->GetRoot(), {}, {{15.522411, 52.5140, 4.252722}});
    A->MainUpdate(A->GetRoot(), {}, {{68.626423, 63.84510, 6.05664}});
    A->CheckValidAll(2);
    p_test->searcher_->PrintCubePoolSize();
    // CubesBuild(A);
    // std::cout << "check place" << std::endl;
    p_test->searcher_->HalfSolve(Point3D(48., 3., 1.2), Point3D(36., 72., 1.8));
    p_test->searcher_->AstarTest();
    p_test->searcher_->PostProcessing();
}

// TEST_F(GlobalPathSearchTest, addition_gtest2) {
//     auto FILE = freopen("stdout6", "w", stdout);
//     std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
//     std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
//     p_test->searcher_ = std::make_shared<Searcher> (A);
//     std::cout << "start construct obs" << std::endl;
//     // double x;
//     // std::cin >> x;
//     // std::cout << x << std::endl;
//     std::vector<Point3D> pt_vec;
//     for (double x = 30.; x <= 70.; x += 0.5) {
//         for (double z = 0.2; z <= 32.; z += 0.5) {
//             pt_vec.push_back(Point3D(x, 50., z));
//             // std::cout << "x = " << x << " , z = " << z << std::endl;
//         }
//     }
//     // A->MainUpdate(A->GetRoot(), {}, pt_vec);
//     std::vector<size_t> level_leaf_node_cnt(13);
//     std::vector<size_t> level_node_cnt(13);
//     auto all_nodes = A->GetAllLeafNodes();
//     for(auto node_pair: all_nodes) {
//         auto oct_node = node_pair.first;
//         if (oct_node->level == 1) {
//             std::cout << "LEVEL 1 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
//             std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
//         }
//         if (oct_node->level == 2) {
//             std::cout << "LEVEL 2 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
//             std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
//             std::cout << "is leaf = " << node_pair.second << std::endl;
//         }
//         if(node_pair.second == true) {
//             level_leaf_node_cnt[node_pair.first->level] ++;
//         } else {

//         }
//         level_node_cnt[node_pair.first->level] ++;
//     }
//     for (size_t i = 0; i < 12; i++) {
//         std::cout << "LEVEL : " << i << " leaf cnt = " << level_leaf_node_cnt[i] << std::endl;
//         std::cout << "LEVEL : " << i << " cnt = " << level_node_cnt[i] << std::endl;
//     }
//     // A->MainUpdate(A->GetRoot(), {}, {{16.528153, 10.591540, 1.772057}});
//     // A->MainUpdate(A->GetRoot(), {}, {{80.638753, 50.984510, 13.205734}});
//     // A->CheckValidAll(1);
//     // A->MainUpdate(A->GetRoot(), {}, {{15.522411, 52.5140, 4.252722}});
//     // A->MainUpdate(A->GetRoot(), {}, {{68.626423, 63.84510, 6.05664}});
//     // A->CheckValidAll(2);
//     p_test->searcher_->PrintCubePoolSize();
//     // CubesBuild(A);
//     // std::cout << "check place" << std::endl;
//     p_test->searcher_->HalfSolve(Point3D(50., 1., 12.), Point3D(50., 99., 15.));
//     p_test->searcher_->AstarTest();
// }

TEST_F(GlobalPathSearchTest, sce1_gtest) {
    auto FILE = freopen("sce1_stdout", "w", stdout);
    nlohmann::json obs_vizer;
    std::ofstream os;
    os.open("obs.json");
    std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
    std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
    std::cout << "start construct obs" << std::endl;
    std::vector<Point3D> pt_vec;
    auto m_pts0 = p_test->DrawMountain(140., 140., 160., 1.35, 0.65, 1, 62.);
    auto m_pts1 = p_test->DrawMountain(80., 248., 100., 1.6, 1., 1., 62.);
    auto m_pts2 = p_test->DrawMountain(100., 200., 100., 1.35, 1., 1., 50.);
    auto pts0 = p_test->DrawObsCube({166., 216., 12.5}, 15., 10., 25., 1.6);
    auto pts1 = p_test->DrawObsCube({195., 226., 12.5}, 15., 10., 25., 1.35);
    auto pts2 = p_test->DrawObsCube({219., 226., 12.5}, 15., 10., 25., 1.35);
    auto pts3 = p_test->DrawWall({207., 226., 18.}, 10., 3.6, 2.8, 2, 1.35);
    auto pts4 = p_test->DrawObsCube({45., 76., 12.5}, 15., 10., 25., 1.35);
    auto pts5 = p_test->DrawObsCube({69., 76., 12.5}, 15., 10., 25., 1.35);
    auto pts6 = p_test->DrawWall({57., 76., 18.}, 10., 3.6, 2.8, 2, 1.35);
    // auto pts4 = p_test->DrawObsCube({80., 40., 20.01}, 6., 2., 40., 1.);
    pt_vec.insert(pt_vec.end(), m_pts0.begin(), m_pts0.end());
    pt_vec.insert(pt_vec.end(), m_pts1.begin(), m_pts1.end());
    pt_vec.insert(pt_vec.end(), m_pts2.begin(), m_pts2.end());
    pt_vec.insert(pt_vec.end(), pts0.begin(), pts0.end());
    pt_vec.insert(pt_vec.end(), pts1.begin(), pts1.end());
    pt_vec.insert(pt_vec.end(), pts2.begin(), pts2.end());
    pt_vec.insert(pt_vec.end(), pts3.begin(), pts3.end());
    pt_vec.insert(pt_vec.end(), pts4.begin(), pts4.end());
    pt_vec.insert(pt_vec.end(), pts5.begin(), pts5.end());
    pt_vec.insert(pt_vec.end(), pts6.begin(), pts6.end());
    obs_vizer[0]["pt"] = {166., 216., 12.5};
    obs_vizer[0]["len"] = {15., 10., 25.};
    // obs_vizer[1]["pt"] = std::make_tuple(40., 40., 20.);
    // obs_vizer[1]["len"] = std::make_tuple(6., 2., 40.);
    obs_vizer[1]["pt"] = {196., 226., 12.5};
    obs_vizer[1]["len"] = {15., 10., 25.};
    obs_vizer[2]["pt"] = {218., 226., 12.5};
    obs_vizer[2]["len"] = {15., 10., 25.};
    obs_vizer[3]["pt"] = {207., 226., 18.};
    obs_vizer[3]["len"] = {10., 3.6, 2.8};

    obs_vizer[4]["pt"] = {46., 76., 12.5};
    obs_vizer[4]["len"] = {15., 10., 25.};
    obs_vizer[5]["pt"] = {68., 76., 12.5};
    obs_vizer[5]["len"] = {15., 10., 25.};
    obs_vizer[6]["pt"] = {57., 76., 18.};
    obs_vizer[6]["len"] = {10., 3.6, 2.8};

    os << obs_vizer << std::endl;
    // std::cout << "check wa ??" << std::endl;
    /*
    for (double x = 30.; x <= 70.; x += 0.5) {
        for (double z = 0.2; z <= 32.; z += 0.5) {
            pt_vec.push_back(Point3D(x, 50., z));
            // std::cout << "x = " << x << " , z = " << z << std::endl;
        }
    }
    */
    A->MainUpdate(A->GetRoot(), {}, pt_vec);
    std::cout << "main update finish" << std::endl;
    std::vector<size_t> level_leaf_node_cnt(13);
    std::vector<size_t> level_node_cnt(13);
    auto all_nodes = A->GetAllLeafNodes();
    for(auto node_pair: all_nodes) {
        auto oct_node = node_pair.first;
        if (oct_node->level == 1) {
            std::cout << "LEVEL 1 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
            std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
        }
        if (oct_node->level == 2) {
            std::cout << "LEVEL 2 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
            std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
            std::cout << "is leaf = " << node_pair.second << std::endl;
        }
        if(node_pair.second == true) {
            level_leaf_node_cnt[node_pair.first->level] ++;
        } else {

        }
        level_node_cnt[node_pair.first->level] ++;
    }
    for (size_t i = 0; i < 12; i++) {
        std::cout << "LEVEL : " << i << " leaf cnt = " << level_leaf_node_cnt[i] << std::endl;
        std::cout << "LEVEL : " << i << " cnt = " << level_node_cnt[i] << std::endl;
    }
    int flag;
    while(true) {
        int user;
        auto FILE = freopen("sce1_stdout", "w", stdout);
        // std::cout << "input flag" << std::endl;
        // std::cin >> flag;
        flag = 1;
        if (flag == 0) {
            break;
        }
        else if (flag == 1) {
            std::cout << "start running" << std::endl;
            auto start_ = std::chrono::high_resolution_clock::now();
            p_test->searcher_ = std::make_shared<Searcher> (A);
            p_test->searcher_->HalfSolve(Point3D(57., 76., 12.), Point3D(207., 226., 12.));
            p_test->searcher_->AstarTest();
            p_test->searcher_->PostProcessing();
            auto end_ = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_cost = end_ - start_;
            std::cout << "searching time: " << time_cost.count() << std::endl;
            std::cout << "running finish, any key to continue" << std::endl;
            // std::cin >> user;
            break;
        } else {
            std::cout << "input paras and running" << std::endl;
            double heu, ha, hr, ga, gr;
            std::cin >> heu >> ha >> hr >> ga >> gr;
            std::cout << "heu = " << heu << "ha = " << ha << "hr = " << hr << "ga = " << ga << "gr = " << gr << std::endl;
            auto start_ = std::chrono::high_resolution_clock::now();
            p_test->searcher_ = std::make_shared<Searcher> (A);
            p_test->searcher_->HalfSolve(Point3D(5., 5., 30.), Point3D(120., 120., 15.));
            p_test->searcher_->SetParas(heu, ha, hr, ga, gr);
            p_test->searcher_->AstarTest();
            auto end_ = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_cost = end_ - start_;
            std::cout << "searching time: " << time_cost.count() << std::endl;
            std::cout << "running finish, any key to continue" << std::endl;
            std::cin >> user;
        }
    }

    // A->MainUpdate(A->GetRoot(), {}, {{16.528153, 10.591540, 1.772057}});
    // A->MainUpdate(A->GetRoot(), {}, {{80.638753, 50.984510, 13.205734}});
    // A->CheckValidAll(1);
    // A->MainUpdate(A->GetRoot(), {}, {{15.522411, 52.5140, 4.252722}});
    // A->MainUpdate(A->GetRoot(), {}, {{68.626423, 63.84510, 6.05664}});
    // A->CheckValidAll(2);
    // p_test->searcher_->PrintCubePoolSize();
    // CubesBuild(A);
    // std::cout << "check place" << std::endl;
}

// TEST_F(GlobalPathSearchTest, sce2_gtest) {
//     auto FILE = freopen("sce2_stdout", "w", stdout);
//     nlohmann::json obs_vizer;
//     std::ofstream os;
//     os.open("obs.json");
//     std::shared_ptr<GlobalPathSearchTest> p_test = std::make_shared<GlobalPathSearchTest>();
//     std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
//     p_test->searcher_ = std::make_shared<Searcher> (A);
//     std::cout << "start construct obs" << std::endl;
//     std::vector<Point3D> pt_vec;
//     auto pts0 = p_test->DrawWall({66., 20., 20.}, 80., 1., 40., 1, 1.44);
//     auto pts1 = p_test->DrawWall({24., 50., 50.}, 46., 1., 100., 1, 1.25);
//     auto pts2 = p_test->DrawWall({82, 50., 50.}, 48., 1., 100., 1, 1.25);
//     auto pts3 = p_test->DrawWall({50., 80., 30.}, 60., 1., 40., 1, 1.25);
//     // auto pts3 = p_test->DrawWall({})
//     // auto pts0 = p_test->DrawMountain(40., 40., 85., 1.35);
//     // auto pts1 = p_test->DrawObsCube({70., 120., 12.51}, 15., 10., 25., 1.6);
//     // auto pts2 = p_test->DrawObsCube({40., 40., 20.01}, 6., 2., 40., 1.);
//     // auto pts3 = p_test->DrawObsCube({102., 102., 12.51}, 15., 10., 25., 1.6);
//     // auto pts4 = p_test->DrawObsCube({80., 40., 20.01}, 6., 2., 40., 1.);
//     pt_vec.insert(pt_vec.end(), pts0.begin(), pts0.end());
//     pt_vec.insert(pt_vec.end(), pts1.begin(), pts1.end());
//     pt_vec.insert(pt_vec.end(), pts2.begin(), pts2.end());
//     pt_vec.insert(pt_vec.end(), pts3.begin(), pts3.end());
//     // pt_vec.insert(pt_vec.end(), pts4.begin(), pts4.end());
//     obs_vizer[0]["pt"] = std::make_tuple(66., 20., 20.);
//     obs_vizer[0]["len"] = std::make_tuple(80., 1., 40.);
//     obs_vizer[1]["pt"] = std::make_tuple(24., 50., 50.);
//     obs_vizer[1]["len"] = std::make_tuple(47., 1., 100.);
//     obs_vizer[2]["pt"] = std::make_tuple(82, 50., 50.);
//     obs_vizer[2]["len"] = std::make_tuple(49., 1., 100.);
//     obs_vizer[3]["pt"] = std::make_tuple(50., 80., 30.);
//     obs_vizer[3]["len"] = std::make_tuple(60., 1., 40.);

//     // obs_vizer[3]["pt"] = std::make_tuple(80., 40., 20.);
//     // obs_vizer[3]["len"] = std::make_tuple(6., 2., 40.);
//     os << obs_vizer << std::endl;
//     A->MainUpdate(A->GetRoot(), {}, pt_vec);
//     std::cout << "main update finished" << std::endl;
//     auto all_nodes = A->GetAllLeafNodes();
//     std::cout << "all nodes size : " << all_nodes.size() << std::endl;
//     p_test->searcher_->HalfSolve(Point3D(50., 1., 6.), Point3D(52., 128., 12.));
//     p_test->searcher_->AstarTest();
//     p_test->searcher_->PostProcessing();
// }


} // namespace global_path_search
