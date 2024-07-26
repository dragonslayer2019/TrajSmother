// created by yuqing.wu on 05/06/2024
# pragma once
#include <memory>
#include <queue>
#include <vector>

#include "zyb_octree.h"
#include "oct_cube.h"
#include "a_star.h"
#include "esdf.h"
#include "vertice.h"

#ifndef GLOBAL_PATH_SEARCH_H
#define GLOBAL_PATH_SEARCH_H

namespace global_path_search {

class GlobalPathSearchTest;

class Searcher {
 public:
    Searcher(std::shared_ptr<OcTree> octree) : octree_(octree) {
      std::cout << "searcher construct" << std::endl;
      std::cout << "oct cube pool size: " << oct_cube_pool_.size() << std::endl;
      // OctCube::max_level = 0;
      /*
      for (size_t i{0U}; i < 8; i++) {
         OctCube::bit[i] = bit[i];
      }
      OctCube::max_level = octree_->MaxLevel;
      OctCube::minimal_half_size = octree_->minimal_half_size;
      OctCube::point3u_bias = octree_->point3u_bias;
      OctCube::minimal_size = octree_->minimal_size;
      OctCube::max_dis = octree_->max_dis;
      OctCube::max_propagate_dis = octree_->max_propagate_dis;
      */
    }
    void Solve();
    void HalfSolve(Point3D start_pt, Point3D end_pt); // only used for test
    void AstarTest(); // only used for test
    void PrintCubePoolSize() { std::cout << "oct cube pool size: " << oct_cube_pool_.size() << std::endl; }
    void SetParas(double heu, double ha, double hr, double ga, double gr) {
      a_star_->SetParas(heu, ha, hr, ga, gr);
    }
    void PostProcessing();

 public:
    void AddOctNode(OcTree::OctNode* oct_node, size_t cur_level = 0, size_t rela_to_node = 0, size_t cur_index = 0);
    NodeId GetCubeIdwithPt(Point3D pt);
    NodeId GetIndex(Point3D pt, int level);
    std::vector<Point3D> GetFinalWaypts() const  { return final_waypts_; }

 private:
    double CalBenifit(std::vector<double>& len_vertice, Point3D l_pt, Point3D r_pt, size_t l_idx, size_t r_idx);

 private:
    std::shared_ptr<OcTree> octree_;
    std::shared_ptr<Astar> a_star_{nullptr};
    std::vector<std::unique_ptr<OctCube>> oct_cube_pool_;
    std::vector<Vertice> res_path_;
    OctCube* end_cube{nullptr};
    OctCube* start_cube{nullptr};
    std::unordered_map<NodeId, OctCube*> id2cube_;
    Point3D src_pt;
    Point3D des_pt;
    std::vector<Point3D> final_waypts_;

   friend class GlobalPathSearchTest;
};

} // namespace global_path_search

#endif // GLOBAL_PATH_SEARCH_H
