// created by yuqing.wu on 03/06/2024
# pragma once
#include <algorithm>
#include <memory>
#include <vector>
#include <unordered_map>

#include "zyb_octree.h"
#include "types.h"
#include "static_config.h"

#ifndef OCT_CUBE_H
#define OCT_CUBE_H

namespace global_path_search {

typedef ull NodeId;

class OctCubeTest;

class OctCube{
 public:
   const static int max_level = 10;
   static Point3U point3u_bias;
   static Point3D minimal_size;
   static Point3D minimal_half_size;
   // static int bit[8];
   static double max_dis;
   static double max_propagate_dis;

 public:
    OctCube(OcTree::OctNode* o_node, size_t level = 0, size_t rela = 0, size_t psu_index = 0) : oct_node(o_node), rela_to_node(rela), cube_level(level), psu_cube_id(psu_index) {
      // SetCubeLevel(o_node->level);
      if (cube_level <= static_config::secondCubeThreashold) {
         cube_type = CubeType::CUBE_TWO;
      } else {
         cube_type = CubeType::CUBE_ONE;
      }
      Point3D pt_node = o_node->CentralPoint();
      size_t node_level = level;
      std::vector<double> edge_lens = {0.24, 0.24, 0.24};
      size_t cur_rela = rela;
      Point3D cur_pt = pt_node;
      // std::cout << "cur_pt : " << cur_pt.ToString() << std::endl;
      while (cur_rela > 1) {
         size_t dir = cur_rela % 8;
         cur_rela /= 8;
         // std::cout << "dir = " << dir << std::endl;
         // size_t dir_flag[] = { cur_rela | 1, cur_rela | 2, cur_rela | 4 };
         cur_pt = cur_pt + (Point3D{0.12, 0.12, 0.12} * (1<<(max_level - node_level)))[dir];
         Point3D bias = (Point3D{0.12, 0.12, 0.12} * (1<<(max_level - node_level)))[dir];
         // std::cout << bias.ToString() << std::endl;
         node_level --;
         // std::cout << "cur_pt : " << cur_pt.ToString() << std::endl;
      }
      central_pt = cur_pt;
      central_ptu = Point3U(ull(cur_pt.x / 0.24), ull(cur_pt.y / 0.24), ull(cur_pt.z / 0.24));
      x_len_ = edge_lens[0] * (1 << (max_level - level));
      y_len_ = edge_lens[1] * (1 << (max_level - level));
      z_len_ = edge_lens[2] * (1 << (max_level - level));
      cube_id = psu_index;
      // std::cout << "cube finished" << std::endl;
    }
    ~OctCube() = default;

 public:
    size_t GetCubeLevel() const { return cube_level; }
    void CalNeighOctcubes(std::unordered_map<NodeId, OctCube*>& id2cube);
    std::vector<OctCube*> GetAllNeighOctcubes() const { return all_neigh_octcubes_; }
    std::vector<OctCube*> GetNeiOctCubesOnOneSurface(std::unordered_map<NodeId, OctCube*>& id2cube ,size_t ori);
    size_t GetCubeId() const { return cube_id; }
    OcTree::OctNode* GetOctNode() const { return oct_node; }
    CubeType GetCubeType() const { return cube_type; }
    double Xlen() const { return x_len_; }
    double Ylen() const { return y_len_; }
    double Zlen() const { return z_len_; }
    Point3D CentralPoint() const { return central_pt; }
    Point3U CentralPointu() const { return central_ptu; }
    void SetCubeLevel(size_t level) { cube_level = level; }
    size_t GetRelaToNode() const { return rela_to_node; }

 private:
    OcTree::OctNode* oct_node;
    CubeType cube_type;
    size_t cube_level;
    std::vector<OctCube*> all_neigh_octcubes_;
    size_t cube_id;
    double x_len_;
    double y_len_;
    double z_len_;
    size_t rela_to_node;
    NodeId psu_cube_id;
    Point3D central_pt;
    Point3U central_ptu;

   friend class OctCubeTest;
};

} // namespace global_path_search

#endif // OCT_CUBE_H
