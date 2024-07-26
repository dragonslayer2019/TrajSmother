//created by yuqing.wu on 04/06/2024

#include <algorithm>
#include <memory>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>

#include "vertice.h"

#ifndef A_STAR_H
#define A_STAR_H

namespace global_path_search {

class AstarTest;

class StateNode {
 public:
    StateNode(Vertice vert, double f_s, double g_s, std::shared_ptr<StateNode> p_node) : ver(vert), f_score(f_s), g_score(g_s), parent_node(p_node) {} 
 public:
    Vertice ver;
    double f_score;
    double g_score;
    std::shared_ptr<StateNode> parent_node;
    NodeState node_state{NodeState::UNVISITED};
};

struct StateNodeCmp {
    bool operator()(std::shared_ptr<StateNode> node_a, std::shared_ptr<StateNode> node_b) {
        return node_a->f_score > node_b->f_score;
    }
};

class StateNodeHash {
 public:
    std::shared_ptr<StateNode> find(Vertice ver) {
      if (state_map_.find(std::make_pair(ver.GetOctCube()->GetCubeId(), ver.GetIdxOnCube())) == state_map_.end()) {
         return nullptr;
      }
      return state_map_[std::make_pair(ver.GetOctCube()->GetCubeId(), ver.GetIdxOnCube())];
        // return nullptr;
    }
    void insert(std::shared_ptr<StateNode> node) {
      state_map_[std::make_pair(node->ver.GetOctCube()->GetCubeId(), node->ver.GetIdxOnCube())] = node;
        return;
    }
 private:
    std::map<std::pair<NodeId, size_t>, std::shared_ptr<StateNode>> state_map_;
};

class Astar {
 public:
    Astar(OctCube* start_cube, OctCube* end_cube, Point3D src, Point3D des) : start_cube(start_cube), end_cube(end_cube) {
      std::cout << "build astar" << std::endl;
      Vertice start_vertice(start_cube, 0);
      start_vertice.SetGeomPt(start_cube->CentralPoint());
      std::shared_ptr<StateNode> start_node_state = std::make_shared<StateNode>(start_vertice, Heuristic(start_vertice), 0., nullptr);
      std::cout << "astar get start" << std::endl;
      start_node_state->node_state = NodeState::IN_OPEN_SET;
      open_set_.push(start_node_state);
      fuses = 0;
      std::cout << "end cube = " << end_cube->CentralPoint().ToString() << std::endl;
      end_pt = des;
    }
    ~Astar() = default;
 public:
    SearchResult Search(std::unordered_map<NodeId, OctCube*>& id2cube);
    double Heuristic(Vertice& ver);
    double CalTransitionCost(Vertice& origin_ver, Vertice& new_ver);
    bool ReachDestination(std::shared_ptr<StateNode> node);
    bool NearDestination(std::shared_ptr<StateNode> node);
    std::vector<Vertice> GetResPath() const { return res_path_; }
    double CalObsCost(Vertice& ver);
    void JsonWrite() const;
    void SetParas(double heu, double ha, double hr, double ga, double gr) {
      heu_ratio = heu;
      h_small_grid_cost_add = ha;
      h_small_grid_cost_ratio = hr;
      g_small_grid_cost_add = ga;
      g_small_grid_cost_ratio = gr;
    }


 private:
    StateNodeHash state_hash_;
    std::priority_queue<std::shared_ptr<StateNode>, std::vector<std::shared_ptr<StateNode> >, StateNodeCmp> open_set_;
    std::vector<Vertice> res_path_;
    OctCube* start_cube;
    OctCube* end_cube;
    Point3D src_pt;
    Point3D end_pt;
    size_t fuses;
    double heu_ratio{1.06};
    double h_small_grid_cost_add{0.2};
    double h_small_grid_cost_ratio{1.08};
    double g_small_grid_cost_add{0.};
    double g_small_grid_cost_ratio{1.02};


 private:
    void CalResPath(std::shared_ptr<StateNode> final_node);
    void ShortCut(); // API remained for res path short cut
};

} // global_path_search

#endif // A_STAR_H