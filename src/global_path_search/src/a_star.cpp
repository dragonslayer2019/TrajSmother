//created by yuqing.wu on 04/06/2024

#include <nlohmann/json.hpp>

#include "a_star.h"

namespace global_path_search {

SearchResult Astar::Search(std::unordered_map<NodeId, OctCube*>& id2cube) {
    std::cout << "start search" << std::endl;
    size_t sp_cnt = 0;
    size_t spp_cnt = 0;
    std::shared_ptr<StateNode> cur_node;
    while(!open_set_.empty()) {
        fuses ++;
        if (fuses >= 800000) {
            break;
        }
        std::cout << "searching!!!" << std::endl;
        cur_node = open_set_.top();
        open_set_.pop();
        std::cout << "OPEN SET view: geom_pt = " << cur_node->ver.GetGeomPt().ToString() <<  
        " ,cube = " << cur_node->ver.GetOctCube()->CentralPoint().ToString() << std::endl;
        std::cout << "cur cube len = " << cur_node->ver.GetOctCube()->Xlen() << " , cube_level = " << cur_node->ver.GetOctCube()->GetCubeLevel() << std::endl;
        cur_node->node_state = NodeState::IN_CLOSET_SET;
        // std::cout << "wha??" << std::endl;
        if(ReachDestination(cur_node)) {
            std::cout << "reach destination !" << std::endl;
            return SearchResult::SUCCESS;
        }
        if(NearDestination(cur_node)) {
            std::cout << "Near Destination!!" << std::endl;
            Vertice des_ver(end_cube, 10, 10 * (1 << static_config::cellsCntEachEdgeOnCube2) * (1 << static_config::cellsCntEachEdgeOnCube2));
            des_ver.SetGeomPt(end_pt);
            double des_g_score = cur_node->g_score + CalTransitionCost(cur_node->ver, des_ver);
            double des_f_score = des_g_score;
            auto des_state_node = std::make_shared<StateNode>(des_ver, des_f_score, des_g_score, cur_node);
            std::shared_ptr<StateNode> pre_node = state_hash_.find(des_ver);
            if (pre_node == nullptr) {
                des_state_node->node_state = NodeState::IN_OPEN_SET;
                state_hash_.insert(des_state_node);
                open_set_.push(des_state_node);
            } else {
                if (pre_node->node_state == NodeState::IN_CLOSET_SET) {
                    continue;
                } else {
                    if (pre_node->f_score > des_f_score) {
                        state_hash_.insert(des_state_node);
                    }
                }
            }
            continue;
        }
        
        auto cur_ver = cur_node->ver;
        // std::cout << "cur ver" << std::endl;
        std::vector<Vertice> cur_neighs = cur_ver.GetAllEdges(id2cube);
        std::cout << "cur_neighs size : " << cur_neighs.size() << std::endl;
        bool sp_flag = false;
        bool spp_flag = false;
        if (fuses == 1) {
            sp_flag = true;
            auto ccube = cur_ver.GetOctCube();
            std::cout << "cube central : " << ccube->CentralPoint().ToString() << std::endl;
            std::cout << "cube level : " << ccube->GetCubeLevel() << " , rela = " << ccube->GetRelaToNode() << std::endl;
            std::cout << "cur ver (ori & idx) = " << cur_ver.GetOriOnCube() << " ," << cur_ver.GetIdxOnCube() << std::endl;
            std::cout << "node centre : " << ccube->GetOctNode()->CentralPoint().ToString() << std::endl;
            std::cout << "cubeId : " << ccube->GetCubeId() << std::endl;
        }
        if (fuses == 2) {
            spp_flag = true;
            auto ccube = cur_ver.GetOctCube();
            std::cout << "cube central : " << ccube->CentralPoint().ToString() << std::endl;
            std::cout << "cube level : " << ccube->GetCubeLevel() << " , rela = " << ccube->GetRelaToNode() << std::endl;
            std::cout << "cur ver (ori & idx) = " << cur_ver.GetOriOnCube() << " ," << cur_ver.GetIdxOnCube() << std::endl;
            std::cout << "node centre : " << ccube->GetOctNode()->CentralPoint().ToString() << std::endl;
            std::cout << "cubeId : " << ccube->GetCubeId() << std::endl;
        }
        for (auto neigh_ver : cur_neighs) {
            double neigh_g_score = cur_node->g_score + CalTransitionCost(cur_ver, neigh_ver);
            double neigh_f_score = neigh_g_score + Heuristic(neigh_ver);
            // std::cout << "ver pt: " << neigh_ver.GetGeomPt().ToString() << std::endl;
            // std::cout << "f(x) = " << neigh_f_score << std::endl;
            // std::cout << "g(x) = " << neigh_g_score << std::endl;
            // std::cout << "h(x) = " << neigh_f_score - neigh_g_score << std::endl;

            if (sp_flag) {
                // vizer["sp"][sp_cnt++] = std::make_tuple(neigh_ver.GetGeomPt().x, neigh_ver.GetGeomPt().y, neigh_ver.GetGeomPt().z);
                std::cout << "SPECIALLLL" << std::endl;
                std::cout << "GEOM PT = " << neigh_ver.GetGeomPt().ToString() << std::endl;
                std::cout << "g_score = " << neigh_g_score << std::endl;
                std::cout << "f_score = " << neigh_f_score << std::endl;
                std::cout << "ori & idx = " << neigh_ver.GetOriOnCube() << " ," << neigh_ver.GetIdxOnCube() << std::endl;
                std::cout << "cubeId : " << neigh_ver.GetOctCube()->GetCubeId() << std::endl;
            }
            if (spp_flag) {
                // vizer["spp"][spp_cnt++] = std::make_tuple(neigh_ver.GetGeomPt().x, neigh_ver.GetGeomPt().y, neigh_ver.GetGeomPt().z);
                std::cout << "SPPPPPPPPP" << std::endl;
                std::cout << "GEOM PT = " << neigh_ver.GetGeomPt().ToString() << std::endl;
                std::cout << "g_score = " << neigh_g_score << std::endl;
                std::cout << "f_score = " << neigh_f_score << std::endl;
                std::cout << "ori & idx = " << neigh_ver.GetOriOnCube() << " ," << neigh_ver.GetIdxOnCube() << std::endl;
                std::cout << "cubeId : " << neigh_ver.GetOctCube()->GetCubeId() << std::endl;
            }
            auto neigh_state_node = std::make_shared<StateNode>(neigh_ver, neigh_f_score, neigh_g_score, cur_node);
            std::shared_ptr<StateNode> pre_node = state_hash_.find(neigh_ver); 
            if (pre_node == nullptr) {
                // std::cout << "NUUL ptr" << std::endl;
                neigh_state_node->node_state = NodeState::IN_OPEN_SET;
                state_hash_.insert(neigh_state_node);
                open_set_.push(neigh_state_node);
            } else {
                if (pre_node->node_state == NodeState::IN_CLOSET_SET) {
                    // std::cout << "closed! " << std::endl;
                    continue;
                } else {
                    if (pre_node->f_score > neigh_f_score) {
                        // std::cout << "Insert! " << std::endl;
                        state_hash_.insert(neigh_state_node);
                    }
                    std::cout << "give up" << std::endl;
                }
            }
        }
    }
    std::cout << "FAILURE !!" << std::endl;
    return SearchResult::FAILURE;
}

double Astar::Heuristic(Vertice& ver) {
    double origin_h_score = GetDis(ver.GetGeomPt(), end_pt) + 0.7 * abs(ver.GetGeomPt().z - end_pt.z) + 0.3 * std::max(0., 10. - ver.GetGeomPt().z);
    if (ver.GetOctCube()->GetCubeLevel() >= 5) {
        origin_h_score = h_small_grid_cost_ratio * (origin_h_score + h_small_grid_cost_add);
    }
    return heu_ratio * origin_h_score;
}

double Astar::CalTransitionCost(Vertice& origin_ver, Vertice& new_ver) {
    double obs_cost = CalObsCost(new_ver);
    double tran_score = GetDis(origin_ver.GetGeomPt(), new_ver.GetGeomPt())  * (1. + obs_cost);
    if (obs_cost >= 1e8) {
        tran_score = 1e9;
    }
    tran_score += abs(origin_ver.GetGeomPt().z - new_ver.GetGeomPt().z) * 0.7;
    if (new_ver.GetOctCube()->GetCubeLevel() >= 5) {
        // tran_score *= 1.18;
        tran_score = (tran_score + g_small_grid_cost_add) * g_small_grid_cost_ratio;
    }
    // std::cout << "transition score = " << tran_score << std::endl;
    return tran_score;
}

double Astar::CalObsCost(Vertice& ver) {
    // return 0.;
    double dis_max = ver.GetOctCube()->GetOctNode()->cost_min;
    if (dis_max <= 0.6) {
        return 1e9;
    }
    double res = 0.6 * pow(std::max(4.8 - dis_max, 0.), 2.) + 60000. * pow(std::max(1.5 - dis_max, 0.), 2.);
    std::cout << "cal obs cost = " << res << std::endl;
    return res;
}

bool Astar::ReachDestination(std::shared_ptr<StateNode> node) {
    if(node->ver.GetOriOnCube() == 10) {
        CalResPath(node);
        return true;
    }
    return false;
}

bool Astar::NearDestination(std::shared_ptr<StateNode> node) {
    if (node->ver.GetOctCube() == end_cube) {
        return true;
    }
    return false;
}

void Astar::CalResPath(std::shared_ptr<StateNode> node) {
    // res_path_.push_back(node->ver);
    res_path_.clear();
    std::cout << "CalResPath" << std::endl;
    std::shared_ptr<StateNode> cur_ptr = node;
    while(cur_ptr->parent_node != nullptr) {
        std::cout << "cur state node: f_score = " << cur_ptr->f_score << " ,g_score = " << cur_ptr->g_score << std::endl;
        std::cout << "cube pt = " << cur_ptr->ver.GetOctCube()->CentralPoint().ToString() << " ,cube level = " << cur_ptr->ver.GetOctCube()->GetCubeLevel() << std::endl;
        CalObsCost(cur_ptr->ver);
        res_path_.push_back(cur_ptr->ver);
        cur_ptr = cur_ptr->parent_node;
    }
    std::cout << "cur state node: f_score = " << cur_ptr->f_score << " ,g_score = " << cur_ptr->g_score << std::endl;
    res_path_.push_back(cur_ptr->ver);
    std::reverse(res_path_.begin(), res_path_.end());
}

void Astar::JsonWrite() const {
    nlohmann::json vizer;
    std::ofstream os;
    os.open("astar_result.json");
    vizer.clear();
    for (size_t i{0U}; i < res_path_.size(); i++) {
        vizer["traj"][i] = {res_path_[i].GetGeomPt().x, res_path_[i].GetGeomPt().y, res_path_[i].GetGeomPt().z};
    }
    vizer["traj"][res_path_.size()] = {end_pt.x, end_pt.y, end_pt.z};
    // vizer["obs"][0] = {50., 50., 16.};
    vizer["end_point"][0] = {start_cube->CentralPoint().x, start_cube->CentralPoint().y, start_cube->CentralPoint().z};
    vizer["end_point"][1] = {end_pt.x, end_pt.y, end_pt.z};
    vizer["ha"] = h_small_grid_cost_add;
    vizer["hr"] = h_small_grid_cost_ratio;
    vizer["ga"] = g_small_grid_cost_add;
    vizer["gr"] = g_small_grid_cost_ratio;
    vizer["her"] = heu_ratio;
    os << vizer << std::endl;
}


void Astar::ShortCut() {
    return;
}

} // namespace global_path_search
