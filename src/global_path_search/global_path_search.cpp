#include <iostream>
#include <unordered_map>
#include <memory>
#include <nlohmann/json.hpp>

#include "global_path_search.h"
#include "zyb_octree.h"
#include "oct_cube.h"
#include "a_star.h"
#include "types.h"
#include "static_config.h"

namespace global_path_search {

static int xx[] = {0, 1, 0 ,0, 1, 1, 0, 1};
static int yy[] = {0, 0, 1, 0, 1, 0, 1, 1};
static int zz[] = {0, 0, 0, 1, 0, 1, 1, 1};

NodeId Searcher::GetIndex(Point3D pt, int level) {
    Point3U ptu{ull(pt.x / octree_->minimal_size.x), ull(pt.y / octree_->minimal_size.y), ull(pt.z / octree_->minimal_size.z)};
    NodeId base = (1ULL << octree_->MaxLevel) - (1ULL << (octree_->MaxLevel - level));
    return ((((ptu.x & base) << (octree_->MaxLevel << 1)) + ((ptu.y & base) << octree_->MaxLevel) + (ptu.z & base)) << 4) + level;
}

NodeId Searcher::GetCubeIdwithPt(Point3D pt) {
    for (size_t i = octree_->MaxLevel; i - 1 > 0; i--) {
        // std::cout << "i = " << i << std::endl;
        NodeId pro_idx = GetIndex(pt, i);
        if (id2cube_.find(pro_idx) != id2cube_.end()) {
            if (id2cube_[pro_idx] != nullptr)
                // std::cout << pro_idx << std::endl;
                return pro_idx;
        } else {
            // std::cout << "not found" << std::endl;
        }
    }
    return 0;
}

void Searcher::AddOctNode(OcTree::OctNode* oct_node, size_t cur_level, size_t rela_to_node, size_t cur_index) {
    // std::cout << "cur level = " << cur_level << " , cur_index = " << cur_index << std::endl;
    // if (cur_level <= 3) {
        // std::cout << "cur_level = " << cur_level << " ,rela_to_node =  " << rela_to_node << std::endl;
    // }
    if (cur_level >= static_config::maximalCubeLevel) {
        // std::cout << "Add immediately" << std::endl;
        std::unique_ptr<OctCube> oct_cube = std::make_unique<OctCube>(oct_node, cur_level, rela_to_node, cur_index);
        // std::cout << "unique ptr make" << std::endl;
        // std::cout << oct_cube->GetCubeId() << std::endl;;
        if (oct_cube == nullptr) {
            // std::cout << "nullptr warn" << std::endl;
        }
        // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        // std::cout << "pool size : " << oct_cube_pool_.size() << std::endl;
        // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        oct_cube_pool_.emplace_back(std::move(oct_cube));
        // std::cout << "push to pool" << std::endl;
        id2cube_[cur_index] = oct_cube_pool_.back().get();
    }
    else {
        // std::cout << "TOO LARGE CUBE" << std::endl;
        id2cube_[cur_index] = nullptr;
        cur_level ++;
        // bool vis_flag = false;
        // if(oct_node->index == 8589934593) {
            // vis_flag = true;
        // }
        for(size_t i{0U}; i < 8; i++) {
            ull next_index = cur_index;
            if (i & 1) {
                next_index += (1ull << (3 * octree_->MaxLevel - cur_level + 4));
            }
            if (i & 2) {
                next_index += (1ull << ((octree_->MaxLevel << 1) - cur_level + 4));
            }
            if (i & 4) {
                next_index += (1ull << (octree_->MaxLevel - cur_level + 4));
            }
            next_index ++;
            // if (vis_flag) {
                // std::cout << "next_index vis: " << next_index << std::endl;
            // }
            AddOctNode(oct_node, cur_level, rela_to_node * 8 + i, next_index);
        }
    }
    return;
}

void Searcher::HalfSolve(Point3D start_pt, Point3D end_pt) {
    std::cout << "Start HalfSolve" << std::endl;
    auto all_leaf_nodes = octree_->GetAllLeafNodes();
    std::cout << "get all leaf nodes" << std::endl;
    for (auto oct_node_pair : all_leaf_nodes) {
        auto oct_node = oct_node_pair.first;
        bool is_leaf = oct_node_pair.second;
        if (!is_leaf) {
            // std::cout << "not a leaf" << std::endl;
            id2cube_[oct_node->index] = nullptr;
        } else {
            AddOctNode(oct_node, oct_node->level, 1, oct_node->index);
        }
    }
    std::cout << "id2cube size : " << id2cube_.size() << std::endl;
    // Point3D start_pt = Point3D(50., 1., 12.);
    // Point3D end_pt = Point3D(50., 99., 15.);
    NodeId start_cube_idx = GetCubeIdwithPt(start_pt);
    NodeId end_cube_idx = GetCubeIdwithPt(end_pt);
    std::cout << "start cube idx = " << start_cube_idx << std::endl;
    std::cout << "end cube idx = " << end_cube_idx << std::endl;
    if (id2cube_.find(start_cube_idx) != id2cube_.end()) {
        start_cube = id2cube_[start_cube_idx];
        std::cout << "start_cube found" << std::endl;
    }
    if (id2cube_.find(end_cube_idx) != id2cube_.end()) {
        end_cube = id2cube_[end_cube_idx];
        std::cout << "end_cube found" << std::endl;
    }
    a_star_ = std::make_shared<Astar>(start_cube, end_cube, start_pt, end_pt);
}

void Searcher::AstarTest() {
    a_star_->Search(id2cube_);
    res_path_ = a_star_->GetResPath();
    std::cout << "print res path" << std::endl;
    for (auto ver : res_path_) {
        std::cout << "pt: " << ver.GetGeomPt().ToString() << " , ver ori & idx:" << ver.GetOriOnCube() << " ," << ver.GetIdxOnCube() << std::endl;
        std::cout << "idx: " << ver.GetOctCube()->GetCubeId() << std::endl;
    }
    a_star_->JsonWrite();
    return;
}

void Searcher::Solve() {
    auto all_leaf_nodes = octree_->GetAllLeafNodes();
    for (auto oct_node_pair : all_leaf_nodes) {
        auto oct_node = oct_node_pair.first;
        bool is_leaf = oct_node_pair.second;
        if (!is_leaf) {
            id2cube_[oct_node->index] = nullptr;
        }
        AddOctNode(oct_node, oct_node->level, 1, oct_node->index);
        /*
        std::unique_ptr<OctCube> oct_cube = std::make_unique<OctCube>(oct_node);
        if(oct_node->node_type == OctNodeType::SRC) {
            start_cube = oct_cube.get();
        } else if(oct_node->node_type == OctNodeType::DES) {
            end_cube = oct_cube.get();
        }
        oct_cube_pool_.push_back(std::move(oct_cube));
        id2cube_[oct_node->index] = oct_cube_pool_.back().get();
        */
    }
    NodeId start_cube_idx = GetCubeIdwithPt(src_pt);
    NodeId end_cube_idx = GetCubeIdwithPt(des_pt);
    if (id2cube_.find(start_cube_idx) != id2cube_.end()) {
        start_cube = id2cube_[start_cube_idx];
    }
    if (id2cube_.find(end_cube_idx) != id2cube_.end()) {
        end_cube = id2cube_[end_cube_idx];
    }

    a_star_ = std::make_shared<Astar>(start_cube, end_cube, src_pt, des_pt);
    a_star_->Search(id2cube_);
    res_path_ = a_star_->GetResPath();
    return;
}

void Searcher::PostProcessing() {
    std::vector<double> len_vertice;
    len_vertice.push_back(0.);
    for (size_t i{1U}; i < res_path_.size(); i++) {
        len_vertice.push_back(len_vertice.back() + GetDis(res_path_[i-1].GetGeomPt(), res_path_[i].GetGeomPt()));
        std::cout << "len_vertice: " << i << " " << len_vertice[i] << std::endl;
    }
    double total_len = len_vertice.back();
    std::cout << "total_len = " << total_len << std::endl;
    double seg_len = total_len / 20;
    std::vector<size_t> loc_segends(20);
    std::vector<double> len_segends(20);
    std::vector<Point3D> inter_pts(20);
    loc_segends[0] = 0;
    len_segends[0] = 0.;
    size_t path_idx = 1;
    for (size_t i{1U}; i < 20; i++) {
        double cur_len = total_len * i / 20;
        len_segends[i] = cur_len;
        while(len_vertice[path_idx] < cur_len) {
            path_idx ++;
        }
        loc_segends[i] = path_idx - 1;
        std::cout << "loc_segends : " << i << "  " << loc_segends[i] << std::endl;
    }
    std::vector<std::vector<double>> merge_dp;
    merge_dp.resize(20);
    for (size_t i = 0; i < 20; i++) {
        merge_dp[i].resize(4);
    }
    std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> dp_prev;
    inter_pts[0] = res_path_[0].GetGeomPt();
    for (size_t i = 1; i < 20; i++) {
        for (size_t j = 0; j < 4; j++) {
            if (j >= i) {
                continue;
            }
            std::cout << "i = " << i << "  , j = " << j << std::endl;
            double benifit;
            size_t left_end = loc_segends[i - j];
            size_t right_end = loc_segends[i];
            double len_extra_left = len_segends[i - 1] - len_vertice[left_end];
            double len_extra_right = len_segends[i] - len_vertice[right_end];
            if (j == 0) {
                benifit = 0.;
                inter_pts[i] = Interpolation(res_path_[right_end].GetGeomPt(), res_path_[right_end + 1].GetGeomPt(), len_extra_right / (len_vertice[i] - len_vertice[i - 1]));
                std::cout << "inter_pts: " << inter_pts[i].ToString() << std::endl;
            } else {
                // Point3D left_pt = Interpolation(res_path_[i - j -1].GetGeomPt(), res_path_[i - j].GetGeomPt(), len_extra_left / (len_vertice[i - j] - len_vertice[i - j - 1]));
                Point3D left_pt = inter_pts[i - j];
                // Point3D right_pt = Interpolation(res_path_[i - 1].GetGeomPt(), res_path_[i].GetGeomPt(), len_extra_right / (len_vertice[i] - len_vertice[i - 1]));
                Point3D right_pt = inter_pts[i];
                benifit = CalBenifit(len_vertice, left_pt, right_pt, left_end, right_end);
            }
            std::cout << " benifit = " << benifit << std::endl;
            dp_prev[std::make_pair(i, j)] = std::make_pair(i - j - 1, 0);
            for (size_t k = 1; k < 4; k++) {    
                if (merge_dp[i-j-1][k] + benifit > merge_dp[i][j]) {
                    dp_prev[std::make_pair(i, j)] = std::make_pair(i-j-1, k);
                    merge_dp[i][j] = merge_dp[i-j-1][k] + benifit;
                }
                // merge_dp[i][j] = std::max(merge_dp[i][j], merge_dp[i-j-1][k] + benifit);
            }
        }
    }
    std::cout << "dp finish" << std::endl;
    size_t cur_idx = 19;
    size_t cur_st = 0;
    for (size_t k = 1; k < 4; k++) {
        if(merge_dp[cur_idx][cur_st] < merge_dp[cur_idx][k]) {
            cur_st = k;
        }
    }
    path_idx = res_path_.size() - 1;
    while(cur_idx != 0) {
        std::cout << "cur_idx = " << cur_idx << " ,cur_st = " << cur_st << std::endl;
        auto prev_pair = dp_prev[std::make_pair(cur_idx, cur_st)];
        std::cout << "prev_pair idx : " << prev_pair.first << std::endl;
        if (cur_idx == prev_pair.first + 1) {
            cur_idx = prev_pair.first;
            cur_st = prev_pair.second;
            continue;
        } else {
            while(path_idx > loc_segends[cur_idx]) {
                std::cout << "path_idx = " << path_idx << std::endl;
                if(path_idx == 0) {
                    std::cout << "check error!!!" << std::endl;
                    break;
                }
                final_waypts_.push_back(res_path_[path_idx].GetGeomPt());
                path_idx --;
            }
            final_waypts_.push_back(inter_pts[cur_idx]);
            final_waypts_.push_back(inter_pts[prev_pair.first]);
            path_idx = loc_segends[prev_pair.first];
            cur_idx = prev_pair.first;
            cur_st = prev_pair.second;
        }
    }
    if (path_idx > 0) {
        final_waypts_.push_back(res_path_[path_idx--].GetGeomPt());
    }
    std::reverse(final_waypts_.begin(), final_waypts_.end());
    nlohmann::json vizer;
    std::ofstream os;
    os.open("final_waypts.json");
    // for (size_t i{0U}; i < final_waypts_.size(); i++) {
    //     vizer["waypt"][i] = std::make_tuple(final_waypts_[i].x, final_waypts_[i].y, final_waypts_[i].z);
    // }
    for (size_t i = 0; i < final_waypts_.size(); ++i) {
        vizer["waypt"][i] = {final_waypts_[i].x, final_waypts_[i].y, final_waypts_[i].z};
    }
    os << vizer << std::endl;
}

double Searcher::CalBenifit(std::vector<double>& len_vertice, Point3D l_pt, Point3D r_pt, size_t l_idx, size_t r_idx) {
    double origin_len = len_vertice[r_idx] - len_vertice[l_idx] + GetDis(res_path_[r_idx].GetGeomPt(), r_pt);
    double new_len = GetDis(res_path_[l_idx].GetGeomPt(), l_pt) + GetDis(l_pt, r_pt);
    double origin_max_obs_dis = 1e9;
    double new_max_obs_dis = 1e9;
    for (size_t i = l_idx; i <= r_idx; i++) {
        origin_max_obs_dis = std::min(origin_max_obs_dis, res_path_[i].GetOctCube()->GetOctNode()->cost_min);
    }
    // return 0.;
    for (double ra = 0.; ra <= 1.; ra += 0.025) {
        Point3D r_pt = Interpolation(l_pt, r_pt, ra);
        Point3U r_ptu = octree_->GetPoint3U(r_pt);
        auto r_ptu_node = octree_->FindFirstNode(r_ptu, octree_->MaxLevel);
        new_max_obs_dis = std::min(new_max_obs_dis, r_ptu_node->cost_min);
    }
    if (new_max_obs_dis < 4.8 && new_max_obs_dis < origin_max_obs_dis) {
        return -1.;
    }
    else {
        return std::max((origin_len - new_len), -1.);
    }
}


} // namespace global_path_search
