// created by yuqing.wu on 04/06/2024

#include "oct_cube.h"

#include <queue>

namespace global_path_search {

static int xx[] = {0, 1, 0, 0, 0, 0, -1};
static int yy[] = {0, 0, 1, 0, 0, -1, 0};
static int zz[] = {0, 0, 0, 1, -1, 0, 0};

static std::vector<std::vector<size_t>> step_metric = {
    { 1, 2 },
    { 0, 2 },
    { 0, 1 },
    { 0, 1 },
    { 0, 2 },
    { 1, 2 }
};

static std::vector<NodeId> GetSonIndex(NodeId idx, size_t level) {
    level ++;
    std::vector<NodeId> res;
    for (size_t i = 0; i < 8; i++) {
        NodeId nxt_idx = idx;
        if (i & 1) {
            nxt_idx += (1ull << (3 * OctCube::max_level - level + 4));
            // nxt_idx += (1ull << (3 * 12 - level + 4));
        }
        if (i & 2) {
            nxt_idx += (1ull << ((OctCube::max_level << 1) - level + 4));
            // nxt_idx += (1ull << ((12 << 1) - level + 4));
        }
        if (i & 4) {
            nxt_idx += (1ull << (OctCube::max_level - level + 4));
            // nxt_idx += (1ull << (12 - level + 4));
        }
        nxt_idx ++;
        if (level <= 4) {
            std::cout << "level = " << level << std::endl;
            std::cout << "nxt_idx = " << nxt_idx << std::endl;
        }
        
        res.push_back(nxt_idx);
    }
    return res;
}

static NodeId GetIndex(Point3U ptu, int level) {
    NodeId base = (1ULL << OctCube::max_level) - (1ULL << (OctCube::max_level - level));
    return ((((ptu.x & base) << (OctCube::max_level << 1)) + ((ptu.y & base) << OctCube::max_level) + (ptu.z & base)) << 4) + level;
    // NodeId base = (1ULL << 12) - (1ULL << (12 - level));
    // return ((((ptu.x & base) << (12 << 1)) + ((ptu.y & base) << 12) + (ptu.z & base)) << 4) + level;

}

void OctCube::CalNeighOctcubes(std::unordered_map<NodeId, OctCube*>& id2cube) {
    // std::vector<OctCube*> res{};
    if (!all_neigh_octcubes_.empty()) {
        return;
    }
    for (size_t i = 0; i < 6; i++) {
        std::cout << "ori = " << i << std::endl;
        
        auto one_surface_res = GetNeiOctCubesOnOneSurface(id2cube, i);
        for(auto sss : one_surface_res) {
            std::cout << "neigh cubes on one face : " << sss->CentralPoint().ToString() << std::endl;
        }
        all_neigh_octcubes_.insert(all_neigh_octcubes_.end(), one_surface_res.begin(), one_surface_res.end());
    }
    /*
    auto nei_oct_nodes = oct_node->GetAllNeighOctNodes();
    for (auto nei_oct_node : nei_oct_nodes) {
        NodeId nei_id = nei_oct_node->index;
        all_neigh_octcubes_.push_back(id2cube[nei_id]);
    }
    */
    return;
}

std::vector<OctCube*> OctCube::GetNeiOctCubesOnOneSurface(std::unordered_map<NodeId, OctCube*>& id2cube, size_t ori) {
    std::vector<OctCube*> res;
    // Point3U next_ptu = central_ptu + (OctCube::point3u_bias * (1 << (OctCube::max_level - cube_level)))[ori];
    // std::cout << "cur cube : " << CentralPoint().ToString() << std::endl;
    // std::cout << " GetNeiOctCubes On one surface , ori = " << ori << std::endl;
    // std::cout << "central_ptu : " << CentralPointu().ToString() << std::endl;
    Point3U next_ptu = central_ptu + (Point3U{1, 1, 1} * (1 << (OctCube::max_level - cube_level)))[ori];
    if (next_ptu.x >= (1 << max_level) * 1) { return {}; }
    if (next_ptu.y >= (1 << max_level) * 1) { return {}; }
    if (next_ptu.z >= (1 << max_level) * 1) { return {}; }
    std::cout << "next_ptu : " << next_ptu.ToString() << std::endl;
    // Point3U next_ptu = central_ptu + (Point3U{1, 1, 1} * (1 << (12 - cube_level)))[ori];
    std::queue<NodeId> q;
    // std::cout << "mmmmmm cube level = " << cube_level << std::endl;
    NodeId nei_id = GetIndex(next_ptu, cube_level);
    std::cout << "nei_id = " << nei_id << std::endl;
    if (id2cube.find(nei_id) == id2cube.end()) {
        // std::cout << "no such id!" << std::endl;
        for (size_t i = cube_level - 1; i >= 0; i--) {
            NodeId t_id = GetIndex(next_ptu, i);
            // std::cout << "t cube level = " << i << std::endl;
            // std::cout << "check t_id = " << t_id << std::endl;
            if (id2cube.find(t_id) != id2cube.end()) {
                // std::cout << "cube level = " << i << std::endl;
                if (id2cube[t_id] != nullptr) {
                    return { id2cube[t_id] };
                } else {
                    return {};
                }
            }
        }
        return {};
    } else {
        std::cout << "get such id" << std::endl;
        q.push(nei_id);
        int cur_level = cube_level;
        while (!q.empty()) {
            std::cout << "queue size: " << q.size() << std::endl;
            size_t q_size = q.size();
            for (size_t i{0U}; i < q_size; i++) {
                NodeId t_id = q.front();
                q.pop();
                if (id2cube.find(t_id) != id2cube.end()) {
                    if (id2cube[t_id] != nullptr) {
                        std::cout << "case A: nei cube valid" << std::endl;
                        res.push_back(id2cube[t_id]);
                        auto c_cube = id2cube[t_id];
                        std::cout << "ccccc center pt : " << c_cube->CentralPoint().ToString() << " ,len: " << c_cube->Xlen() << std::endl;
                    } else {
                        std::cout << "case B: nei cube is smaller" << std::endl;
                        if (cur_level == OctCube::max_level) {
                            std::cout << "the max_level" << std::endl;
                            continue;
                        }
                        std::vector<NodeId> q_son = GetSonIndex(t_id, cur_level);
                        for (int j = 0; j < 8; j++) {
                            if (bit[j] & (1 << ori)) {
                                q.push(q_son[j]);
                            }
                        }
                    }
                }
            }
            cur_level += 1;
        }
    }
    return res;
}

} // namespace global_path_search
