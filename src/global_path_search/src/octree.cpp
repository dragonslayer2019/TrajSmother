// created by yuqing.wu on 06/06/2024
//
// a temp file for running global_path_search algorithm before the real octree module been prepared
//

/*
#include "octree.h"

#include <vector>
#include <memory>

namespace global_path_search {


std::vector<OctNode*> OctNode::GetAllNeighOctNodes() {
    return {};
}

std::vector<std::vector<OctNode*>> OctNode::GetNeiOctNodesOnOneSurface(size_t ori, int partition_cnt) {
    return {};
}


OctNode* GetNodePointFromHashMap(ull index) {
    if (hash.find(index) != hash.end()) {
        // printf("find %llu\n", index);
        return hash[index];
    }
    return nullptr;
}


ull GetIndex(Point3U p, int level) {
    ull base = (1ULL << maxlevel) - (1ULL << (maxlevel - level));
    return ((((p.x & base) << (maxlevel << 1)) + ((p.y & base) << maxlevel) + (p.z & base)) << 4) + level;
}

ull GetIndex(Point3D p, int level) {
    Point3U uindex = GetPoint3U(p);
    return GetIndex(uindex, level);
}

bool CheckIn(Point3D p) {
    if (p.x < 0. || p.x > (1 << maxlevel) * minimal_size.x) return false;
    if (p.y < 0. || p.y > (1 << maxlevel) * minimal_size.y) return false;
    if (p.z < 0. || p.z > (1 << maxlevel) * minimal_size.z) return false;
    return true;
}

bool CheckIn(Point3U p) {
    if (p.x > (1 << maxlevel) * point3u_bias.x) return false;
    if (p.y > (1 << maxlevel) * point3u_bias.y) return false;
    if (p.z > (1 << maxlevel) * point3u_bias.z) return false;
    return true;
}


} // namespace global_path_search
*/