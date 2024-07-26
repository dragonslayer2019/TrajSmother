#include <bits/stdc++.h>
#include "zyb_octree.h"
int tot_flag = 0;

double start, end;
double tmp_start, tmp_end, tmp_tot, tmp_tot1, tmp_tot2;
ull OcTree::GetIndex(Point3U p, int level) {
    ull base = (1ULL << MaxLevel) - (1ULL << (MaxLevel - level));
    return ((((p.x & base) << (MaxLevel << 1)) + ((p.y & base) << MaxLevel) + (p.z & base)) << 4) + level;
}

Point3U OcTree::GetPoint3U(Point3D p) {
    // todo  eps???
    Point3U index{ull(p.x / minimal_size.x), ull(p.y / minimal_size.y), ull(p.z / minimal_size.z)};
    return index;
}

ull OcTree::GetIndex(Point3D p, int level) {
    Point3U uindex = GetPoint3U(p);
    return GetIndex(uindex, level);
}

ull OcTree::GetNextIndex(OctNode* node, Point3U pointu) {
    ull next_index = node->index;
    int level = node->level + 1;
    if (pointu.x > node->pointu.x) {
        next_index += (1ull << (3 * MaxLevel - level + 4));
    }
    if (pointu.y > node->pointu.y) {
        next_index += (1ull << ((MaxLevel << 1) - level + 4));
    }
    if (pointu.z > node->pointu.z) {
        next_index += (1ull << (MaxLevel - level + 4));
    }
    next_index++;
    return next_index;
}

OcTree::OctNode* OcTree::NextNode(OctNode* node, Point3U pointu) {
    // return GetNodePoint(GetNextIndex(node, pointu));
    int idx = 0;
    if (pointu.x > node->pointu.x) {
        idx |= 1;
    }
    if (pointu.y > node->pointu.y) {
        idx |= 2;
    }
    if (pointu.z > node->pointu.z) {
        idx |= 4;
    }
    return node->son[idx];
}

// void SetListRoot(ull index, ListNode* ln) {
//     ListHash[index] = ln;
// }

// ListNode* GetListRoot(ull index) {
//     if (ListHash.find(index) != ListHash.end()) {
//         return ListHash[index];
//     }
//     return nullptr;
// }

void OcTree::DeleteCocList(OctNode* node, ull index, double cost, int idx) {
    // return ;
    if (cost == inf) {
        return ;
    }
    ListNode* ln = node->ln[idx];
    if (ln->pre != nullptr) {
        ln->pre->nxt = ln->nxt;
    }
    else {
        list_hash.Insert(index, ln->nxt);
    }
    if (ln->nxt != nullptr) {
        ln->nxt->pre = ln->pre;
    }
    ln->nxt = ln->pre = nullptr;
}

void OcTree::AddCocList(OctNode* node, ull index, int idx) {
    // return ;
    if (node->cost[idx] == inf) {
        puts("AddCocList Error");
        exit(0);
    }
    ListNode* ln = list_hash.Get(index);
    node->ln[idx]->nxt = ln;
    if (ln != nullptr) {
        ln->pre = node->ln[idx];
    }
    node->ln[idx]->pre = nullptr;
    list_hash.Insert(index, node->ln[idx]);
}

bool OcTree::CheckIn(Point3D p) {
    if (p.x < 0. || p.x > (1 << MaxLevel) * minimal_size.x) return false;
    if (p.y < 0. || p.y > (1 << MaxLevel) * minimal_size.y) return false;
    if (p.z < 0. || p.z > (1 << MaxLevel) * minimal_size.z) return false;
    return true;
}

bool OcTree::CheckIn(Point3U p) {
    if (p.x >= (1 << MaxLevel) * point3u_bias.x) return false;
    if (p.y >= (1 << MaxLevel) * point3u_bias.y) return false;
    if (p.z >= (1 << MaxLevel) * point3u_bias.z) return false;
    return true;
}

void OcTree::UpdateToQueueSingal(OctNode* node, Point3D point, int idx, ull point_idx) {
    if (node == nullptr) {
        puts("nullptr12");
        exit(0);
    }
    Point3D p = node->point + ((minimal_half_size * (1<<(MaxLevel - node->level)))[idx]);
    /*
    if (!CheckIn(p)) {
        std::cout << "p :" << p.ToString() << std::endl;
        puts("!!!!!!!");
        exit(0);
        return ;
    }
    */
    double cost = fd(p, point);
    if (cost + eps < node->cost[idx]) {
        
        // DeleteCocList(node, GetIndex(GetPoint3U(node->coc[idx]), MaxLevel), node->cost[idx], idx);
        DeleteCocList(node, node->coc_idx[idx], node->cost[idx], idx);
        node->cost[idx] = cost;
        node->coc[idx] = point;
        node->coc_idx[idx] = point_idx;
        AddCocList(node, point_idx, idx);
        if (cost > max_propagate_dis * (1ll << (MaxLevel - node->level)) + eps) {
            return ;
        }
        queue[node->level].push(Pair{node, idx, cost, 0});
    }
}

void OcTree::UpdateToQueue(OctNode* node, Point3D point, ull point_idx) {
    if (node == nullptr) {
        puts("nullptr13");
        exit(0);
    }
    double cost_min = node->son_cost_min;
    double cost_max = node->son_cost_max;
    for (int i = 0; i < 8; i++) {
        UpdateToQueueSingal(node, point, i, point_idx);
        cost_min = std::min(cost_min, node->cost[i]);
        if (node->cost[i] != inf) {
            cost_max = std::max(cost_max, node->cost[i]);
        }
    }
    node->cost_min = cost_min;
    node->cost_max = cost_max;
}

void OcTree::UpdateFromNeighborToQueue(OctNode* node, int idx) {
    if (node == nullptr) {
        puts("nullptr13");
        exit(0);
    }
    for (int i = 0; i < 6; i++) {
        auto neighbor_pointu = node->pointu + ((point3u_bias * (1 << (MaxLevel - node->level)))[i]);
        if (!CheckIn(neighbor_pointu)) {
            continue;
        }
        auto neighbor_index = GetIndex(neighbor_pointu, node->level);
        auto neighbor_node = oct_hash.Get(neighbor_index);
        if (neighbor_node == nullptr) {
            continue;
        }
        for (int j = 0; j < 8; j++) {
            if (neighbor_node->cost[j] == inf) {
                continue;
            }
            UpdateToQueueSingal(node, neighbor_node->coc[j], idx, neighbor_node->coc_idx[j]);
        }
    }
    if (node->cost[idx] == inf) {
        queue[node->level].push(Pair{node, idx, node->cost[idx], 0});
    }
}

void OcTree::CheckValidAll(int id) {
    for (auto p : list_hash.UMap) {
        auto node = p.second;
        for (; node!=nullptr; node = node->nxt) {
            // if (GetIndex(node->node->coc[node->idx], MaxLevel) != p.first) {
            if (node->node->coc_idx[node->idx] != p.first) {
                puts("list error");
                exit(0);
            }
            
        }
    }
    for (auto p : oct_hash.UMap) {
        auto node = p.second;
        if (node != root) {
            if (node->level != node->fa->level + 1) {
                printf("check valid%d fail %d %d\n", id, node->level, node->fa->level);
                exit(0);
            }
        }
        for (int i = 0; i < 8; i++) {
            if ((node->son[i] != nullptr) != (node->son[0] != nullptr)) {
                printf("error son %d %lld\n", node->level, node->index);
                for (int j = 0; j < 8; j++) {
                    if (node->son[j] == nullptr) {
                        printf("0");
                    }
                    else {
                        printf("1");
                    }
                    puts("");
                }
                exit(0);
            }
            if (node->son[i] != nullptr) {
                if (node->son[i]->fa != node) {
                    printf("error1 %lld\n", node->index);
                    exit(0);
                }
                
            }
        }
    }
}

OcTree::OctNode* OcTree::SplitToRoot(Point3U pointu, int level) {
    ull index = GetIndex(pointu, level);
    auto node = oct_hash.Get(index);
    if (node != nullptr) {
        return node;
    }
    auto fa = SplitToRoot(pointu, level - 1);
    Split(fa, 1, 1);
    return oct_hash.Get(index);
}

void OcTree::Split(OctNode* node, int flag, int insert_flag) {
    
    if (node->level == MaxLevel) {
        printf("error level!!!\n");
        exit(0);
    }
    if (node->in_merge_queue) {
        printf("error in_merge_queue!!!\n");
        exit(0);
    }
    if (flag) {
        node->in_merge_queue = true;
        merge_queue[node->level + 1].emplace_back(node);
    }
    for (int i = 0; i < 8; i++) {
        if (node->son[i] != nullptr) {
            puts("split error!");
            exit(0);
        }
        tmp_start = clock();
        node->son[i] = new(OctNode);
        tmp_end = clock();
        tmp_tot2 += (tmp_end - tmp_start);
        
        OctNode* x = node->son[i];
        ClearOctNode(x);
        x->fa = node;
        x->level = node->level + 1;
        x->point = node->point + ((minimal_half_size * (1<<(MaxLevel - x->level)))[i]);
        x->pointu = GetPoint3U(x->point);
        x->cost_min = x->son_cost_min = inf;
        x->cost_max = x->son_cost_max = -inf;
        tmp_start = clock();
        x->index = GetIndex(x->pointu, x->level);
        // if (x->index == 5024563242) {
        //     puts("cc");
        // }
        x->in_split_queue = true;
        
        split_queue[x->level].insert(x);
        oct_hash.Insert(x->index, x);
        tmp_end = clock();
        tmp_tot1 += (tmp_end - tmp_start);
        
        tmp_start = clock();
        for (int j = 0; j < 8; j++) {
            // x->ln[j] = new(ListNode);
            // x->ln[j]->node = x;
            // x->ln[j]->idx = j;
            // x->cost[j] = inf;
            // x->coc[j] = Point3D(0., 0., 0.);
            Point3D pj = x->point + ((minimal_half_size * (1<<(MaxLevel - x->level)))[j]);
            for (int k = 0; k < 8; k++) {
                if (node->cost[k] != inf) {
                    double cost = fd(node->coc[k], pj);
                    if (cost < x->cost[j]) {
                        x->cost[j] = cost;
                        x->coc[j] = node->coc[k];
                        x->coc_idx[j] = node->coc_idx[k];
                    }
                }
                
            }
            if (x->cost[j] != inf) {
                // AddCocList(x, GetIndex(GetPoint3U(x->coc[j]), MaxLevel), j);
                AddCocList(x, x->coc_idx[j], j);
                x->cost_min = std::min(x->cost_min, x->cost[j]);
                x->cost_max = std::max(x->cost_max, x->cost[j]);
                // if (insert_flag) {
                //     queue[x->level].push(Pair{x, j, x->cost[j]});
                // }
            }
        }
        tmp_end = clock();
        tmp_tot += (tmp_end - tmp_start);
    }
    
}

void OcTree::UpdateToLeaf(OctNode* node, Point3D point, Point3U pointu, ull point_idx) {
    // printf("level %d\n", node->level);
    if (node->level < MaxLevel) {
        if (node->son[0] == nullptr) {
            Split(node, 1, 1);
            // return ;
        }
    }
    // if (node->index == 5024563242) {
    //     puts("bb");
    // }
    // if (node->level > 0) {
    UpdateToQueue(node, point, point_idx);
    // }
    if (node->level < MaxLevel) {
        for (int i = 0; i < 8; i++) {
            if (node->son[i]->fa != node) {
                printf("AAAA\n");
                exit(0);
            }
            if (node->son[i]->level != node->level + 1) {
                printf("BBBB\n");
                exit(0);
            }
        }
        UpdateToLeaf(NextNode(node, pointu), point, pointu, point_idx);
    }
}

void OcTree::DeletedCoc(std::vector<Point3D> deleted_grid) {
    // struct tmp_pair {
    //     OctNode* node;
    //     int idx;
    // };
    // std::vector<tmp_pair> tt;
    // for (auto grid : deleted_grid) {
    //     for (auto aa : oct_hash.UMap) {
    //         auto node = aa.second;
    //         for (int i = 0; i < 8; i++) {
    //             if (node->coc[i] == grid || GetIndex(GetPoint3U(grid), MaxLevel) == GetIndex(GetPoint3U(node->coc[i]), MaxLevel)) {
    //                 node->coc[i] = {0., 0., 0.};
    //                 node->cost[i] = inf;
    //                 tt.emplace_back(tmp_pair{node, i});
    //             }
    //         }
    //     }
    // }
    // for (auto t : tt) {
    //     UpdateFromNeighborToQueue(t.node, t.idx);
    // }
    // return ;
    for (auto grid : deleted_grid) {
        auto ln = list_hash.Get(GetIndex(GetPoint3U(grid), MaxLevel));
        for (;ln != nullptr; ln = ln->nxt) {
            auto node = ln->node;
            node->coc[ln->idx] = {0., 0., 0.};
            node->coc_idx[ln->idx] = 0;
            node->cost[ln->idx] = inf;
        }
    }
    for (auto grid : deleted_grid) {
        auto ln = list_hash.Get(GetIndex(GetPoint3U(grid), MaxLevel));
        for (;ln != nullptr;) {
            auto nxt = ln->nxt;
            auto node = ln->node;
            ln->nxt = ln->pre = nullptr;
            
            UpdateFromNeighborToQueue(node, ln->idx);
            ln = nxt;
        }
        list_hash.Insert(GetIndex(GetPoint3U(grid), MaxLevel), nullptr);
    }
    // fprintf(stderr, "ccccc\n");
}

inline double Sqr(double x) {
    return x * x;
}

double OcTree::F(double x) {
    double res = 0.6 * pow(std::max(4.8 - x, 0.), 2.) + 60000. * pow(std::max(1.5 - x, 0.), 2.);
    return res;
    // x = std::max(0., std::min(x, max_dis));
    // return Sqr(std::max(0., 1. - x / max_dis));
}

bool OcTree::CanMerge(OctNode* node) {
    if (node == nullptr) {
        puts("nullptr11");
        exit(0);
    }
    // return true;
    double x1 = F(node->cost_max);
    double x2 = F(node->cost_min);
    if ((x2 - x1 < 2.8 or x1 >= 100.)) {
        return true;
    }
    return false;
}

bool OcTree::CanSplit(OctNode* node) {
    return !CanMerge(node);
}

void OcTree::FreeSon(OctNode* node, int idx) {
    if (node == nullptr) {
        puts("nullptr10");
        exit(0);
    }
    auto son = node->son[idx];
    if (son == nullptr) {
        puts("son nullptr");
        exit(0);
    }
    for (int i = 0; i < 8; i++) {
        if (son->level < MaxLevel && son->son[i] != nullptr) {
            FreeSon(son, i);
            // puts("error have son");
            // printf("index %lld %d\n", node->index, node->level);
            // printf("son min max %lf %lf\n", son->cost_min, son->cost_max);
            // printf("min max %lf %lf\n", node->cost_min, node->cost_max);
            // exit(0);
        }
    }
    for (int i = 0; i < 8; i++) {
        
        if (son->cost[i] != inf) {
            // DeleteCocList(son, GetIndex(GetPoint3U(son->coc[i]), MaxLevel), son->cost[i], i);
            DeleteCocList(son, son->coc_idx[i], son->cost[i], i);
        }
        else {
            if (son->coc[i] != Point3D(0., 0., 0.)) {
                puts("delete error2");
                exit(0);    
            }
        }
        if (son->ln[i]->nxt != nullptr || son->ln[i]->pre != nullptr) {
            puts("delete error");
            exit(0);
        }
        delete(son->ln[i]);
        son->ln[i] = nullptr;
    }
    if (son->in_split_queue) {
        split_queue[son->level].erase(son);
    }
    oct_hash.Delete(son->index);
    delete(son);
    node->son[idx] = nullptr;
}

void OcTree::MergeNode(OctNode* node) {
    if (node == nullptr) {
        puts("nullptr9");
        exit(0);
    }
    for (int i = 0; i < 8; i++) {
        if (node->son[i] == nullptr) {
            printf("Error in CanMerge : son[%d] is nullptr!!!\n", i);
            exit(0);
            continue;
        }
        if (!CanMerge(node->son[i])) {
            printf("Error in CanMerge : son[%d] can't merge!!!, node->cost_min / max = %lf / %lf, node->son_cost_min / max = %lf / %lf, son[%d]->cost_min / max = %lf / %lf\n", i, node->cost_min, node->cost_max, node->son_cost_min, node->son_cost_max, i, node->son[i]->cost_min, node->son[i]->cost_max);
            exit(0);
            continue;
        }
        FreeSon(node, i);
    }
    node->son_cost_min = inf;
    node->son_cost_max = -inf;
    UpdateSelf(node);
}

void OcTree::Merge(int level) {
    for (auto node : merge_queue[level]) {
        if (node == nullptr) {
            puts("nullptr8");
            exit(0);
        }
        node->in_merge_queue = false;
        if (CanMerge(node)) {
            MergeNode(node);
        }
    }
}

inline bool SegmentIn(double l, double r, double L, double R) {
    return (l >= L && r <= R);
}

int OcTree::GetBoxStatus(OctNode* node, Point3D lower, Point3D upper) {
    Point3D node_lower = node->point + ((minimal_half_size * (1<<(MaxLevel - node->level)))[0]);
    Point3D node_upper = node->point + ((minimal_half_size * (1<<(MaxLevel - node->level)))[7]);
    if (node_lower.x > upper.x || node_upper.x < lower.x) return 0;
    if (node_lower.y > upper.y || node_upper.y < lower.y) return 0;
    if (node_lower.z > upper.z || node_upper.z < lower.z) return 0;
    if (SegmentIn(node_lower.x, node_upper.x, lower.x, upper.x) && SegmentIn(node_lower.y, node_upper.y, lower.y, upper.y) && SegmentIn(node_lower.z, node_upper.z, lower.z, upper.z)) {
        return 1;
    }
    return 2;
}

inline bool InBox(Point3D point, Point3D lower, Point3D upper) {
    if (point.x < lower.x || point.x > upper.x) return false;
    if (point.y < lower.y || point.y > upper.y) return false;
    if (point.z < lower.z || point.z > upper.z) return false;
    return true;
}

void OcTree::QueryNearestObsToLeaf(OctNode* node, Point3D lower, Point3D upper, Point3D point, Point3U pointu, std::pair<Point3D, double>& best) {
    if (node == nullptr) return ;
    for (int i = 0; i < 8; i++) {
        if (node->cost[i] != inf) {
            double dis = fd(point, node->coc[i]);
            if (!InBox(node->coc[i], lower, upper)) {
                continue;
            }
            if (dis < best.second) {
                best = std::make_pair(node->coc[i], dis);
            }
        }
    }
    if (node->level == MaxLevel) {
        return ;
    }
    QueryNearestObsToLeaf(NextNode(node, pointu), lower, upper, point, pointu, best);
}

void OcTree::QueryNearestObsOnTree(OctNode* node, Point3D lower, Point3D upper, Point3D point, Point3U pointu, std::pair<Point3D, double>& best) {
    if (node == nullptr) return ;
    int status = GetBoxStatus(node, lower, upper);
    if (!status) {
        return ;
    }
    if (status == 1) {
        for (int i = 0; i < 8; i++) {
            if (node->cost[i] != inf) {
                double dis = fd(point, node->coc[i]);
                if (!InBox(node->coc[i], lower, upper)) {
                    continue;
                }
                if (dis < best.second) {
                    best = std::make_pair(node->coc[i], dis);
                }
            }
        }
        if (node->level == MaxLevel) {
            return ;
        }
        QueryNearestObsToLeaf(NextNode(node, pointu), lower, upper, point, pointu, best);
        return ;
    }
    if (node->level == MaxLevel) {
        return ;
    }
    for (int i = 0; i < 8; i++) {
        if (node->son[i] != nullptr) {
            QueryNearestObsOnTree(node->son[i], lower, upper, point, pointu, best);
        }
    }
}


Point3D OcTree::GetNearestObs(double length, double width, double height, Point3D start, Point3D end) {
    Point3D lower{std::min(start.x, end.x) - length, std::min(start.y, end.y) - width, std::min(start.z, end.z) - height};
    Point3D upper{std::max(start.x, end.x) + length, std::max(start.y, end.y) + width, std::max(start.z, end.z) + height};
    std::pair<Point3D, double> best = std::make_pair(Point3D(0., 0., 0.), inf);
    QueryNearestObsOnTree(root, lower, upper, (start + end) * 0.5, GetPoint3U((start + end) * 0.5), best);
    printf("point %lf %lf %lf %lf\n", best.first.x, best.first.y, best.first.z, best.second);
    return best.first;
}

void OcTree::GetObsToLeaf(OctNode* node, Point3D lower, Point3D upper, std::vector<Point3D>& obs) {
    if (node == nullptr) return ;
    if (node->cost_min > cost_eps) return ;
    if (node->level == MaxLevel) {
        obs.emplace_back(node->point);
        return ;
    }
    for (int i = 0; i < 8; i++) {
        if (node->son[i] != nullptr) {
            GetObsToLeaf(node->son[i], lower, upper, obs);
        }
    }
}

void OcTree::GetObsOnTree(OctNode* node, Point3D lower, Point3D upper, std::vector<Point3D>& obs) {
    if (node == nullptr) return ;
    if (node->cost_min > cost_eps) return ;
    int status = GetBoxStatus(node, lower, upper);
    if (!status) {
        return ;
    }
    if (status == 1) {
        if (node->level == MaxLevel) {
            obs.emplace_back(node->point);
            return ;
        }
        GetObsToLeaf(node, lower, upper, obs);
        return ;
    }
    for (int i = 0; i < 8; i++) {
        if (node->son[i] != nullptr) {
            GetObsOnTree(node->son[i], lower, upper, obs);
        }
    }
}

std::vector<Point3D> OcTree::GetObs(double length, double width, double height, Point3D start, Point3D end) {
    Point3D lower{std::min(start.x, end.x) - length, std::min(start.y, end.y) - width, std::min(start.z, end.z) - height};
    Point3D upper{std::max(start.x, end.x) + length, std::max(start.y, end.y) + width, std::max(start.z, end.z) + height};
    std::vector<Point3D> obs{};
    GetObsOnTree(root, lower, upper, obs);
    return obs;
}

void OcTree::Split(int level) {
    // TODO
    for (auto node : split_queue[level]) {
        if (node == nullptr) {
            continue;
        }
        node->in_split_queue = false;
        if (level < MaxLevel) {
            if (CanSplit(node) && node->son[0] == nullptr) {
                Split(node, 0);
            }
        }
    }
}

void OcTree::MergeAndSplit() {
    // fprintf(stderr, "before merge\n");
    for (int i = MaxLevel; i; i--) {
        Merge(i);
        merge_queue[i].clear();
        
    }
    // fprintf(stderr, "after merge\n");
    for (int i = 1; i <= MaxLevel; i++) {
        // printf("split %d\n", i);
        Split(i);
        split_queue[i].clear();
    }
    return ;
}

void OcTree::UpdateSelf(OctNode* node) {
    if (node == nullptr) {
        puts("nullptr7");
        exit(0);
    }
    node->cost_min = node->son_cost_min;
    node->cost_max = node->son_cost_max;
    for (int i = 0; i < 8; i++) {
        node->cost_min = std::min(node->cost_min, node->cost[i]);
        if (node->cost[i] != inf) {
            node->cost_max = std::max(node->cost_max, node->cost[i]);
        }
    }
}

void OcTree::Up(OctNode* node) {
    if (node == nullptr) {
        puts("nullptr6");
        exit(0);
    }
    node->son_cost_min = inf;
    node->son_cost_max = -inf;
    node->cost_min = inf;
    node->cost_max = -inf;

    for (int i = 0; i < 8; i++) {
        auto son = node->son[i];
        if (son == nullptr) {
            printf("error level %d %lld\n", node->level, node->index);
            for (int j = 0; j < 8; j++) {
                if (node->son[j] == nullptr) {
                    printf("0");
                }
                else {
                    printf("1");
                }
                puts("");
            }
            exit(0);
        }
        node->son_cost_min = std::min(node->son_cost_min, son->cost_min);
        node->son_cost_max = std::max(node->son_cost_max, son->cost_max);
        node->cost_min = std::min(node->cost_min, node->cost[i]);
        if (node->cost[i] != inf) {
            node->cost_max = std::max(node->cost_max, node->cost[i]);
        }
    }
    node->cost_min = std::min(node->cost_min, node->son_cost_min);
    node->cost_max = std::max(node->cost_max, node->son_cost_max);
}

void OcTree::UpAll(OctNode* node) {
    node->son_cost_min = inf;
    node->son_cost_max = -inf;
    node->cost_min = inf;
    node->cost_max = -inf;

    for (int i = 0; i < 8; i++) {
        auto son = node->son[i];
        node->son_cost_min = std::min(node->son_cost_min, son->cost_min);
        node->son_cost_max = std::max(node->son_cost_max, son->cost_max);
        for (int j = 0; j < 8; j++) {
            Point3D pj = node->point + ((minimal_half_size * (1<<(MaxLevel - node->level)))[j]);
            for (int k = 0; k < 8; k++) {
                if (son->cost[k] != inf) {
                    double cost = fd(son->coc[k], pj);
                    if (cost + eps < node->cost[j]) {
                        // DeleteCocList(node, GetIndex(GetPoint3U(node->coc[j]), MaxLevel), node->cost[j], j);
                        DeleteCocList(node, node->coc_idx[j], node->cost[j], j);
                        node->cost[j] = cost;
                        node->coc[j] = son->coc[k];
                        node->coc_idx[j] = son->coc_idx[k];
                        AddCocList(node, node->coc_idx[j], j);
                        queue[node->level].push(Pair{node, j, cost, 0});
                    }
                }
            }
            node->cost_min = std::min(node->cost_min, node->cost[j]);
            if (node->cost[j] != inf) {
                node->cost_max = std::max(node->cost_max, node->cost[j]);
            }
        }
    }
    node->cost_min = std::min(node->cost_min, node->son_cost_min);
    node->cost_max = std::max(node->cost_max, node->son_cost_max);
}

void OcTree::UpSingal(OctNode* node, OctNode* son, int idx) {
    if (node == nullptr) {
        puts("nullptr5");
        exit(0);
    }

    if (son == nullptr) {
        puts("son nullptr1");
        exit(0);
    }
    
    for (int j = 0; j < 8; j++) {
        Point3D pj = node->point + ((minimal_half_size * (1 << (MaxLevel - node->level)))[j]);
        double cost = fd(son->coc[idx], pj);
        if (cost + eps < node->cost[j]) {
            // DeleteCocList(node, GetIndex(GetPoint3U(node->coc[j]), MaxLevel), node->cost[j], j);
            DeleteCocList(node, node->coc_idx[j], node->cost[j], j);
            node->cost[j] = cost;
            node->coc[j] = son->coc[idx];
            node->coc_idx[j] = son->coc_idx[idx];
            // AddCocList(node, GetIndex(GetPoint3U(node->coc[j]), MaxLevel), j);
            AddCocList(node, node->coc_idx[j], j);
            queue[node->level].push(Pair{node, j, cost, 0});
        }
    }
}

void OcTree::MainUpdate(OctNode* root, std::vector<Point3D> deleted_grid, std::vector<Point3D> added_grid) {
    //CheckValidAll(0);
    // fprintf(stderr, "aaa\n");
    DeletedCoc(deleted_grid);
    fprintf(stderr, "AAA\n");
    //CheckValidAll(0);
    for (auto grid : added_grid) {
        // printf("added %lf %lf %lf %llu %llu %llu\n", grid.x, grid.y, grid.z, GetPoint3U(grid).x, GetPoint3U(grid).y, GetPoint3U(grid).z);
        UpdateToLeaf(root, grid, GetPoint3U(grid), GetIndex(GetPoint3U(grid), MaxLevel));
    }
    fprintf(stderr, "BBB\n");
    
    //CheckValidAll(1);
    for (int i = MaxLevel; i >= 0; i--) {
        fprintf(stderr, "dji %d\n", i);
        Dij(i);
        // end = clock();
        fprintf(stderr, "time %lf\n", (end - start) / CLOCKS_PER_SEC);
    }
    // fprintf(stderr, "CCC\n");
    //CheckValidAll(3);
    MergeAndSplit();
    // fprintf(stderr, "DDD\n");
}

void OcTree::PrintBinary(int x, int level) {
    if (level > 1) {
        PrintBinary(x, level - 1);
    }
    printf("%d", (x & (1 << (MaxLevel - level))) > 0);
}

void OcTree::Dij(int level) {
    int tot = 0;
    // std::vector<Pair> tmp{};
    update_queue[level].clear();
    fprintf(stderr, "before dij %d %d\n", level, queue[level].size());
    Pair Last;
    int tot1 = 0;
    while (!queue[level].empty()) {
        tot++;
        auto pa = queue[level].top();
        
        queue[level].pop();
        if (pa == Last) {
            tot1++;
            continue;
        }
        Last = pa;
        auto node = pa.node;
        if (node == nullptr) {
            puts("nullptr1");
            exit(0);
        }
        if (node->level != (node->index & 15)) {
            printf("error in index %d %lld\n", node->level, (node->index & 15));
            exit(0);
        }
        // if (pa.node->index == 5024563242 && pa.idx == 3) {
        //     printf("aa %lf\n", pa.cost);
        // }
        if (pa.cost > node->cost[pa.idx] + eps) {
            tot1++;
            continue;
        }
        // tmp.emplace_back(pa);
        if (!pa.node->in_update_queue) {
            pa.node->in_update_queue = true;
            update_queue[level].emplace_back(pa.node);
        }
        if (pa.cost == inf) continue;
        if (pa.cost > max_propagate_dis * (1ll << (MaxLevel - level)) + eps) continue;
        if (pa.cost > max_dis + eps) continue;
        
        auto pointu = node->pointu;
        // for (int i = 0; i < 6; i++) {
        //     auto bias = ((point3u_bias * (1 << (MaxLevel - level)))[i]);
        //     auto next_pointu = pointu + bias;
        //     if (!CheckIn(next_pointu)) {
        //         continue;
        //     }
        //     auto next_index = GetIndex(next_pointu, level);
        //     auto next_node = oct_hash.Get(next_index);
        //     if (next_node == nullptr) {
        //         // continue;
        //         next_node = SplitToRoot(next_pointu, level);
        //         // if (next_node->index == 5024563242) {
        //         //     printf("cc\n");
        //         // }
        //     }
        // }
        // if (!queue[level].empty()) {
        //     if (queue[level].top().cost + eps < pa.cost) {
        //         fprintf(stderr, "%lf %lf %llu %llu\n", queue[level].top().cost, pa.cost, queue[level].top().node->index, pa.node->index);
        //         queue[level].push(pa);
        //         fprintf(stderr, "after %lf %lf %llu %llu\n", queue[level].top().cost, pa.cost, queue[level].top().node->index, pa.node->index);
        //         continue;
        //     }
            
        // }
        for (int i = 0; i < 6; i++) {
            auto bias = ((point3u_bias * (1 << (MaxLevel - level)))[i]);
            auto next_pointu = pointu + bias;
            int nxt_idx = pa.idx;
            if (next_pointu.x != pointu.x) {
                 if (!((next_pointu.x > pointu.x) == ((pa.idx & 1) > 0))) {
                    continue;
                 }
                 nxt_idx ^= 1;
            }
            if (next_pointu.y != pointu.y) {
                 if (!((next_pointu.y > pointu.y) == ((pa.idx & 2) > 0))) {
                    continue;
                 }
                 nxt_idx ^= 2;
            }
            if (next_pointu.z != pointu.z) {
                 if (!((next_pointu.z > pointu.z) == ((pa.idx & 4) > 0))) {
                    continue;
                 }
                 nxt_idx ^= 4;
            }
            
            if (!CheckIn(next_pointu)) {
                continue;
            }
            auto next_index = GetIndex(next_pointu, level);
            auto next_node = oct_hash.Get(next_index);
            if (next_node == nullptr) {
                // continue;
                // puts("error1");
                // exit(0);
                next_node = SplitToRoot(next_pointu, level);
            }
            if (next_node->level != level) {
                printf("level error !!!!! %d %d\n", next_node->level, level);
                exit(0);
            }
            // if (pa.cost > next_node->cost[nxt_idx] + eps) {
            //     puts("error");
            //     printf("%.30lf %.30lf\n", pa.cost, next_node->cost[nxt_idx]);
            //     printf("%llu %llu %llu %llu\n", next_node->index, nxt_idx, pa.node->index, pa.idx);
            //     exit(0);
            // }
            
            UpdateToQueueSingal(next_node, node->coc[pa.idx], nxt_idx, node->coc_idx[pa.idx]);
        }
        for (int i = 0; i < 3; i++) {
            int nxt_idx = (pa.idx ^ (1 << i));
            UpdateToQueueSingal(node, node->coc[pa.idx], nxt_idx, node->coc_idx[pa.idx]);
        }
    }
    fprintf(stderr, "in dij %d %d %d %d\n", level, tot, tot1, queue[level].MaxSize());
    // end = clock();
    // fprintf(stderr, "dij1\n");
    //CheckValidAll(100 + level);
    for (auto node : update_queue[level]) {
        UpdateSelf(node);
        node->in_update_queue = false;
    }
    end = clock();
    // fprintf(stderr, "dij2\n");
    if (level > 0) {
        // fprintf(stderr, "in1\n");
        for (auto node : update_queue[level]) {
            std::set<Point3D> se;
            for (int i = 0; i < 8; i++) {
                if (node->cost[i] != inf) {
                    se.insert(node->coc[i]);
                }
            }
            for (auto coc : se) {
                UpdateToQueue(node->fa, coc, GetIndex(GetPoint3U(coc), MaxLevel));
            }
            if (!node->fa->in_merge_queue) {
                node->fa->in_merge_queue = true;
                merge_queue[level].emplace_back(node->fa);
            }
            if (!node->in_split_queue) {
                node->in_split_queue = true;
                split_queue[level].insert(node);
            }
        }
        // for (auto tmp_pa : tmp) {
        //     if (tmp_pa.node == nullptr) {
        //         puts("nullptr2");
        //         exit(0);
        //     }
        //     // UpdateSelf(tmp_pa.node);
        //     if (tmp_pa.node->cost[tmp_pa.idx] != inf) {
        //         UpSingal(tmp_pa.node->fa, tmp_pa.node, tmp_pa.idx);
        //     }
        
        //     if (tmp_pa.node->fa == nullptr) {
        //         puts("nullptr3");
        //         exit(0);
        //     }   
        //     if (!tmp_pa.node->fa->in_merge_queue) {
        //         tmp_pa.node->fa->in_merge_queue = true;
        //         merge_queue[level].emplace_back(tmp_pa.node->fa);
        //         // if (tmp_pa.node->fa->index == 5) {
        //         //     puts("insert b");
        //         //     printf("node level %d %d %lld %lld\n", tmp_pa.node->fa->index, tmp_pa.node->level, tmp_pa.node->fa, oct_hash.Get(5));
        //         //     if (tmp_pa.node->fa->son[0] != nullptr) {
        //         //         puts("have son");
        //         //     }
        //         //     else {
        //         //         puts("no son");
        //         //     }
        //         // }
        //         if (tmp_pa.node->fa->level == MaxLevel) {
        //             printf("error !!! %d %d\n", tmp_pa.node->level, tmp_pa.node->fa->level);
        //             exit(0);
        //         }
        //     }
        //     if (!tmp_pa.node->in_split_queue) {
        //         tmp_pa.node->in_split_queue = true;
        //         split_queue[level].insert(tmp_pa.node);
        //     }
        // }
        //CheckValidAll(200 + level);
        // fprintf(stderr, "in3\n");
        for (auto node : merge_queue[level]) {
            // node->in_merge_queue = false;
            if (node == nullptr) {
                puts("nullptr4");
                exit(0);
            }
            Up(node);
        }
        //CheckValidAll(300 + level);
    }
    else {
        // fprintf(stderr, "in2\n");
        // for (auto tmp_pa : tmp) {
        //     if (tmp_pa.node == nullptr) {
        //         puts("nullptr40");
        //         exit(0);
        //     }
        //     if (!tmp_pa.node->in_split_queue) {
        //         tmp_pa.node->in_split_queue = true;
        //         split_queue[level].insert(tmp_pa.node);
        //     }
        // }
        for (auto node : update_queue[level]) {
            if (!node->in_split_queue) {
                node->in_split_queue = true;
                split_queue[level].insert(node);
            }
        }
    }
    // fprintf(stderr, "dij3\n");
    //CheckValidAll(300 + level);
}

void OcTree::PrintIndex(ull index) {
    int x = 0, y = 0, z = 0, level;
    level = int(index & 15);
    index >>= 4;
    for (int i = 1; i <= MaxLevel; i++) {
        x <<= 1;
        y <<= 1;
        z <<= 1;
        if (i <= level) {
            x += int((index & (1ull << (3 * MaxLevel - i))) > 0ull);
            y += int((index & (1ull << (2 * MaxLevel - i))) > 0ull);
            z += int((index & (1ull << (MaxLevel - i))) > 0ull);
        }
    }
    // printf("xyz %d %d %d\n", x, y, z);
    printf("x : "); 
    PrintBinary(x, level);
    printf(" , ");
    printf("y : ");
    PrintBinary(y, level);
    printf(" , ");
    printf("z : ");
    PrintBinary(z, level);
    printf(" , ");
    printf("level : %d, index %llu\n", level, index);
}

void OcTree::PrintTree(OctNode* x) {
    PrintIndex(x->index);
    auto p1 = x->point + ((minimal_half_size * (1<<(MaxLevel - x->level)))[0]);
    auto p2 = x->point + ((minimal_half_size * (1<<(MaxLevel - x->level)))[7]);
    printf("x : [ %lf , %lf ], y : [ %lf , %lf ], z : [ %lf , %lf ]\n", p1.x, p2.x, p1.y, p2.y, p1.z, p2.z);
    printf("cost min : %lf , cost max : %lf\n", x->cost_min, x->cost_max);
    for (int i = 0; i < 8; i++) {
        if (x->son[i] != nullptr) {
            PrintTree(x->son[i]);
        }
    }
}

std::vector<std::pair<OcTree::OctNode*, bool>> OcTree::GetAllLeafNodes() {
    std::vector<std::pair<OctNode*, bool>> res;
    PushBackLeafNodes(root ,res);
    return res;
}

void OcTree::PushBackLeafNodes(OctNode* node, std::vector<std::pair<OctNode*, bool>>& res) {
    bool isLeaf = true;
    for (int i = 0; i < 8; i++) {
        if (node->son[i] != nullptr) {
            isLeaf = false;
            PushBackLeafNodes(node->son[i], res);
        }
    }
    // if (isLeaf) {
    res.push_back(std::make_pair(node, isLeaf));
    // }
}

void OcTree::SplitToLeaf(OctNode* node) {
    if (node->level == MaxLevel) {
        return ;
    }

}

double OcTree::QueryCost(OctNode* node, Point3D point, Point3U pointu) {
    double cost = inf;
    for (int i = 0; i < 8; i++) {
        if (node->cost[i] != inf) {
            cost = std::min(cost, fd(point, node->coc[i]));
        }
    }
    if (node->level == MaxLevel) {
        return cost;
    }
    auto next_index = GetIndex(pointu, node->level + 1);
    auto next_node = oct_hash.Get(next_index);
    if (next_node != nullptr) {
        cost = std::min(cost, QueryCost(next_node, point, pointu));
    }
    return cost;
    
}

void OcTree::Query(std::vector<Point3D> query_queue) {
    for (auto point : query_queue) {
        printf("point : %lf %lf %lf , dis : %lf\n", point.x, point.y, point.z, QueryCost(root, point, GetPoint3U(point)));
    }
}

OcTree::OctNode* OcTree::FindFirstNode(Point3U pointu, int level) {
    for (int i = level - 1; i >= 0; i--) {
        auto node = oct_hash.Get(GetIndex(pointu, i));
        if (node != nullptr) {
            return node;
        }
    }
    printf("Error! not find!!! point : (%llu, %llu, %llu), level : %d\n", pointu.x, pointu.y, pointu.z, level);
    return nullptr;
}

std::vector<OcTree::OctNode*> OcTree::GetNeiOctNodesOnOneSurface(OctNode* node, size_t ori) {
    auto next_pointu = node->pointu + (point3u_bias * (1 << (MaxLevel - node->level)))[ori];
    if (!CheckIn(next_pointu)) {
        return {};
    }
    std::vector<OctNode*> nodes{};
    std::queue<OctNode*> tmp;
    auto t = oct_hash.Get(GetIndex(next_pointu, node->level));
    if (t == nullptr) {
        return {FindFirstNode(next_pointu, node->level - 1)};
    }
    tmp.push(t);
    while (!tmp.empty()) {
        auto p = tmp.front();
        tmp.pop();
        if (p->level == MaxLevel || p->son[0] == nullptr) {
            nodes.emplace_back(p);
            continue;
        }
        for (int i = 0; i < 8; i++) {
            if (bit[i] & (1 << ori)) {
                tmp.push(p->son[i]);
            }
        }
    }
    return nodes;
}

std::vector<OcTree::OctNode*> OcTree::GetAllNeighOctNodes(OctNode* node) {
    std::vector<OctNode*> nodes{};
    for (int i = 0; i < 6; i++) {
        auto tmp = GetNeiOctNodesOnOneSurface(node, i);
        nodes.insert(nodes.end(), tmp.begin(), tmp.end());
    }
    return nodes;
}

// int main() {
//     // freopen("c1.in", "r", stdin);
//     start=clock();
//     OcTree A(10, {0.1, 0.1, 0.1}, {0.05, 0.05, 0.05}, {1, 1, 1}, 1e20, 0.6, 4.8, 0.0866);
//     int tot;
//     std::vector<Point3D> InsertQueue, DeleteQueue;
//     // scanf("%d", &tot);
//     // for (int i = 0; i < tot; i++) {
//     //     tot_flag = i;
//     //     char c[10];
//     //     double x, y, z;
//     //     scanf("%s%lf%lf%lf", c, &x, &y, &z);
//     //     if (c[0] == 'I') {
//     //         InsertQueue.emplace_back(x, y, z);
//     //         A.MainUpdate(A.GetRoot(), {}, {{x, y, z}});
//     //     }
//     //     else if (c[0] == 'D') {
//     //         // continue;
//     //         DeleteQueue.emplace_back(x, y, z);
//     //         A.MainUpdate(A.GetRoot(), {{x, y, z}}, {});
//     //     }
//     //     else if (c[0] == 'Q') {
//     //         // A.MainUpdate(A.GetRoot(), DeleteQueue, InsertQueue);
//     //         InsertQueue.clear();
//     //         DeleteQueue.clear();
//     //         A.Query({{x, y, z}});
//     //     }
//     //     // fprintf(stderr, "finish %d\n", i);
//     //     // printf("finish %d\n", i);
//     //     fprintf(stderr, "hash size %d %d\n", (int)(A.oct_hash.Size()), i);
//     //     A.CheckValidAll(i);
//     // }
//     // MainUpdate(root, {}, {{0, 0, 0}});
//     // A.MainUpdate(A.GetRoot(), {}, {{4.26, 1.34, 3.12}});
//     std::vector<Point3D> pt_vec;
//     for (double x = 30.; x <= 70.; x += 0.5) {
//         for (double z = 0.2; z <= 32.; z += 0.5) {
//             pt_vec.push_back(Point3D(x, 50., z));
//             // std::cout << "x = " << x << " , z = " << z << std::endl;
//         }
//     }
//     A.MainUpdate(A.GetRoot(), {}, pt_vec);
//     end = clock();
//     puts("------");
//     printf("%f\n", (end - start) / CLOCKS_PER_SEC);
//     printf("split time %f %f %f\n", tmp_tot / CLOCKS_PER_SEC, tmp_tot1 / CLOCKS_PER_SEC, tmp_tot2 / CLOCKS_PER_SEC);
//     // MainUpdate(root, {{0, 0, 0}}, {});
//     A.Query({{30., 49., 0.2}, {30., 49., 0.1}});
//     // Query({{10, 0.2, 0.2}});

//     // PrintTree(root);
// }