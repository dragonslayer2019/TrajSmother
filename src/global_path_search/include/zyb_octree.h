# pragma once
#include <bits/stdc++.h>

#include "geom.h"
#include "types.h"

#ifndef ZYB_OCTREE_H
#define ZYB_OCTREE_H

template<class Node>
class HashMap {
    public:
        Node* GetNodePointFromHashMap(ull index);
        void Insert(ull index, Node* node) {
            // if (index == 729595946) {
            //     fprintf(stderr, "dd\n");
            // }
            if (node == nullptr) {
                Delete(index);
                return ;
            }
            UMap[index] = node;
        }
        Node* Get(ull index) {
            if (UMap.find(index) != UMap.end()) {
                return UMap[index];
            }
            return nullptr;
        }
        void Delete(ull index) {
            // printf("delete! %d\n", (int)(hash.size()));
            UMap.erase(index);
        }
        size_t Size() {
            return UMap.size();
        }
    std::unordered_map<ull, Node*> UMap;
    // private:
    //     std::unordered_map<ull, Node*> UMap;
};




// const int MaxLevel = 12;
// const Point3D minimal_size = {0.1, 0.1, 0.1};
// const Point3D minimal_half_size = {0.05, 0.05, 0.05};
// const Point3U point3u_bias = {1, 1, 1};
// const double inf = 1e20;
// const double max_propagate_dis = 0.6;
// const double max_dis = 3.;
// const double cost_eps = fd(minimal_half_size, {0., 0., 0.});
const int bit[8] = {7 /* 000111 */, 38 /* 100110 */, 21 /* 010101 */, 52 /* 110100 */, 11 /* 001011 */, 42 /* 101010 */, 25 /* 011001 */, 56 /* 111000 */};




class OcTree {
    
    public:
        const int MaxLevel = 12;
        const Point3D minimal_size = {0.24, 0.24, 0.24};
        const Point3D minimal_half_size = {0.12, 0.12, 0.12};
        const Point3U point3u_bias = {1, 1, 1};
        const double inf = 1e20;
        double max_propagate_dis = 0.4;
        double max_dis = 3.;
        const double cost_eps = 0.207846; // 0.0866 = 0.05 * sqrt(3);
        static constexpr double eps = 0;

    public :
        OcTree(int MaxLevel_input = 12, Point3D minimal_size_input = {0.24, 0.24, 0.24}, Point3D minimal_half_size_input = {0.12, 0.12, 0.12},
               Point3U point3u_bias_input = {1, 1, 1},  double inf_input = 1e20, double max_propagate_dis_input = 0.4,
               double max_dis_input = 3., double cost_eps_input = 0.207846):
               MaxLevel(MaxLevel_input), minimal_size(minimal_size_input), minimal_half_size(minimal_half_size_input),
               point3u_bias(point3u_bias_input), inf(inf_input), max_propagate_dis(max_propagate_dis_input), max_dis(max_dis_input),
               cost_eps(cost_eps_input) {
            root = new(OctNode);
            ClearOctNode(root);
            root->point = minimal_half_size * (1 << (MaxLevel));
            root->pointu = GetPoint3U(root->point);
            root->level = 0;
            root->index = GetIndex(root->pointu, root->level);
            oct_hash.Insert(root->index, root);
            queue.resize(MaxLevel + 1);
            merge_queue.resize(MaxLevel + 1);
            split_queue.resize(MaxLevel + 1);
            update_queue.resize(MaxLevel + 1);
            oct_hash.UMap.reserve(10000000);
            list_hash.UMap.reserve(10000000);
            merge_queue[MaxLevel].reserve(10000000);
            split_queue[MaxLevel].reserve(10000000);
            update_queue[MaxLevel].reserve(10000000);
        }
        struct OctNode;
        struct ListNode {
            ListNode* pre{nullptr};
            ListNode* nxt{nullptr};
            OctNode* node;
            int idx;
        };
        struct Pair {
            OctNode* node;
            int idx;
            double cost;
            int id;
            // friend bool operator < (const Pair &x, const Pair &y) {
            //     if  (x.cost != y.cost) {
            //         return x.cost < y.cost;
            //     }
            //     if (x.node->index != y.node->index) {
            //         return x.node->index < y.node->index;
            //     }
            //     return x.idx < y.idx;
            // }
            friend bool operator > (const Pair &x, const Pair &y) {
                if (fabs(x.cost - y.cost) >= OcTree::eps) {
                    return x.cost > y.cost;
                }
                if (x.node->index != y.node->index) {
                    return x.node->index > y.node->index;
                }
                return x.idx > y.idx;
            }

            friend bool operator == (const Pair &x, const Pair &y) {
                if (fabs(x.cost - y.cost) >= OcTree::eps) return false;
                if (x.node->index != y.node->index) return false;
                if (x.idx != y.idx) return false;
                return true;
            }
        };
        struct OctNode {
            ListNode* ln[8];
            OctNode* son[8];
            OctNode* fa;
            double cost[8], cost_min, cost_max, son_cost_min, son_cost_max;
            Point3D coc[8];
            ull coc_idx[8];
            ull index;
            int level;
            Point3D point;
            Point3U pointu;
            bool in_merge_queue = false;
            bool in_split_queue = false;
            bool in_update_queue = false;
            Pair* pa[8];
            Point3D CentralPoint() const {
                return point;
            }
        };

        bool IsOctNodeValid(OctNode* node) const {
            return node->cost_min >= cost_eps; // TODO
        }
        void ClearOctNode(OctNode* node) {
            for (int i = 0; i < 8; i++) {
                node->ln[i] = new(ListNode);
                node->ln[i]->node = node;
                node->ln[i]->idx = i;
                node->ln[i]->nxt = node->ln[i]->pre = nullptr;
                node->son[i] = nullptr;
                node->coc[i] = Point3D(0., 0., 0.);
                node->coc_idx[i] = 0;
                node->cost[i] = inf;
                node->pa[i] = nullptr;
            }
            node->fa = nullptr;
            node->cost_min = node->son_cost_min = inf;
            node->cost_max = node->son_cost_max = -inf;
        }

        OctNode* GetRoot() {
            return root;
        }
        ull GetIndex(Point3U p, int level);
        Point3U GetPoint3U(Point3D p);
        ull GetIndex(Point3D p, int level);
        ull GetNextIndex(OctNode* node, Point3U pointu);
        bool CheckIn(Point3D p);
        bool CheckIn(Point3U p);
        void MainUpdate(OctNode* root, std::vector<Point3D> deleted_grid, std::vector<Point3D> added_grid);
        void PrintBinary(int x, int level);
        void PrintIndex(ull index);
        void PrintTree(OctNode* x);
        void Query(std::vector<Point3D> query_queue);
        OctNode* FindFirstNode(Point3U pointu, int level);
        std::vector<OctNode*> GetNeiOctNodesOnOneSurface(OctNode* node, size_t ori);
        std::vector<OctNode*> GetAllNeighOctNodes(OctNode* node);

        HashMap<OctNode> oct_hash;
        HashMap<ListNode> list_hash;

        void CheckValidAll(int id);
        std::vector<std::pair<OctNode*, bool>> GetAllLeafNodes();
        Point3D GetNearestObs(double length, double width, double height, Point3D start, Point3D end);
        std::vector<Point3D> GetObs(double length, double width, double height, Point3D start, Point3D end);

    private :
        void PushBackLeafNodes(OctNode* node, std::vector<std::pair<OctNode*, bool>>& res);
        int GetBoxStatus(OctNode* node, Point3D lower, Point3D upper);
        void QueryNearestObsToLeaf(OctNode* node, Point3D lower, Point3D upper, Point3D point, Point3U pointu, std::pair<Point3D, double>& best);
        void QueryNearestObsOnTree(OctNode* node, Point3D lower, Point3D upper, Point3D point, Point3U pointu, std::pair<Point3D, double>& best);
        void GetObsToLeaf(OctNode* node, Point3D lower, Point3D upper, std::vector<Point3D>& obs);
        void GetObsOnTree(OctNode* node, Point3D lower, Point3D upper, std::vector<Point3D>& obs);

        OctNode* root;
        
        struct Heap {

            // void up(int x) {
            //     x++;
            //     while (x > 1) {
            //         if (*qu[x / 2 - 1] > *qu[x - 1]) {
            //             std::swap(qu[x - 1], qu[x / 2 - 1]);
            //             qu[x - 1]->id = x - 1;
            //             qu[x / 2 - 1]->id = x / 2 - 1;
            //             x = x / 2;
            //         }
            //         else break;
            //     }
            // }
            // void down(int x) {
            //     x++;
            //     while (x * 2 <= qu.size()) {
            //         int j = x * 2;
            //         if (j + 1 < qu.size() && *qu[j - 1] > *qu[j]) {
            //             j++;
            //         }
            //         if (*qu[x - 1] > *qu[j - 1]) {
            //             std::swap(qu[x - 1], qu[j - 1]);
            //             qu[x - 1]->id = x - 1;
            //             qu[j - 1]->id = j - 1;
            //             x = j;
            //         }
            //         else break;
            //     }
            // }

            // void pop() {
            //     std::swap(qu[0], qu.back());
            //     qu[0]->id = 0;
            //     qu.back()->node->pa[qu.back()->idx] = nullptr;
            //     delete(qu.back());
            //     qu.pop_back();
            //     if (qu.size() > 0) {
            //         down(0);
            //     }
            // }
            // bool empty() {
            //     return qu.size() == 0;
            // }

            // Pair top() {
            //     return *qu[0];
            // }

            // void push(Pair x) {
            //     if (x.node->pa[x.idx] == nullptr) {
            //         Pair* pa = new(Pair);
            //         *pa = x;
            //         x.node->pa[x.idx] = pa;
            //         qu.emplace_back(pa);
            //         pa->id = qu.size() - 1;
            //         // up(qu.size() - 1);
            //     }
            //     else {
            //         x.node->pa[x.idx]->cost = x.cost;
            //         // up(x.node->pa[x.idx]->id);
            //     }
            //     max_size = std::max(max_size, (int)qu.size());
            // }
            // void clear() {
            //     qu.clear();
            // }
            // int size() {
            //     return qu.size();
            // }
            void push(Pair x) {
                if (x.node->pa[x.idx] == nullptr) {
                    Pair* pa = new(Pair);
                    *pa = x;
                    x.node->pa[x.idx] = pa;
                    de.push_back(pa);
                    // pa->id = de.size() - 1;
                    // up(qu.size() - 1);
                }
                else {
                    x.node->pa[x.idx]->cost = x.cost;
                    // up(x.node->pa[x.idx]->id);
                }
                max_size = std::max(max_size, (int)de.size());
            }
            void clear() {
                de.clear();
            }
            int size() {
                return de.size();
            }
            int MaxSize() {
                return max_size;
            }
            bool empty() {
                return de.size() == 0;
            }

            Pair top() {
                return *de[0];
            }
            void pop() {
                de.front()->node->pa[de.front()->idx] = nullptr;
                delete(de.front());
                de.pop_front();
            }
            std::vector<Pair*> qu;
            std::deque<Pair*> de;
            int max_size = 0;
        };
        std::vector<Heap> queue;
        std::vector<std::vector<OctNode*>> merge_queue;
        std::vector<std::unordered_set<OctNode*>> split_queue;
        std::vector<std::vector<OctNode*>> update_queue;

        void Dij(int level);
        void Split(OctNode* node, int flag = 1, int insert_flag = 0);
        OctNode* NextNode(OctNode* node, Point3U pointu);
        void DeleteCocList(OctNode* node, ull index, double cost, int idx);
        void AddCocList(OctNode* node, ull index, int idx);
        void UpdateToQueueSingal(OctNode* node, Point3D point, int idx, ull point_idx);
        void UpdateToQueue(OctNode* node, Point3D point, ull point_idx);
        void UpdateFromNeighborToQueue(OctNode* node, int idx);
        OctNode* SplitToRoot(Point3U pointu, int level);
        void UpdateToLeaf(OctNode* node, Point3D point, Point3U pointu, ull point_idx);
        void DeletedCoc(std::vector<Point3D> deleted_grid);
        double F(double x);
        bool CanMerge(OctNode* node);
        bool CanSplit(OctNode* node);
        void FreeSon(OctNode* node, int idx);
        void MergeNode(OctNode* node);
        void Merge(int level);
        void Split(int level);
        void MergeAndSplit();
        void Up(OctNode* node);
        void UpAll(OctNode* node);
        void UpSingal(OctNode* node, OctNode* son, int idx);
        void SplitToLeaf(OctNode* node);
        double QueryCost(OctNode* node, Point3D point, Point3U pointu);
        void UpdateSelf(OctNode* node);
        

};

#endif // ZYB_OCTREE_H
