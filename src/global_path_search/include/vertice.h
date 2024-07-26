// created by yuqing.wu on 04/06/2024

#include <vector>
#include <queue>
#include <memory>

#include "types.h"
#include "oct_cube.h"
#include "static_config.h"
#include "geom.h"

#ifndef VERTICE_H
#define VERTICE_H

namespace global_path_search {

class Vertice {
 public:
    Vertice(OctCube* cube, size_t ori_idx, size_t idx = 0) : oct_cube(cube), ori_on_cube(ori_idx), idx_on_cube(idx) {
      if (ori_on_cube == 0) {
         vertice_type = VerticeType::CENTER;
      } else {
         vertice_type = VerticeType::SURFACE;
      }
    }
    ~Vertice() = default;

    std::vector<Vertice> GetAllEdges(std::unordered_map<NodeId, OctCube*>& id2cube);
    size_t GetIdxOnCube() const { return idx_on_cube; }
    std::vector<std::vector<Vertice>> GetVerticesOnOneSurface(size_t ori, int partition_cnt);
    OctCube* GetOctCube() const { return oct_cube; }
    size_t GetOriOnCube() const { return ori_on_cube; }
    Point3D GetGeomPt() const { return geom_pt; }
    void SetGeomPt(Point3D pt) { geom_pt = pt; }
    void SetOri(size_t ori) { ori_on_cube = ori; }
    void SetIdx(size_t idx) { idx_on_cube = idx; }

 private:
    OctCube* oct_cube;
    VerticeType vertice_type;
    size_t idx_on_cube;
    size_t ori_on_cube; // 0 for center, 1~6 for surfaces
    Point3D geom_pt;

 private:
    Point3D CalRegionDiff(size_t ori, size_t idx);
};

} // namespace global_path_search

#endif // VERTICE_H
