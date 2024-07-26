// created by yuqing.wu on 04/06/2024

#include "vertice.h"

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

namespace global_path_search {

Point3D Vertice::CalRegionDiff(size_t ori, size_t idx) {
    size_t seg_cnt = (1 << static_config::cellsCntEachEdgeOnCube2);
    // std::cout << "ori = " << ori << " , idx = " << idx << std::endl;
    size_t dir1 = idx / seg_cnt;
    size_t dir2 = idx % seg_cnt;
    
    // std::cout << "dir1 = " << dir1 << std::endl;
    // std::cout << "dir2 = " << dir2 << std::endl;

    std::vector<double> edge_lens = {oct_cube->Xlen(), oct_cube->Ylen(), oct_cube->Zlen()};
    std::vector<double> ori_diff = {0., 0., 0.};
    ori_diff[step_metric[ori - 1][0]] = -edge_lens[step_metric[ori - 1][0]] / 2. + edge_lens[step_metric[ori - 1][0]] / seg_cnt * (dir1 + 0.5) ;
    ori_diff[step_metric[ori - 1][1]] = -edge_lens[step_metric[ori - 1][1]] / 2. + edge_lens[step_metric[ori - 1][1]] / seg_cnt * (dir2 + 0.5) ;
    return Point3D(ori_diff[0], ori_diff[1], ori_diff[2]);
}

std::vector<Vertice> Vertice::GetAllEdges(std::unordered_map<NodeId, OctCube*>& id2cube) {
    std::vector<Vertice> res;
    // std::cout << "Get All Edges" << std::endl;
    size_t cur_cube_level = oct_cube->GetCubeLevel();
    if (ori_on_cube == 0) {
        // std::cout << "  ver 1" << std::endl;
        oct_cube->CalNeighOctcubes(id2cube);
        auto neigh_cubes = oct_cube->GetAllNeighOctcubes();
        // std::cout << "neigh_cubes size : " << neigh_cubes.size() << std::endl;
        for (auto neigh_cube : neigh_cubes) {
            // std::cout << "neight_cube cpt: " << neigh_cube->CentralPoint().ToString() << std::endl;
            Vertice nver(neigh_cube, 0, 0);
            nver.geom_pt = neigh_cube->CentralPoint();
            res.push_back(nver);
        }
        if (oct_cube->GetCubeType() == CubeType::CUBE_TWO) {
            // std::cout << "origin node: " << oct_cube->GetOctNode()->CentralPoint().ToString() << std::endl;
            // std::cout << "I am CUBE TWO" << std::endl;
            // std::cout << "rela to node = " << oct_cube->GetRelaToNode() << std::endl;
            std::vector<double> cube_centre = {oct_cube->CentralPoint().x, oct_cube->CentralPoint().y, oct_cube->CentralPoint().z};
            // std::cout << "cube center: " << cube_centre[0] << " ," << cube_centre[1] << " ," << cube_centre[2] << std::endl;
            // std::cout << "cube len :" << oct_cube->Xlen() << std::endl;
            // std::cout << "cube id: " << oct_cube->GetCubeId() << std::endl;
            std::vector<double> half_edge_lens = {oct_cube->Xlen() / 2., oct_cube->Ylen() / 2., oct_cube->Zlen() / 2.};
            for (size_t i = 1; i <= 6; i++) {
                std::vector<OctCube*> region_cube((1 << static_config::cellsCntEachEdgeOnCube2) * (1 << static_config::cellsCntEachEdgeOnCube2));
                for (size_t k{0U}; k < region_cube.size(); k++) {
                    region_cube[k] = nullptr;
                }
                // std::cout << "init region_cube " << std::endl;
                auto neigh_cubes = oct_cube->GetNeiOctCubesOnOneSurface(id2cube, i - 1);
                // std::cout << "neigh_cubes get" << std::endl;
                double ori1_l = cube_centre[step_metric[i - 1][0]] - half_edge_lens[step_metric[i - 1][0]];
                double ori1_r = cube_centre[step_metric[i - 1][0]] + half_edge_lens[step_metric[i - 1][0]];
                double ori2_l = cube_centre[step_metric[i - 1][1]] - half_edge_lens[step_metric[i - 1][1]];
                double ori2_r = cube_centre[step_metric[i - 1][1]] + half_edge_lens[step_metric[i - 1][1]];
                size_t seg_cnt = (1 << static_config::cellsCntEachEdgeOnCube2);
                double seg_len_1 = half_edge_lens[step_metric[i - 1][0]] * 2. / (1 << static_config::cellsCntEachEdgeOnCube2);
                double seg_len_2 = half_edge_lens[step_metric[i - 1][1]] * 2. / (1 << static_config::cellsCntEachEdgeOnCube2);
                // std::cout << "neigh cubes size: " << neigh_cubes.size() << std::endl;
                for (auto neigh_cube : neigh_cubes) {
                    size_t neigh_level = neigh_cube->GetCubeLevel();
                    // std::cout << "neigh rela = " << neigh_cube->GetRelaToNode() << std::endl;
                    // std::cout << "neigh node pt = " << neigh_cube->GetOctNode()->CentralPoint().ToString() << std::endl;
                    // std::cout << "cur cube level = " << cur_cube_level << " ,neigh_level = " << neigh_level << std::endl;
                    // std::cout << "neigh len = " << neigh_cube->Xlen() << std::endl;
                    std::vector<double> neigh_centre = {neigh_cube->CentralPoint().x, neigh_cube->CentralPoint().y, neigh_cube->CentralPoint().z};
                    // std::cout << "neigh_centre: " << neigh_centre[0] << ", " << neigh_centre[1] << ", " << neigh_centre[2] << std::endl;
                    double neigh_ori1 = neigh_centre[step_metric[i - 1][0]];
                    double neigh_ori2 = neigh_centre[step_metric[i - 1][1]];
                    if (cur_cube_level + static_config::cellsCntEachEdgeOnCube2 <= neigh_level) {
                        // std::cout << "small neigh cubes" << std::endl;
                        // std::cout << "seg_len = " << seg_len_1 << std::endl;
                        // std::cout << "ori1_l = " << ori1_l << std::endl;
                        // std::cout << "neigh_ori1 = " << neigh_ori1 << "  ,neigh_ori2 = " << neigh_ori2 << std::endl;
                        size_t dir1 = (neigh_ori1 - ori1_l) / seg_len_1; 
                        size_t dir2 = (neigh_ori2 - ori2_l) / seg_len_2;
                        // std::cout << "dir1 = " << dir1 << ", dir2 = " << dir2 << std::endl;
                        size_t region_idx = dir1 * seg_cnt + dir2;
                        // std::cout << "region cube size: " << region_cube.size() << std::endl;
                        // std::cout << "region_idx = " << region_idx << std::endl;
                        // std::cout << "region cube size: " << region_cube.size() << std::endl;
                        if (region_cube[region_idx] == nullptr) {
                            // std::cout << "case 1" << std::endl;
                            region_cube[region_idx] = neigh_cube;
                        } else if (region_cube[region_idx]->GetCubeLevel() > neigh_cube->GetCubeLevel()) {
                            // std::cout << "case 2" << std::endl;
                            region_cube[region_idx] = neigh_cube;
                        }
                        // std::cout << "wha ??" << std::endl;
                    } else {
                        // std::cout << "large neigh cubes" << std::endl;
                        std::vector<double> neigh_half_elen{neigh_cube->Xlen() / 2., neigh_cube->Ylen() / 2., neigh_cube->Zlen() / 2.};
                        size_t dir1l = std::max(neigh_ori1 - neigh_half_elen[step_metric[i - 1][0]] - ori1_l + 1e-2, 0.) / seg_len_1;
                        size_t dir1r = std::max(neigh_ori1 + neigh_half_elen[step_metric[i - 1][0]] - ori1_l + 1e-2, 0.) / seg_len_1;
                        size_t dir2l = std::max(neigh_ori2 - neigh_half_elen[step_metric[i - 1][1]] - ori2_l + 1e-2, 0.) / seg_len_2;
                        size_t dir2r = std::max(neigh_ori2 + neigh_half_elen[step_metric[i - 1][1]] - ori2_l + 1e-2, 0.) / seg_len_2;
                        // std::cout << "dir1l = " << dir1l << " , dir1r = " << dir1r << std::endl;
                        // std::cout << "dir2l = " << dir2l << " , dir2r = " << dir2r << std::endl;
                        dir1r = std::min(dir1r, seg_cnt);
                        dir2r = std::min(dir2r, seg_cnt);
                        for (size_t pt1 = dir1l; pt1 < dir1r; pt1++) {
                            for (size_t pt2 = dir2l; pt2 < dir2r; pt2++) {
                                size_t region_idx = pt1 * seg_cnt + pt2;
                                if (region_cube[region_idx] == nullptr) {
                                    region_cube[region_idx] = neigh_cube;
                                } else if (region_cube[region_idx]->GetCubeLevel() > neigh_cube->GetCubeLevel()) {
                                    region_cube[region_idx] = neigh_cube;
                                }
                            }
                        }
                    }
                }
                // std::cout << "travel region cubes, region size: " << region_cube.size() << std::endl;
                for (size_t j{0U}; j < region_cube.size(); j++) {
                //for(auto neigh_cubes: neigh_cubes_regions) {
                    auto cube_selected = region_cube[j];
                    if (cube_selected == nullptr) {
                        continue;
                        // std::cout << "waaaarning !!" << std::endl;
                    }
                    // std::cout << "check selected cube" << cube_selected->CentralPoint().ToString() << std::endl;
                    // auto cube_selected = neigh_cubes.front();
                    Vertice sur_ver(cube_selected, i, (i << (static_config::cellsCntEachEdgeOnCube2 * static_config::cellsCntEachEdgeOnCube2)) + j);
                    if (cube_selected->GetCubeType() == CubeType::CUBE_ONE) {
                        // std::cout << "CUBE ONE" << std::endl;
                        sur_ver.SetIdx(0);
                        sur_ver.SetOri(0);
                        sur_ver.geom_pt = cube_selected->CentralPoint();
                    } else {
                        // std::cout << "CUBE TWO" << std::endl;
                        // TODO
                        size_t cube_level = oct_cube->GetCubeLevel();
                        size_t selected_cube_level = cube_selected->GetCubeLevel();
                        // std::cout << "wha ??" << std::endl;
                        if (selected_cube_level - cube_level <= static_config::cellsCntEachEdgeOnCube2) {
                            // std::cout << "i = " << i << std::endl;
                            // std::cout << "wh !!" << std::endl;
                            sur_ver.geom_pt = oct_cube->CentralPoint() + Point3D(oct_cube->Xlen() * xx[i] / 2., oct_cube->Ylen() * yy[i] / 2., oct_cube->Zlen() * zz[i] / 2.) + CalRegionDiff(i, j);
                            sur_ver.SetOri(7 - i);
                            sur_ver.SetIdx(j + oct_cube->GetCubeId());
                            // std::cout << "now specail case , checking !!!" << std::endl;
                            // std::cout << "oct_cube len: " << oct_cube->Xlen() << std::endl;
                            // std::cout << "central point : " << oct_cube->CentralPoint().ToString() << std::endl;
                            // std::cout << "base bias : " << Point3D(oct_cube->Xlen() * xx[i], oct_cube->Ylen() * yy[i], oct_cube->Zlen() * zz[i]).ToString() << std::endl;
                            // std::cout << "region diff : " << CalRegionDiff(i,j).ToString() << std::endl;
                            /// / std::cout << "geom pt : " << sur_ver.geom_pt.ToString() << std::endl;
                            // std::cout << "selecte cube : " << cube_selected->CentralPoint().ToString() << std::endl;
                        } else {
                            sur_ver.SetOri(7 - i);
                            sur_ver.SetIdx(j + oct_cube->GetCubeId());
                            size_t ori_inv = 7 - i;
                            // std::cout << "ori inv = " << ori_inv << std::endl;
                             sur_ver.geom_pt = cube_selected->CentralPoint() + Point3D(cube_selected->Xlen() * xx[ori_inv] / 2., cube_selected->Ylen() * yy[ori_inv] / 2., cube_selected->Zlen() * zz[ori_inv] / 2.);
                        }
                    }
                    res.push_back(sur_ver);
                }
            }
        }

        return res;
    } else {
        std::cout << "  ver 2" << std::endl;
        std::vector<double> cube_centre = {oct_cube->CentralPoint().x, oct_cube->CentralPoint().y, oct_cube->CentralPoint().z};
            std::vector<double> half_edge_lens = {oct_cube->Xlen() / 2., oct_cube->Ylen() / 2., oct_cube->Zlen() / 2.};
            for (size_t i = 1; i <= 6; i++) {
                if (i == ori_on_cube) {
                    continue;
                }
                std::vector<OctCube*> region_cube((1 << static_config::cellsCntEachEdgeOnCube2) * (1 << static_config::cellsCntEachEdgeOnCube2));
                for (size_t k{0U}; k < region_cube.size(); k++) {
                    region_cube[k] = nullptr;
                }
                auto neigh_cubes = oct_cube->GetNeiOctCubesOnOneSurface(id2cube, i - 1);
                double ori1_l = cube_centre[step_metric[i - 1][0]] - half_edge_lens[step_metric[i - 1][0]];
                double ori1_r = cube_centre[step_metric[i - 1][0]] + half_edge_lens[step_metric[i - 1][0]];
                double ori2_l = cube_centre[step_metric[i - 1][1]] - half_edge_lens[step_metric[i - 1][1]];
                double ori2_r = cube_centre[step_metric[i - 1][1]] + half_edge_lens[step_metric[i - 1][1]];
                size_t seg_cnt = (1 << static_config::cellsCntEachEdgeOnCube2);
                double seg_len_1 = half_edge_lens[step_metric[i - 1][0]] * 2. / (1 << static_config::cellsCntEachEdgeOnCube2);
                double seg_len_2 = half_edge_lens[step_metric[i - 1][1]] * 2. / (1 << static_config::cellsCntEachEdgeOnCube2);
                // std::cout << "neigh size: " << neigh_cubes.size() << std::endl;
                for (auto neigh_cube : neigh_cubes) {
                    // std::cout << "neigh cube id: " << neigh_cube->GetCubeId() << std::endl;
                    // std::cout << "neigh cube centre: " << neigh_cube->CentralPoint().ToString() << std::endl;
                    if(neigh_cube == nullptr) {
                        // std::cout << "warning !! nullptr" << std::endl;
                    }
                    // std::cout << "neigh cube level: " << neigh_cube->GetCubeLevel() << std::endl;
                    // std::cout << "neigh cube xlen: " << neigh_cube->Xlen() << std::endl;
                    double  neigh_level = neigh_cube->GetCubeLevel();
                    std::vector<double> neigh_centre = {neigh_cube->CentralPoint().x, neigh_cube->CentralPoint().y, neigh_cube->CentralPoint().z};
                    double neigh_ori1 = cube_centre[step_metric[i - 1][0]];
                    double neigh_ori2 = cube_centre[step_metric[i - 1][1]];
                    if (cur_cube_level + static_config::cellsCntEachEdgeOnCube2 <= neigh_level) {
                        size_t dir1 = (neigh_ori1 - ori1_l) / seg_len_1; 
                        size_t dir2 = (neigh_ori2 - ori2_l) / seg_len_2;
                        // std::cout << "dir1 = " << dir1 << ", dir2 = " << dir2 << std::endl;
                        size_t region_idx = dir1 * seg_cnt + dir2;
                        // std::cout << "region_idx = " << region_idx << std::endl;
                        if (region_cube[region_idx] == nullptr) {
                            region_cube[region_idx] = neigh_cube;
                        } else if (region_cube[region_idx]->GetCubeLevel() > neigh_cube->GetCubeLevel()) {
                            region_cube[region_idx] = neigh_cube;
                        }
                        // std::cout << "whall ??" << std::endl;
                    } else {
                        // std::cout << "large neigh cubes" << std::endl;
                        std::vector<double> neigh_half_elen{neigh_cube->Xlen() / 2., neigh_cube->Ylen() / 2., neigh_cube->Zlen() / 2.};
                        size_t dir1l = std::max(neigh_ori1 - neigh_half_elen[step_metric[i - 1][0]] - ori1_l + 1e-2, 0.) / seg_len_1;
                        size_t dir1r = std::max(neigh_ori1 + neigh_half_elen[step_metric[i - 1][0]] - ori1_l + 1e-2, 0.) / seg_len_1;
                        size_t dir2l = std::max(neigh_ori2 - neigh_half_elen[step_metric[i - 1][1]] - ori2_l + 1e-2, 0.) / seg_len_2;
                        size_t dir2r = std::max(neigh_ori2 + neigh_half_elen[step_metric[i - 1][1]] - ori2_l + 1e-2, 0.) / seg_len_2;
                        dir1r = std::min(dir1r, seg_cnt);
                        dir2r = std::min(dir2r, seg_cnt);
                        for (size_t pt1 = dir1l; pt1 < dir1r; pt1++) {
                            for (size_t pt2 = dir2l; pt2 < dir2r; pt2++) {
                                size_t region_idx = pt1 * seg_cnt + pt2;
                                if (region_cube[region_idx] == nullptr) {
                                    region_cube[region_idx] = neigh_cube;
                                } else if (region_cube[region_idx]->GetCubeLevel() > neigh_cube->GetCubeLevel()) {
                                    region_cube[region_idx] = neigh_cube;
                                }
                            }
                        }
                    }
                    // std::cout << "sss" << std::endl;
                }
                std::cout << "travel region cubes!" << std::endl;
                for (size_t j{0U}; j < region_cube.size(); j++) {
                //for(auto neigh_cubes: neigh_cubes_regions) {
                    auto cube_selected = region_cube[j];
                    if (cube_selected == nullptr) {
                        // std::cout << "warrrrr !!!!!!" << std::endl;
                        continue;
                    }
                    // std::cout << "check cube selcted idx = " << cube_selected->GetCubeId() << std::endl;
                    // auto cube_selected = neigh_cubes.front();
                    Vertice sur_ver(cube_selected, i, (i << (static_config::cellsCntEachEdgeOnCube2 * static_config::cellsCntEachEdgeOnCube2)) + j);
                    if (cube_selected->GetCubeType() == CubeType::CUBE_ONE) {
                        // std::cout << "CUBE ONE" << std::endl;
                        sur_ver.geom_pt = cube_selected->CentralPoint();
                        sur_ver.SetIdx(0);
                        sur_ver.SetOri(0);
                    } else {
                        // std::cout << "CUBE TWO" << std::endl;
                        // TODO
                        size_t cube_level = oct_cube->GetCubeLevel();
                        size_t selected_cube_level = cube_selected->GetCubeLevel();
                        if (selected_cube_level - cube_level <= static_config::cellsCntEachEdgeOnCube2) {
                            // std::cout << "wtf ??" << std::endl;
                            sur_ver.geom_pt = oct_cube->CentralPoint() + Point3D(oct_cube->Xlen() * xx[i] / 2., oct_cube->Ylen() * yy[i] / 2., oct_cube->Zlen() * zz[i] / 2.) + CalRegionDiff(i, j);
                            // std::cout << "gpt check : " << sur_ver.geom_pt.ToString() << std::endl;
                            sur_ver.SetOri(7 - i);
                            sur_ver.SetIdx(j + oct_cube->GetCubeId());
                        } else {
                            // std::cout << "ERERER" << std::endl;
                            size_t ori_inv = 7 - i;
                            sur_ver.geom_pt = cube_selected->CentralPoint() + Point3D(cube_selected->Xlen() * xx[ori_inv] / 2., cube_selected->Ylen() * yy[ori_inv] / 2., cube_selected->Zlen() * zz[ori_inv] / 2.);
                            sur_ver.SetOri(7 - i);
                            sur_ver.SetIdx(j + oct_cube->GetCubeId());
                        }
                    }
                    res.push_back(sur_ver);
                }
            }
        Vertice cver(oct_cube, 0, 0);
        cver.geom_pt = oct_cube->CentralPoint();
        res.push_back(cver);
        return res;
    }
}

/*
std::vector<std::vector<Vertice>> GetVerticesOnOneSurface(size_t ori, int partition_cnt) {
    std::vector<std::vector<Vertice>> res;
    if (ori_on_cube)
}
*/

} // namespace global_path_search
