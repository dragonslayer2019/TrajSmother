#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>
#include <sys/time.h>
#include "BlockMatrix.h"
#include "MPC.h"
#include "FunctionG.h"
#include "json.hpp"
#include <decomp_ros_utils/data_ros_utils.h>
#include <decomp_util/ellipsoid_decomp.h>
#include <nav_msgs/Path.h>
using namespace std;

/*
const int SizeX = 4, SizeU = 1;
const int SizeEqx = 4, SizeEqu = 1;
const int NumEllx = 1, NumEllu = 0;
const int SizeG = SizeEqx + NumEllx + SizeEqu + NumEllu;
const int SizeYx = 6, SizeYu = 1;
const int HorizonNum = 49;
const float pi = M_PI;
int SizeEllx[NumEllx] = {2}, SizeEllu[NumEllu] = {};

typedef Eigen::Matrix<float, SizeX, SizeX> MatrixX;
typedef Eigen::Matrix<float, SizeU, SizeU> MatrixU;
typedef Eigen::Matrix<float, SizeX, SizeU> MatrixB;
typedef Eigen::Matrix<float, SizeX, 1> VectorX;
typedef Eigen::Matrix<float, SizeU, 1> VectorU;
typedef Eigen::Matrix<float, SizeYx, SizeX> MatrixHx;
typedef Eigen::Matrix<float, SizeYu, SizeU> MatrixHu;
VectorX x_init;
MatrixX Ai, Qi;
MatrixB Bi;
MatrixHx Hxi;
MatrixU Ri, Hui;
VectorX ci, Li;
VectorU Wi;
*/


// MPC步数、状态变量个数、控制变量个数、5个切平面+6个bbox约束+1个椭球约束、1个最小化曲率+1个防止曲率过大的约束
const int SizeX = 6, SizeU = 3;
const int SizeEqx = 11, SizeEqu = 0;
const int NumEllx = 1, NumEllu = 2;
const int SizeG = SizeEqx + NumEllx + SizeEqu + NumEllu;
const int SizeYx = 14, SizeYu = 4;
const int HorizonNum = 49;
const float pi = M_PI;
int KAcc = 150;
const float inf = 1e5;
int SizeEllx[NumEllx] = {3}, SizeEllu[2] = {2,2};

typedef Eigen::Matrix<float, SizeX, SizeX> MatrixX;
typedef Eigen::Matrix<float, SizeU, SizeU> MatrixU;
typedef Eigen::Matrix<float, SizeX, SizeU> MatrixB;
typedef Eigen::Matrix<float, SizeX, 1> VectorX;
typedef Eigen::Matrix<float, SizeU, 1> VectorU;
typedef Eigen::Matrix<float, SizeYx, SizeX> MatrixHx;
typedef Eigen::Matrix<float, SizeYu, SizeU> MatrixHu;


bool isPointOnSegment(const Eigen::Vector3f& A, const Eigen::Vector3f& B, const Eigen::Vector3f& P) {
    Eigen::Vector3f AB = B - A;
    Eigen::Vector3f AP = P - A;
    
    // 计算叉积
    Eigen::Vector3f crossProduct = AB.cross(AP);
    
    // 判断叉积是否为零向量
    if (!crossProduct.isZero(1e-1)) { // 允许一定的误差范围
        return false;
    }
    
    // 计算点积
    float dotProduct = AP.dot(AB);
    float squaredLengthAB = AB.dot(AB);
    
    // 判断点积是否在规定范围内
    if (0 <= dotProduct && dotProduct <= squaredLengthAB) {
        return true;
    } else {
        return false;
    }
}


void solveunit3D(vector<float> dt, vector<float> Px, vector<float> Py, vector<float> Pz, vector<float> v_norm,
                 vector<Eigen::Matrix<float, 3, 3>> Rk, vector<float> lamb1, vector<float> lamb2, vector<float> lamb3, vector<float> lamb4,
                 vector<float> lamb5, vector<vector<Hyperplane<3>>> CorridorP, std::array<Eigen::Matrix<float, SizeYx - SizeEqx, 1>, HorizonNum + 1> new_centerX,
                 std::array<Eigen::Matrix<float, SizeYu - SizeEqu, 1>, HorizonNum + 1> new_centerU, std::array<Eigen::Matrix<float, 3, 3>, HorizonNum + 1> elliE, int K = 250) {

    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    std::array<MatrixHx, HorizonNum + 1> Hx;
    std::array<MatrixHu, HorizonNum + 1> Hu;

    MatrixX Ai, Qi;
    MatrixHx Hxi;
    MatrixB Bi;
    MatrixU Ri;
    MatrixHu Hui;
    VectorX ci, Li;
    VectorU Wi;

    std::array<std::array<FunctionG<float>, SizeG>, HorizonNum + 1> g;
    std::vector<vector<float>> DistC(HorizonNum+1, std::vector<float>(HorizonNum+1));
    std::vector<float> aaa, bbb, ccc;
    std::vector<float> ppp;
    // 状态变量(px, py, pz, vx, vy, vz)^T
    // 控制输入(mu1, mu2, mu3)^T
    VectorX x_init;
    x_init << 5, 9.5, 0.5, 3, 3, 3;



    for(int i = 0; i <= HorizonNum; ++i) {
        // 此处为m个切平面约束+一个椭球约束
        // 计算优化路点到切平面之间距离的结果中的常数转移到gxi函数中               
        int M = CorridorP[i].size();
        for (int j = 0; j < M; ++j) {
            // m行到切平面距离   
            Hxi(j,0) = CorridorP[i][j].n_.x();
            Hxi(j,1) = CorridorP[i][j].n_.y();
            Hxi(j,2) =  CorridorP[i][j].n_.z();
            Hxi(j,3) = 0;
            Hxi(j,4) = 0;
            Hxi(j,5) = 0;
            // 计算点到切平面距离中的const向量
            DistC[i][j] = -CorridorP[i][j].n_.x()*CorridorP[i][j].p_.x() - CorridorP[i][j].n_.y()*CorridorP[i][j].p_.y() - CorridorP[i][j].n_.z()*CorridorP[i][j].p_.z();
            std::cout << "DistC[i][j]" <<  i << "," << j << ":" <<  DistC[i][j] << std::endl;      
        }
                
        // 1个到椭球圆心距离约束
        for (int k = 0; k < 3; k++) {
            Hxi(M+k, 0) = elliE[i](k, 0);
            Hxi(M+k, 1) = elliE[i](k, 1);
            Hxi(M+k, 2) = elliE[i](k, 2);
            Hxi(M+k, 3) = 0;
            Hxi(M+k, 4) = 0;
            Hxi(M+k, 5) = 0;
        }

        // 此处为一个曲率平方约束+一个曲率罚函数   ck = sqrt(mu2^2+mu3^2)
        Hui << 0, 1, 0,
               0, 0, 1,
               0, 1, 0,
               0, 0, 1;
        Hx[i] = Hxi;
        Hu[i] = Hui;
    }

    for(int i = 0; i < HorizonNum; ++i) {
        // 6*6
        Ai << 1, 0, 0, dt[i], 0, 0,
              0, 1, 0, 0,  dt[i], 0,
              0, 0, 1 , 0, 0, dt[i],
              0, 0, 0 , 1, 0, 0,
              0, 0, 0 , 0, 1, 0,
              0, 0, 0 , 0, 0, 1;
        // 6*3
        Bi << 0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](0,0), 0.5*v_norm[i]*v_norm[i]*dt[i]*Rk[i](0,1), 0.5*v_norm[i]*v_norm[i]*dt[i]*Rk[i](0,2),
              0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](1,0), 0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](1,1), 0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](1,2),
              0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](2,0), 0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](2,1), 0.5*v_norm[i]*v_norm[i]*dt[i]*dt[i]*Rk[i](2,2),
              v_norm[i]*v_norm[i]*dt[i]*Rk[i](0,0), v_norm[i]*v_norm[i]*dt[i]*Rk[i](0,1), 0.5*v_norm[i]*v_norm[i]*dt[i]*Rk[i](0,2),
              v_norm[i]*v_norm[i]*dt[i]*Rk[i](1,0), v_norm[i]*v_norm[i]*dt[i]*Rk[i](1,1), 0.5*v_norm[i]*v_norm[i]*dt[i]*Rk[i](1,2),
              v_norm[i]*v_norm[i]*dt[i]*Rk[i](2,0), v_norm[i]*v_norm[i]*dt[i]*Rk[i](2,1), 0.5*v_norm[i]*v_norm[i]*dt[i]*Rk[i](2,2);
        // 6*1
        ci << 0, 0, 0, 0, 0, 0;

        A[i] = Ai; B[i] = Bi; c[i] = ci;
    }

    for(int i = 0; i <= HorizonNum; ++i) {
        // 追参考位置
        Qi << lamb1[i], 0, 0, 0, 0, 0,
              0, lamb1[i], 0, 0, 0, 0,
              0, 0, lamb1[i], 0, 0, 0,
              0, 0, 0, 0.01, 0, 0,
              0, 0, 0, 0, 0.01, 0,
              0, 0, 0, 0, 0, 0.01;
        Li << -Px[i] * lamb1[i],
              -Py[i] * lamb1[i],
              -Pz[i] * lamb1[i],
              0,
              0,
              0;
        // 限制纵向加速度
        Ri << lamb2[i], 0, 0,
              0, 0.1, 0,
              0, 0, 0.1;
        Wi << 0,
              0,
              0;
        Q[i] = Qi; R[i] = Ri; L[i] = Li; W[i] = Wi;
    }

    // 设置状态变量相关的罚函数形状
    for(int i = 0;i <= 3; ++i) {
        aaa.push_back(0);
        bbb.push_back(0);
        ccc.push_back(0);
        ppp.push_back(0);
    }
    for(int i = 0; i <= HorizonNum; ++i) {
        // 凸走廊切平面罚函数
        int M = CorridorP[i].size();
        for (int j = 0; j < M; j++) {
            ppp[0] = -inf;ppp[1] = 0.2;ppp[2] = 1;ppp[3] = inf;
            aaa[0] = 10;bbb[0] = -22 + 20*DistC[i][j];ccc[0] = 10*DistC[i][j]*DistC[i][j] - 22*DistC[i][j] + 5;
            aaa[1] = 1;bbb[1] = 2*DistC[i][j] - 2.45;ccc[1] = DistC[i][j]*DistC[i][j] - 2.45*DistC[i][j] + 1.45;
            aaa[2] = 0;bbb[2] = 0;ccc[2] = 0;
            g[i][j].AddQuadratic(3, aaa, bbb, ccc, ppp);
        }
        
        // 凸走廊椭球罚函数
        ppp[0] = -inf;ppp[1] = 0.5;ppp[2] = 1.0;ppp[3] = inf;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = 1; bbb[1] = 0; ccc[1] = -0.25; 
        aaa[2] = 10; bbb[2] = 0; ccc[2] = -9.25;
        g[i][M].AddQuadratic(3, aaa, bbb, ccc, ppp);

        // 曲率平方约束
        ppp[0] = -inf;ppp[1] = 0;ppp[2] = inf;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = 1; bbb[1] = 0; ccc[1] = 0;
        g[i][M+1].AddQuadratic(2, aaa, bbb, ccc, ppp);

        // 曲率罚函数
        ppp[0] = -inf;ppp[1] = 0.5;ppp[2] = 2.0;ppp[3] = inf;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = 1; bbb[1] = 0; ccc[1] = -0.25;
        aaa[2] = 10; bbb[2] = 0; ccc[2] = -36.25;
        g[i][M+2].AddQuadratic(3, aaa, bbb, ccc, ppp);
    }

    std::cout << "in solveunit3D" << std::endl;
    // Q, R, L, W, A, B, c, x_init, Hx, Hu, new_centerX, new_centerU, g, 100, KAcc, K
    // std::cout << "print A" << std::endl;
    // for (auto& a : A) {
    //     std::cout << a << std::endl;
    // }
    // std::cout << "print B" << std::endl;
    // for (auto& b : B) {
    //     std::cout << b << std::endl;
    // }
    // std::cout << "print Hx" << std::endl;
    // for (auto& hxi : Hx) {
    //     std::cout << hxi << std::endl;
    // }
    std::cout << "print Hu" << std::endl;
    for (auto& hui : Hu) {
        std::cout << hui << std::endl;
    }
    std::cout << "print centerX" << std::endl;
    for (auto& centerX : new_centerX) {
        std::cout << centerX.transpose() << std::endl;
    }
    std::cout << "print centerU" << std::endl;
    for (auto& centerU : new_centerU) {
        std::cout << centerU.transpose() << std::endl;
    }    
    

    MPC_ADMMSolver<float, HorizonNum, SizeX, SizeU, SizeYx, SizeYu, SizeEqx, SizeEqu, NumEllx, NumEllu, SizeG, SizeEllx, SizeEllu> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, new_centerX, new_centerU, g, 100, KAcc, K);
    // MPC_ADMMSolver<float, HorizonNum, SizeX, SizeU, SizeYx, SizeYu, SizeEqx, SizeEqu, NumEllx, NumEllu, SizeG, SizeEllx, SizeEllu> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, new_centerX, new_centerU, g, 100, KAcc, K);
    // BlockVector<float, HorizonNum + 1, SizeX + SizeU> res;
    // mpc.solve(res);
    BlockVector<float, HorizonNum + 1, SizeX + SizeU> res = mpc.solve();
     std::cout << "finish solveunit3D" << std::endl;
    

    using json = nlohmann::json;
    json j;
    json jx = json::array();
    json jy = json::array();
    json jz = json::array();
    json jvx = json::array();
    json jvy = json::array();
    json jvz = json::array();
    json ju = json::array();
    // vector<float> t(HorizonNum), x(HorizonNum), v(HorizonNum), a(HorizonNum), u(HorizonNum);
    for(int i = 0; i <= HorizonNum; ++i) {
        jx.push_back(res.v[i](0, 0));
        jy.push_back(res.v[i](1, 0));
        jz.push_back(res.v[i](2, 0));
        jvx.push_back(res.v[i](3, 0));
        jvy.push_back(res.v[i](4, 0));
        jvz.push_back(res.v[i](5, 0));
        ju.push_back(res.v[i](6, 0));
    }
    std::cout << "px: ";
    for(int i = 0; i <= HorizonNum; ++i) {
        std::cout << " " <<res.v[i](0, 0) << std::endl;
    }
    j["x"] = jx;
    j["y"] = jy;
    j["z"] = jz;
    j["vx"] = jvx;
    j["vy"] = jvy;
    j["vz"] = jvz;
    j["u"] = ju;
    j["Px"] = Px;
    ofstream out("test1.out");
    if(out.is_open()) {
        out << j.dump(4) << endl;
        out.close();
    }
    
    return;
}

int solveMpc(vec_E<Polyhedron<3>> mpc_polyhedrons, std::array<Eigen::Matrix<float, SizeYx - SizeEqx, 1>, HorizonNum + 1> new_centerX, std::array<Eigen::Matrix<float, SizeYu - SizeEqu, 1>, HorizonNum + 1> new_centerU, std::array<Eigen::Matrix<float, 3, 3>, HorizonNum + 1> elliE, vector<float> dt, vector<Eigen::Vector3f> ref_points, vector<float> v_norm, vector<Eigen::Matrix<float, 3, 3>> Rk) {
    float q=0.35, st=0, wei=10, weig=50; int K=250;
    K = 3;
    struct timeval T1,T2;
    float timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<float> Px, Py, Pz, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5;
    vector<vector<float>> center;
    vector<vector<Hyperplane<3>>> CorridorP;//分段切平面约束
    CorridorP.resize(HorizonNum+1);
    Px.resize(HorizonNum+1);
    Py.resize(HorizonNum+1);
    Pz.resize(HorizonNum+1);

    // std::cout << "here1" << std::endl;
    // std::cout << "ref_points size: " << ref_points.size() << std::endl;
    // std::cout << "mpc_polyhedrons size: " << mpc_polyhedrons.size() << std::endl;
    for(int i = 0;i <= HorizonNum; ++i) {
        std::cout << "loop: " << i << std::endl;
        // 切平面约束与椭球约束
        int plane_size = mpc_polyhedrons[i].hyperplanes().size();
        CorridorP[i].resize(plane_size);
        for (int j = 0; j < plane_size; j++) {
            CorridorP[i][j] = mpc_polyhedrons[i].hyperplanes()[j];
        }
        std::cout << "loop: " << i << std::endl;
        // 用ref_points为Px、Py、Pz赋值
        Px[i] = ref_points[i].x();
        Py[i] = ref_points[i].y();
        Pz[i] = ref_points[i].z();
        // 代价权重赋值
        lamb1.push_back(1);
        lamb2.push_back(1);
        lamb3.push_back(1);
        lamb4.push_back(1);
        lamb5.push_back(1);
        std::cout << "loop: " << i << std::endl;
    }
    std::cout << "start solveunit3D" << std::endl;
    solveunit3D(dt, Px, Py, Pz, v_norm, Rk, lamb1, lamb2, lamb3, lamb4, lamb5, CorridorP, new_centerX,  new_centerU, elliE, K);
    gettimeofday(&T2,NULL);
    timeuse = (T2.tv_sec - T1.tv_sec) + (float)(T2.tv_usec - T1.tv_usec)/1000000.0;
    std::cout<<"time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
    return 0;
}



// 计算两个三维点之间的欧氏距离
float distance(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2) {
    return sqrt(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2) + pow(p1(2) - p2(2), 2));
}

// 插值路径点，步长为max_length
vector<Eigen::Vector3f> interpolatePoints(const vector<Eigen::Vector3f>& points, float max_length) {
    vector<Eigen::Vector3f> interpolated_points;
    interpolated_points.push_back(points[0]);

    for (size_t i = 1; i < points.size(); ++i) {
        Eigen::Vector3f p1 = points[i - 1];
        Eigen::Vector3f p2 = points[i];
        float segment_length = distance(p1, p2);
        size_t num_segments = static_cast<size_t>(ceil(segment_length / max_length));
        Eigen::Vector3f direction(p2(0) - p1(0), p2(1) - p1(1), p2(2) - p1(2));
        direction.normalize();

        for (size_t j = 1; j < num_segments; ++j) {
            Eigen::Vector3f new_point = {
                p1.x() + direction.x() * max_length * j,
                p1.y() + direction.y() * max_length * j,
                p1.z() + direction.z() * max_length * j
            };
            interpolated_points.push_back(new_point);
        }
        interpolated_points.push_back(p2);  // Ensure the last point of the segment is added
    }
    return interpolated_points;
}

// 计算曲率
vector<float> computeCurvature(const vector<Eigen::Vector3f>& points, float resample_dist) {
    int n = points.size();
    vector<float> curvature(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        Eigen::Vector3f v1(points[i](0) - points[i - 1](0), points[i](1) - points[i - 1](1), points[i](2) - points[i - 1](2));
        Eigen::Vector3f v2(points[i + 1](0) - points[i](0), points[i + 1](1) - points[i](1), points[i + 1](2) - points[i](2));
        float angle = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
        curvature[i] = angle / resample_dist;
    }
    return curvature;
}

// 高斯平滑
vector<float> smoothCurvature(const vector<float>& curvature, int window_size) {
    int n = curvature.size();
    vector<float> smoothed_curvature(n, 0.0);
    int half_window = window_size / 2;
    for (int i = 0; i < n; ++i) {
        float sum = 0.0;
        int count = 0;
        for (int j = max(0, i - half_window); j <= min(n - 1, i + half_window); ++j) {
            sum += curvature[j];
            ++count;
        }
        smoothed_curvature[i] = sum / count;
    }
    return smoothed_curvature;
}

// 重新采样路径点
vector<Eigen::Vector3f> resamplePath(const vector<Eigen::Vector3f>& points, const vector<float>& smoothed_curvature, float max_length, float max_angle) {
    vector<Eigen::Vector3f> resampled_points;
    resampled_points.push_back(points[0]);
    float cum_length = 0.0;
    float cum_angle = 0.0;

    for (size_t i = 1; i < points.size(); ++i) {
        float segment_length = distance(points[i], points[i - 1]);
        cum_length += segment_length;

        if (i > 1 && i <= smoothed_curvature.size()) {
            float delta_angle = smoothed_curvature[i] * segment_length;
            cum_angle += delta_angle;
        }

        if (cum_length > max_length || cum_angle > (max_angle * M_PI / 180.0)) {
            Eigen::Vector3f direction(points[i](0) - points[i - 1](0), points[i](1) - points[i - 1](1), points[i](2) - points[i - 1](2));
            direction.normalize();
            Eigen::Vector3f new_point = {
                points[i].x() - direction.x() * (cum_length - max_length),
                points[i].y() - direction.y() * (cum_length - max_length),
                points[i].z() - direction.z() * (cum_length - max_length)
            };
            resampled_points.push_back(new_point);
            cum_length = distance(points[i], new_point);
            cum_angle = smoothed_curvature[i] * cum_length;
        }
    }
    resampled_points.push_back(points.back());
    return resampled_points;
}

vector<Eigen::Vector3f> get_sample_point(vector<Eigen::Vector3f>& path) {
    // 输入路径点
    // vector<Eigen::Vector3f> points = {
    //     {5, 9.5, 0.5}, {13, 11.5, 3.0}, {15, 9.5, 1.5}, {14, 5, 2.5}
    // };

    // 步长
    float max_length = 0.5;

    // 插值路径点，确保每段长度不超过max_length
    vector<Eigen::Vector3f> interpolated_points = interpolatePoints(path, max_length);

    // 计算曲率
    vector<float> curvature = computeCurvature(interpolated_points, max_length);

    // 平滑曲率
    vector<float> smoothed_curvature = smoothCurvature(curvature, 8);

    // 重新采样路径点
    float max_angle = 20.0;
    vector<Eigen::Vector3f> resampled_points = resamplePath(interpolated_points, smoothed_curvature, max_length, max_angle);

    // 输出重新采样后的路径点

    
    // 标准化参考路径点数为1 + HorizonNum
    vector<Eigen::Vector3f> ref_points(HorizonNum + 1);
    int back_id = resampled_points.size() - 1;
    if (back_id < HorizonNum) {
        std::copy(resampled_points.begin(), resampled_points.end(), ref_points.begin());
        // 计算最后的path方向
        Eigen::Vector3f direction = (resampled_points[back_id] - resampled_points[back_id - 1]).normalized();

        // std::cout << "back id: " << back_id << "dir: " << direction << std::endl;
        // 延长point
        for (int i = 1; i < HorizonNum+1-back_id; i++) {
            Eigen::Vector3f new_point = resampled_points[back_id] + direction * i * 0.5f;
            ref_points[back_id + i] = new_point;
        }
    } else {
        std::cout << "path length enough" << std::endl;
        std::copy(resampled_points.begin(), resampled_points.begin() + HorizonNum + 1, ref_points.begin());
    }
    // std::cout << "Resampled Path Points:" << std::endl;
    // for (const auto& point : ref_points) {
    //     std::cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    // }
    std::cout << "ref_points size: " << ref_points.size() << std::endl;
    for (const auto& point : ref_points) {
        std::cout << point.transpose() << std::endl;
    }
    return ref_points;
}


// 计算变换矩阵的函数
Eigen::Matrix<float, 3, 3> calculateTransformationMatrix(const Eigen::Vector3f& p0, const Eigen::Vector3f& p1, const Eigen::Vector3f& p2) {
    // 计算向量a
    Eigen::Vector3f a = p2 - p1;
    a.normalize();

    // 固定一个向量b的xy分量
    Eigen::Vector3f b(1, 1, 0);
    b.z() = - (a.x() + a.y()) / a.z();

    b.normalize();

    // 计算第三个正交向量
    Eigen::Vector3f c = a.cross(b);
    c.normalize();

    // 构造正交单位矩阵
    Eigen::Matrix3f transformationMatrix;
    transformationMatrix.col(0) = a;  // x轴
    transformationMatrix.col(1) = b;  // y轴
    transformationMatrix.col(2) = c;  // z轴

    return transformationMatrix;
}


void get_param(const vector<Eigen::Vector3f>& ref_points, vector<float>& v_norm, vector<Eigen::Matrix<float, 3, 3>>& Rk, vector<float>& dt) {
  float v_max = 5;
  float s_max = 0.5;
  v_norm.resize(HorizonNum);
  dt.resize(HorizonNum);
  Rk.resize(HorizonNum);
  Eigen::Vector3f p0;
  Eigen::Vector3f p1;
  Eigen::Vector3f p2;

  for (int i = 0; i < HorizonNum; ++i) {
    
    float dis = distance(ref_points[i], ref_points[i+1]);
    if (fabs(dis) < 0.01) {
        std::cout << "hard code v_norm!!!!!!!!!!" << std::endl;
        v_norm[i] = 1.0;
    } else {
        // v_norm[i] = dis/s_max*v_max;
        v_norm[i] = 1.0;
    } 
    
    dt[i] = (dis/v_norm[i]);
    // std::cout << " dis/ v_norm/ dt " << i << " : " <<  dis << ", " << v_norm[i] << ", " << dt[i] << std::endl;

    // Eigen::Matrix3f  
    // 计算Rk[i]，使得其i第一列列向量与dir平行，其余两列单位列向量与其组成正交矩阵
    if (i == 0) {
        p0 = ref_points[0];
        p1 = ref_points[1];
        p2 = ref_points[2];
        Rk[i] = calculateTransformationMatrix(p0, p1, p2);
    } else if (i == HorizonNum - 1) {
        p0 = ref_points[HorizonNum - 2];
        p1 = ref_points[HorizonNum - 1];
        p2 = ref_points[HorizonNum];
        Rk[i] = calculateTransformationMatrix(p0, p1, p2);
    } else {
        p0 = ref_points[i-1];
        p1 = ref_points[i];
        p2 = ref_points[i+1];
        Rk[i] = calculateTransformationMatrix(p0, p1, p2);
    }
    
      
  }

}
