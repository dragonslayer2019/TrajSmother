#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <limits>
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
#include <algorithm>
#include <numeric>
#include <unsupported/Eigen/Splines>
#include "presteps.h"
using namespace std;

// MPC steps, number of state variables, number of control variables,
// 5 cutting planes + 4 bbox constraints + 1 ellipsoid constraint,
// 1 constraint to minimize curvature + 1 constraint to prevent excessive curvature
const int SizeX2d = 4, SizeU2d = 2;
const int SizeEqx2d = 9, SizeEqu2d = 0;
const int NumEllx2d = 2, NumEllu2d = 1;
const int SizeG2d = SizeEqx2d + NumEllx2d + SizeEqu2d + NumEllu2d;
const int SizeYx2d = 13, SizeYu2d = 1;
int SizeEllx2d[NumEllx2d] = {2,2}, SizeEllu2d[NumEllu2d] = {1};

typedef Eigen::Matrix<float, SizeX2d, SizeX2d> MatrixX2d;
typedef Eigen::Matrix<float, SizeU2d, SizeU2d> MatrixU2d;
typedef Eigen::Matrix<float, SizeX2d, SizeU2d> MatrixB2d;
typedef Eigen::Matrix<float, SizeX2d, 1> VectorX2d;
typedef Eigen::Matrix<float, SizeU2d, 1> VectorU2d;
typedef Eigen::Matrix<float, SizeYx2d, SizeX2d> MatrixHx2d;
typedef Eigen::Matrix<float, SizeYu2d, SizeU2d> MatrixHu2d;

float distance_2d(const Eigen::Vector2f& p1, const Eigen::Vector2f& p2) {
    return sqrt(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2));
}

bool isPointOnSegment_2d(const Eigen::Vector2f& A, const Eigen::Vector2f& B, const Eigen::Vector2f& P) {
    Eigen::Vector2f AB = B - A;
    Eigen::Vector2f AP = P - A;
    
    // Calculate the cross product
    // Eigen::Vector2f crossProduct = AB.cross(AP);
    float crossProduct = AB.x() * AP.y() - AB.y() * AP.x();
    
    if (std::fabs(crossProduct) > 0.1) { // 允许一定的误差范围
        return false;
    }
    
    // Calculating the dot product
    float dotProduct = AP.dot(AB);
    float squaredLengthAB = AB.dot(AB);
    
    // Determine whether the dot product is within the specified range
    if (0 <= dotProduct && dotProduct <= squaredLengthAB) {
        return true;
    } else {
        return false;
    }
}

// Interpolate path points with a step size of max_length
vector<Eigen::Vector2f> interpolatePoints_2d(const vector<Eigen::Vector2f>& points, float max_length) {
    vector<Eigen::Vector2f> interpolated_points;
    interpolated_points.push_back(points[0]);

    for (size_t i = 1; i < points.size(); ++i) {
        Eigen::Vector2f p1 = points[i - 1];
        Eigen::Vector2f p2 = points[i];
        float segment_length = distance_2d(p1, p2);
        size_t num_segments = static_cast<size_t>(ceil(segment_length / max_length));
        Eigen::Vector2f direction(p2(0) - p1(0), p2(1) - p1(1));
        direction.normalize();

        for (size_t j = 1; j < num_segments; ++j) {
            Eigen::Vector2f new_point = {
                p1.x() + direction.x() * max_length * j,
                p1.y() + direction.y() * max_length * j
            };
            interpolated_points.push_back(new_point);
        }
        interpolated_points.push_back(p2);  // Ensure the last point of the segment is added
    }
    return interpolated_points;
}

// Calculating curvature
vector<float> computeCurvature_2d(const vector<Eigen::Vector2f>& points, float resample_dist) {
    int n = points.size();
    vector<float> curvature(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        Eigen::Vector2f v1(points[i](0) - points[i - 1](0), points[i](1) - points[i - 1](1));
        Eigen::Vector2f v2(points[i + 1](0) - points[i](0), points[i + 1](1) - points[i](1));
        float angle = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
        curvature[i] = angle / resample_dist;
    }
    return curvature;
}

// Gaussian smoothing
vector<float> smoothCurvature_2d(const vector<float>& curvature, int window_size) {
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

// Resample waypoints
vector<Eigen::Vector2f> resamplePath_2d(const vector<Eigen::Vector2f>& points, const vector<float>& smoothed_curvature, float max_length, float max_angle) {
    vector<Eigen::Vector2f> resampled_points;
    resampled_points.push_back(points[0]);
    float cum_length = 0.0;
    float cum_angle = 0.0;

    for (size_t i = 1; i < points.size(); ++i) {
        float segment_length = distance_2d(points[i], points[i - 1]);
        cum_length += segment_length;

        if (i > 1 && i <= smoothed_curvature.size()) {
            float delta_angle = smoothed_curvature[i] * segment_length;
            cum_angle += delta_angle;
        }

        if (cum_length > max_length || cum_angle > (max_angle * M_PI / 180.0)) {
            Eigen::Vector2f direction(points[i](0) - points[i - 1](0), points[i](1) - points[i - 1](1));
            direction.normalize();
            Eigen::Vector2f new_point = {
                points[i].x() - direction.x() * (cum_length - max_length),
                points[i].y() - direction.y() * (cum_length - max_length)
            };
            resampled_points.push_back(new_point);
            cum_length = distance_2d(points[i], new_point);
            cum_angle = smoothed_curvature[i] * cum_length;
        }
    }
    resampled_points.push_back(points.back());
    return resampled_points;
}

vector<Eigen::Vector2f> get_sample_point_2d(vector<Eigen::Vector2f>& path) {
    // Interpolate path points to ensure that each segment length does not exceed max_length
    vector<Eigen::Vector2f> interpolated_points = interpolatePoints_2d(path, max_length);

    // Calculating curvature
    vector<float> curvature = computeCurvature_2d(interpolated_points, max_length);

    // Smooth Curvature
    vector<float> smoothed_curvature = smoothCurvature_2d(curvature, 8);

    // Resample waypoints
    float max_angle = 20.0;
    vector<Eigen::Vector2f> resampled_points = resamplePath_2d(interpolated_points, smoothed_curvature, max_length, max_angle);

    // optimizePath(resampled_points);

    vector<Eigen::Vector2f> ref_points(HorizonNum + 1);
    int back_id = resampled_points.size() - 1;
    if (back_id < HorizonNum) {
        std::copy(resampled_points.begin(), resampled_points.end(), ref_points.begin());

        Eigen::Vector2f direction = (resampled_points[back_id] - resampled_points[back_id - 1]).normalized();

        for (int i = 1; i < HorizonNum+1-back_id; i++) {
            Eigen::Vector2f new_point = resampled_points[back_id] + direction * i * max_length;
            ref_points[back_id + i] = new_point;
        }
    } else {
        std::cout << "path length enough" << std::endl;
        std::copy(resampled_points.begin(), resampled_points.begin() + HorizonNum + 1, ref_points.begin());
    }
    return ref_points;
}

Eigen::Matrix<float, 2, 2> calculateTransformationMatrix_2d(const Eigen::Vector2f& p1, const Eigen::Vector2f& p2) {
    Eigen::Vector2f a = p2 - p1;
    a.normalize();

    Eigen::Vector2f b(0, 0);
    if (a.y() != 0) {
        b.x() = 1;
        b.y() = - a.x() / a.y();
    } else if (a.x() != 0) {
        b.y() = 1;
        b.x() = - a.y() / a.x();
    } else {
        std::cout << "p2 - p1 is same point, ERROR" << std::endl;
    }

    b.normalize();

    // 构造正交单位矩阵
    Eigen::Matrix2f transformationMatrix;
    transformationMatrix.col(0) = a;  // x轴
    transformationMatrix.col(1) = b;  // y轴

    return transformationMatrix;
}

void get_param_2d(const vector<Eigen::Vector2f>& ref_points, vector<float>& v_norm, vector<Eigen::Matrix<float, 2, 2>>& Rk, vector<float>& dt) {
  float v_max = 5;
//   float s_max = 0.3;
  v_norm.resize(HorizonNum);
  dt.resize(HorizonNum);
  Rk.resize(HorizonNum);
  Eigen::Vector2f p0;
  Eigen::Vector2f p1;
  Eigen::Vector2f p2;

  for (int i = 0; i < HorizonNum; ++i) {
    
    float dis = distance_2d(ref_points[i], ref_points[i+1]);
    if (fabs(dis) < 0.01) {
        std::cout << "hard code v_norm!!!!!!!!!!" << std::endl;
        v_norm[i] = 1.0;
    } else {
        v_norm[i] = 1.0;
    } 
    
    dt[i] = (dis/v_norm[i]);  

    if (i == 0) {
        p0 = ref_points[0];
        p1 = ref_points[1];
        Rk[i] = calculateTransformationMatrix_2d(p0, p1);
    } else if (i == HorizonNum - 1) {
        p0 = ref_points[HorizonNum - 1];
        p1 = ref_points[HorizonNum];
        Rk[i] = calculateTransformationMatrix_2d(p0, p1);
    } else {
        p0 = ref_points[i];
        p1 = ref_points[i+1];
        Rk[i] = calculateTransformationMatrix_2d(p0, p1);
    }
  }
  v_norm.push_back(v_norm.back());
}

/*
void solveunit2D(vector<float> dt, vector<float> Px, vector<float> Py, vector<float> Pz, vector<float> v_norm,
                 vector<Eigen::Matrix<float, 2, 2>> Rk, vector<float> lamb1, vector<float> lamb2, vector<float> lamb3, vector<float> lamb4,
                 vector<float> lamb5, vector<float> lamb6, vector<vector<Hyperplane<2>>> CorridorP, std::array<Eigen::Matrix<float, SizeYx2d - SizeEqx2d, 1>, HorizonNum + 1> new_centerX,
                 std::array<Eigen::Matrix<float, SizeYu2d - SizeEqu2d, 1>, HorizonNum + 1> new_centerU, std::array<Eigen::Matrix<float, 2, 2>, HorizonNum + 1> elliE, BlockVector<float, HorizonNum + 1, SizeX2d + SizeU2d>& res, int K = 250) {

    std::array<MatrixX2d, HorizonNum + 1> Q;
    std::array<MatrixU2d, HorizonNum + 1> R;
    std::array<VectorX2d, HorizonNum + 1> L;
    std::array<VectorU2d, HorizonNum + 1> W;
    std::array<MatrixX2d, HorizonNum> A;
    std::array<MatrixB2d, HorizonNum> B;
    std::array<VectorX2d, HorizonNum> c;
    std::array<MatrixHx2d, HorizonNum + 1> Hx;
    std::array<MatrixHu2d, HorizonNum + 1> Hu;

    MatrixX2d Ai, Qi;
    MatrixHx2d Hxi;
    MatrixB2d Bi;
    MatrixU2d Ri;
    MatrixHu2d Hui;
    VectorX2d ci, Li;
    VectorU2d Wi;

    std::array<std::array<FunctionG<float>, SizeG2d>, HorizonNum + 1> g;
    std::vector<vector<float>> DistC(HorizonNum+1, std::vector<float>(HorizonNum+1));
    std::vector<float> aaa, bbb, ccc;
    std::vector<float> ppp;
    // State variables(px, py, pz, vx, vy, vz)^T
    // Control Input(mu1, mu2, mu3)^T
    VectorX2d x_init;
    x_init <<  Px[0], Py[0], Pz[0], 0.0, 0.0, 0.0;


    for(int i = 0; i <= HorizonNum; ++i) {
        // Here are m tangent plane constraints + one ellipsoid constraint           
        int M = CorridorP[i].size();
        Hxi.setZero();
        for (int j = 0; j < M; ++j) {   
            Hxi(j,0) = -CorridorP[i][j].n_.x();
            Hxi(j,1) = -CorridorP[i][j].n_.y();
            Hxi(j,2) = 0;
            Hxi(j,3) = 0;
            DistC[i][j] = CorridorP[i][j].n_.x()*CorridorP[i][j].p_.x() + CorridorP[i][j].n_.y()*CorridorP[i][j].p_.y() + CorridorP[i][j].n_.z()*CorridorP[i][j].p_.z();
        }
                
        // 1 distance constraint to the ellipsoid center
        for (int k = 0; k < 3; k++) {
            Hxi(SizeEqx+k, 0) = elliE[i](k, 0);
            Hxi(SizeEqx+k, 1) = elliE[i](k, 1);
            Hxi(SizeEqx+k, 2) = elliE[i](k, 2);
            Hxi(SizeEqx+k, 3) = 0;
            Hxi(SizeEqx+k, 4) = 0;
            Hxi(SizeEqx+k, 5) = 0;
        }

        // 1 trust region constraint
        Hxi.block(SizeEqx+3, 0, 3, 6).setZero();
        Hxi(SizeEqx+3, 3) = 1;
        Hxi(SizeEqx+4, 4) = 1;
        Hxi(SizeEqx+5, 5) = 1;

        // Here is a curvature penalty function   ck = sqrt(mu2^2+mu3^2)
        Hui << 0, 1, 0,
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
        // Tracking reference position
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
        // Limiting longitudinal acceleration
        Ri << lamb2[i], 0, 0,
              0, lamb6[i], 0,
              0, 0, lamb6[i];
        Wi << 0,
              0,
              0;
        Q[i] = Qi; R[i] = Ri; L[i] = Li; W[i] = Wi;
    }

    // Set the shape of the penalty function associated with the state variable
    for(int i = 0;i <= 3; ++i) {
        aaa.push_back(0);
        bbb.push_back(0);
        ccc.push_back(0);
        ppp.push_back(0);
    }
    for(int i = 0; i <= HorizonNum; ++i) {
        // Convex Corridor Tangent Plane Penalty Function
        int M = CorridorP[i].size();
        for (int j = 0; j < M; j++) {
            ppp[0] = -inf;ppp[1] = 0.2-DistC[i][j];ppp[2] = 1-DistC[i][j];ppp[3] = inf;
            aaa[0] = lamb5[i]*10;bbb[0] = lamb5[i]*(-22 + 20*DistC[i][j]);ccc[0] = lamb5[i]*(10*DistC[i][j]*DistC[i][j] - 22*DistC[i][j] + 5);
            aaa[1] = lamb5[i];bbb[1] = lamb5[i]*(2*DistC[i][j] - 2.45);ccc[1] = lamb5[i]*(DistC[i][j]*DistC[i][j] - 2.45*DistC[i][j] + 1.45);
            aaa[2] = 0;bbb[2] = 0;ccc[2] = 0;
            g[i][j].AddQuadratic(3, aaa, bbb, ccc, ppp);
        }
        
        // Convex corridor ellipsoid penalty function
        ppp[0] = -inf;ppp[1] = 0.9;ppp[2] = 1.0;ppp[3] = inf;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = lamb5[i]; bbb[1] = 0; ccc[1] = lamb5[i]*(-0.81); 
        aaa[2] = lamb5[i]*2; bbb[2] = 0; ccc[2] = lamb5[i]*(-1.81);
        g[i][SizeEqx].AddQuadratic(3, aaa, bbb, ccc, ppp);

        // Trust Region Constraints
        ppp[0] = -inf;ppp[1] = 0;ppp[2] = 0.3*v_norm[i];ppp[3] = inf;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = lamb3[i]; bbb[1] = 0; ccc[1] = 0;
        aaa[2] = lamb3[i]*10; bbb[2] = 0; ccc[2] = lamb3[i]*(-9)*(0.3*v_norm[i])*(0.3*v_norm[i]);
        g[i][SizeEqx+1].AddQuadratic(3, aaa, bbb, ccc, ppp);

        // Curvature Penalty Function
        ppp[0] = -inf;ppp[1] = 0.5;ppp[2] = 2.0;ppp[3] = inf;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = lamb4[i]; bbb[1] = 0; ccc[1] = lamb4[i]*(-0.25);
        aaa[2] = lamb4[i]*10; bbb[2] = 0; ccc[2] = lamb4[i]*(-36.25);
        g[i][SizeEqx+2].AddQuadratic(3, aaa, bbb, ccc, ppp);
    }

    MPC_ADMMSolver<float, HorizonNum, SizeX2d, SizeU2d, SizeYx2d, SizeYu2d, SizeEqx2d, SizeEqu2d, NumEllx2d, NumEllu2d, SizeG2d, SizeEllx2d, SizeEllu2d> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, new_centerX, new_centerU, g, 100, KAcc, K);
    auto start_time = std::chrono::high_resolution_clock::now();
    res = mpc.solve();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Time taken by ADMM solver: " << duration.count() << " seconds" << std::endl;

    return;
}
*/


int solveMpc2D(vec_E<Polyhedron<2>>& mpc_polyhedrons, std::array<Eigen::Matrix<float, SizeYx2d - SizeEqx, 1>, HorizonNum + 1>& new_centerX, std::array<Eigen::Matrix<float, SizeYu2d - SizeEqu, 1>, HorizonNum + 1>& new_centerU, std::array<Eigen::Matrix<float, 2, 2>, HorizonNum + 1>& elliE, vector<float> dt, vector<Eigen::Vector2f>& ref_points, vector<float>& v_norm, vector<Eigen::Matrix<float, 2, 2>>& Rk, BlockVector<float, HorizonNum + 1, SizeX2d + SizeU2d>& res) {
    float q=0.35, st=0, wei=10, weig=50; int K=250;
    // std::cout << "please enter the inner iteration num" << std::endl;
    // cin >> K;
    struct timeval T1,T2;
    float timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<float> Px, Py, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6;
    vector<vector<float>> center;
    vector<vector<Hyperplane<2>>> CorridorP;
    CorridorP.resize(HorizonNum+1);
    Px.resize(HorizonNum+1);
    Py.resize(HorizonNum+1);

    for(int i = 0;i <= HorizonNum; ++i) {
        int plane_size = mpc_polyhedrons[i].hyperplanes().size();
        CorridorP[i].resize(plane_size);
        for (int j = 0; j < plane_size; j++) {
            CorridorP[i][j] = mpc_polyhedrons[i].hyperplanes()[j];
        }
        // 用ref_points为Px、Py、Pz赋值
        Px[i] = ref_points[i].x();
        Py[i] = ref_points[i].y();

        // 2d
        // if (i == 0) {
        //     lamb1.push_back(1);
        // } else if (i >= HorizonNum-2) {
        //     lamb1.push_back(1);
        // } else {
        //     lamb1.push_back(0.01);
        // }
        
        // lamb2.push_back(0.1);//纵向加速度
        // lamb3.push_back(0.1);//信赖域约束
        // lamb4.push_back(100.0);//曲率过大
        // lamb5.push_back(0.00000001);//凸走廊约束
        // lamb6.push_back(5.0);//曲率平方

        // 3d
        if (i == 0) {
            lamb1.push_back(1);
        } else if (i == HorizonNum) {
            lamb1.push_back(1.0);
        } else if (i > HorizonNum-6) {
            lamb1.push_back(0.04);
        } else {
            lamb1.push_back(0.04);
        }
        
        lamb2.push_back(0.0001);//纵向加速度
        lamb3.push_back(0.001);//信赖域约束
        lamb4.push_back(1000);//曲率过大
        lamb5.push_back(0.1);//凸走廊约束 0.1
        lamb6.push_back(0.001);//曲率平方
    }
    // for (int i = 0; i < 100; i++) {
        // solveunit2D(dt, Px, Py, Pz, v_norm, Rk, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, CorridorP, new_centerX,  new_centerU, elliE, res, K);
    // }

    gettimeofday(&T2,NULL);
    timeuse = (T2.tv_sec - T1.tv_sec) + (float)(T2.tv_usec - T1.tv_usec)/1000000.0;
    std::cout<<"time taken by mpc problem: "<<timeuse<< " seconds" << std::endl;  //输出时间（单位：ｓ）
    // std::cout<<"average time taken by mpc problem: "<<timeuse / 100.0f<< " seconds" << std::endl;  //输出时间（单位：ｓ）
    return 0;
}