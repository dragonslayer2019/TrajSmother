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
const int HorizonNum = 49;



//double shift test
const int SizeX = 3, SizeU = 1;


typedef Eigen::Matrix<double, SizeX, SizeX> MatrixX;
typedef Eigen::Matrix<double, SizeU, SizeU> MatrixU;
typedef Eigen::Matrix<double, SizeX, SizeU> MatrixB;
typedef Eigen::Matrix<double, SizeX, 1> VectorX;
typedef Eigen::Matrix<double, SizeU, 1> VectorU;
typedef Eigen::Matrix<T, SizeYx, SizeX> MatrixHx;
VectorX x_init;
MatrixX Ai, Qi, Hxi;
MatrixB Bi;
MatrixU Ri, Hui;
VectorX ci, Li;
VectorU Wi;
const double inf = 1e5;




bool isPointOnSegment(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& P) {
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AP = P - A;
    
    // 计算叉积
    Eigen::Vector3d crossProduct = AB.cross(AP);
    
    // 判断叉积是否为零向量
    if (!crossProduct.isZero(1e-10)) { // 允许一定的误差范围
        return false;
    }
    
    // 计算点积
    double dotProduct = AP.dot(AB);
    double squaredLengthAB = AB.dot(AB);
    
    // 判断点积是否在规定范围内
    if (0 <= dotProduct && dotProduct <= squaredLengthAB) {
        return true;
    } else {
        return false;
    }
}

void solveunit3D(vector<double> dt, vector<double> Px, vector<double> Py, vector<double> Pz, vector<double> lamb1, vector<double> lamb2, vector<double> lamb3, vector<double> lamb4, vector<double> lamb5, vector<vector<Hyperplane<3>>> CorridorP, vector<Ellipsoid<3>> CorridorE, int K = 250) {
    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    std::array<MatrixHx, HorizonNum + 1> Hx;
    std::array<MatrixU, HorizonNum + 1> Hu;

    std::array<MatrixU, HorizonNum> Rk;
    std::array<double, HorizonNum+1> Vnorm;
    std::vector<vector<double>> DistC;
    std::vector<double> aaa, bbb, ccc;
    std::vector<double> ppp;
    // 状态变量(px, py, pz, vx, vy, vz)^T
    // 控制输入(mu1, mu2, mu3)^T
    // std::array<Eigen::Matrix<double, SizeYx - SizeX, 1>, HorizonNum + 1> new_center;
    VectorX x_init; // 这里是一维的pvaj，要改成3维的
    std::array<std::array<FunctionG<double>, SizeYx + SizeYu>, HorizonNum + 1> g;
    x_init << 0, 0, 0, 0;
    
    for(int i = 0; i <= HorizonNum; ++i) {
        // 此处为m个切平面约束+一个椭球约束
        // 计算优化路点到切平面之间距离的结果中的常数转移到gxi函数中
        Hxi << 1, 0, 0, 0, 0, 0,
               0, 1, 0, 0, 0, 0,
               0, 0, 1, 0, 0, 0,
               0, 0, 0, 1, 0, 0;
               
            
        for (int j = 0; j < M; ++j) {
            // m行到切平面距离   
            Hxi << CorridorP[i][j].n().x, CorridorP[i][j].n().y, CorridorP[i][j].n().z, 0, 0, 0;
            // 计算点到切平面距离中的const向量
            DistC[i][j] = -CorridorP[i][j].n().x*CorridorP[i][j].p().x - CorridorP[i][j].n().y*CorridorP[i][j].p().y - CorridorP[i][j].n().y*CorridorP[i][j].p().y;
        }
        // 1行椭球到圆心约束！！！！！！！！！！！！！！！

        Hxi << 1, 0, 0, 0, 0, 0;
        //    cos(Theta[i]) / Ella[i], sin(Theta[i]) / Ella[i], 0, 0,
        //    -sin(Theta[i]) / Ellb[i], cos(Theta[i]) / Ellb[i], 0, 0;

        // 此处为一个曲率平方约束+一个曲率罚函数   ck = sqrt(mu2^2+mu3^2)
        Hui << 0, 1, 1,
               0, 1, 1;
        Hx[i] = Hxi;
        Hu[i] = Hui;
        // new_center[i] << cos(Theta[i]) / Ella[i] * center[i][0] + sin(Theta[i]) / Ella[i] * center[i][1],
        //                  -sin(Theta[i]) / Ellb[i] * center[i][0] + cos(Theta[i]) / Ellb[i] * center[i][1];
    }
    // 根据采样参考点计算n-1个Rk，Rk第一列单位向量与vk平行，且为正交矩阵

    // 根据Vx， Vy， Vz计算Vnorm[i]

    for(int i = 0; i < HorizonNum; ++i) {
        // 6*6
        Ai << 1, 0, 0, dt[i], 0, 0,
              0, 1, 0, 0,  dt[i], 0,
              0, 0, 1 , 0, 0, dt[i],
              0, 0, 0 , 1, 0, 0,
              0, 0, 0 , 0, 1, 0,
              0, 0, 0 , 0, 0, 1;
        // 6*3
        Bi << 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](0,0), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](0,1), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](0,2),
              0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](1,0), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](1,1), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](1,2),
              0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](2,0), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](2,1), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i](2,2),
              Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](0,0), Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](0,1), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](0,2),
              Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](1,0), Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](1,1), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](1,2),
              Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](2,0), Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](2,1), 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i](2,2);
        // 6*1
        ci << 0, 0, 0, 0, 0, 0;
        A[i] = Ai; B[i] = Bi; c[i] = ci;
    }

    for(int i = 0; i <= HorizonNum; ++i) {
        // 追参考位置
        Qi << lamb1[i], 0, 0, 0, 0, 0,
              0, lamb1[i], 0, 0, 0, 0,
              0, 0, lamb1[i], 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0;
        Li << -Px[i] * lamb1[i],
              -Py[i] * lamb1[i],
              -Pz[i] * lamb1[i],
              0,
              0,
              0;
        // 限制纵向加速度
        Ri << 0, 0, 0,
              0, 0, 0,
              0, 0, lamb2[i];
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
            ppp[0] = -inf;ppp[1] = 0.1;ppp[2] = 1;ppp[3] = inf;
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


    MPC_ADMMSolver<double, HorizonNum, SizeX, SizeU, SizeYx, SizeYu> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, new_center, g, 100, K);
    BlockVector<double, HorizonNum + 1, SizeX + SizeU> res = mpc.solve();
    

    using json = nlohmann::json;
    json j;
    json jx = json::array();
    json jy = json::array();
    json jv = json::array();
    json ja = json::array();
    json ju = json::array();
    // vector<double> t(HorizonNum), x(HorizonNum), v(HorizonNum), a(HorizonNum), u(HorizonNum);
    for(int i = 0; i <= HorizonNum; ++i) {
        jx.push_back(res.v[i](0, 0));
        jy.push_back(res.v[i](1, 0));
        jv.push_back(res.v[i](2, 0));
        ja.push_back(res.v[i](3, 0));
        ju.push_back(res.v[i](4, 0));
    }
    j["x"] = jx;
    j["y"] = jy;
    j["v"] = jv;
    j["a"] = ja;
    j["u"] = ju;
    j["Px"] = Px;
    ofstream out("test1.out");
    if(out.is_open()) {
        out << j.dump(4) << endl;
        out.close();
    }
    return;
}

int solveMpc(vec_E<Polyhedron<3>> mpc_polyhedrons, vec_E<Ellipsoid<3>> mpc_ellipsoids, vector<double> dt, vec_Vec3f ref_points) {
    double q=0.35, st=0, wei=10, weig=50; int K=250;
    cin >> K;
    struct timeval T1,T2;
    double timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<double> Px, Py, Pz, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5;
    vector<vector<double>> center;
    vector<vector<Hyperplane<3>>> CorridorP;//分段切平面约束
    vector<Ellipsoid<3>> CorridorE;//分段椭球
    CorridorP.resize(HorizonNum+1);
    CorridorE.resize(HorizonNum+1);
    Px.resize(HorizonNum+1);
    Py.resize(HorizonNum+1);
    Pz.resize(HorizonNum+1);


    for(int i = 0;i <= HorizonNum; ++i) {
        // 切平面约束与椭球约束
        int plane_size = mpc_polyhedrons[i].hyperplanes().size();
        CorridorP[i].resize(plane_size);
        for (int j = 0; j < plane_size; j++) {
            CorridorP[i][j] = mpc_polyhedrons[i].hyperplanes()[j];
        }
        CorridorE[i] = mpc_ellipsoids[i];

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
    }
    solveunit3D(dt, Px, Py, Pz,lamb1, lamb2, lamb3, lamb4, lamb5, CorridorP, CorridorE, K);
    gettimeofday(&T2,NULL);
    timeuse = (T2.tv_sec - T1.tv_sec) + (double)(T2.tv_usec - T1.tv_usec)/1000000.0;
    cout<<"time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
    return 0;
}




// 计算两个三维点之间的欧氏距离
double distance(const Vec3f& p1, const Vec3f& p2) {
    return sqrt(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2) + pow(p1(2) - p2(2), 2));
}

// 插值路径点，步长为max_length
vec_Vec3f interpolatePoints(const vec_Vec3f& points, double max_length) {
    vec_Vec3f interpolated_points;
    interpolated_points.push_back(points[0]);

    for (size_t i = 1; i < points.size(); ++i) {
        Vec3f p1 = points[i - 1];
        Vec3f p2 = points[i];
        double segment_length = distance(p1, p2);
        size_t num_segments = static_cast<size_t>(ceil(segment_length / max_length));
        Eigen::Vector3d direction(p2(0) - p1(0), p2(1) - p1(1), p2(2) - p1(2));
        direction.normalize();

        for (size_t j = 1; j < num_segments; ++j) {
            Vec3f new_point = {
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
vector<double> computeCurvature(const vec_Vec3f& points, double resample_dist) {
    int n = points.size();
    vector<double> curvature(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        Eigen::Vector3d v1(points[i](0) - points[i - 1](0), points[i](1) - points[i - 1](1), points[i](2) - points[i - 1](2));
        Eigen::Vector3d v2(points[i + 1](0) - points[i](0), points[i + 1](1) - points[i](1), points[i + 1](2) - points[i](2));
        double angle = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
        curvature[i] = angle / resample_dist;
    }
    return curvature;
}

// 高斯平滑
vector<double> smoothCurvature(const vector<double>& curvature, int window_size) {
    int n = curvature.size();
    vector<double> smoothed_curvature(n, 0.0);
    int half_window = window_size / 2;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
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
vec_Vec3f resamplePath(const vec_Vec3f& points, const vector<double>& smoothed_curvature, double max_length, double max_angle) {
    vec_Vec3f resampled_points;
    resampled_points.push_back(points[0]);
    double cum_length = 0.0;
    double cum_angle = 0.0;

    for (size_t i = 1; i < points.size(); ++i) {
        double segment_length = distance(points[i], points[i - 1]);
        cum_length += segment_length;

        if (i > 1 && i <= smoothed_curvature.size()) {
            double delta_angle = smoothed_curvature[i] * segment_length;
            cum_angle += delta_angle;
        }

        if (cum_length > max_length || cum_angle > (max_angle * M_PI / 180.0)) {
            Eigen::Vector3d direction(points[i](0) - points[i - 1](0), points[i](1) - points[i - 1](1), points[i](2) - points[i - 1](2));
            direction.normalize();
            Vec3f new_point = {
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

vec_Vec3f get_sample_point(vec_Vec3f& path) {
    // 输入路径点
    // vec_Vec3f points = {
    //     {5, 9.5, 0.5}, {13, 11.5, 3.0}, {15, 9.5, 1.5}, {14, 5, 2.5}
    // };

    // 步长
    double max_length = 0.5;

    // 插值路径点，确保每段长度不超过max_length
    vec_Vec3f interpolated_points = interpolatePoints(path, max_length);

    // 计算曲率
    vector<double> curvature = computeCurvature(interpolated_points, max_length);

    // 平滑曲率
    vector<double> smoothed_curvature = smoothCurvature(curvature, 8);

    // 重新采样路径点
    double max_angle = 20.0;
    vec_Vec3f resampled_points = resamplePath(interpolated_points, smoothed_curvature, max_length, max_angle);

    // 输出重新采样后的路径点
    cout << "Resampled Path Points:" << endl;
    for (const auto& point : resampled_points) {
        cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << endl;
    }

    return resampled_points;
}