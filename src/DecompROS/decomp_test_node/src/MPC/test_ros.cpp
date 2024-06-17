#include "/home/alan/catkin_ws/src/DecompROS/decomp_test_node/src/bag_reader.hpp"
#include "/home/alan/catkin_ws/src/DecompROS/decomp_test_node/src/txt_reader.hpp"
#include <decomp_ros_utils/data_ros_utils.h>
#include <ros/ros.h>
#include <decomp_util/ellipsoid_decomp.h>
#include <sensor_msgs/point_cloud_conversion.h>
#include <nav_msgs/Path.h>
#include <chrono>

#include <bits/stdc++.h>
#include <sys/time.h>
#include "BlockMatrix.h"
#include "MPC.h"
#include "FunctionG.h"
#include "json.hpp"

using namespace std;
//double shift test
const int SizeX = 3, SizeU = 1;
const int HorizonNum = 49;

typedef Eigen::Matrix<double, SizeX, SizeX> MatrixX;
typedef Eigen::Matrix<double, SizeU, SizeU> MatrixU;
typedef Eigen::Matrix<double, SizeX, SizeU> MatrixB;
typedef Eigen::Matrix<double, SizeX, 1> VectorX;
typedef Eigen::Matrix<double, SizeU, 1> VectorU;
VectorX x_init;
MatrixX Ai, Qi, Hxi;
MatrixB Bi;
MatrixU Ri, Hui;
VectorX ci, Li;
VectorU Wi;
const double inf = 1e5;

std::vector<double> aaa, bbb, ccc;
std::vector<double> ppp;


void solveunit3D(vector<double> dt, vector<double> Px, vector<double> Py, vector<double> Pz, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> lamb1, vector<double> lamb2, vector<double> lamb3, vector<double> lamb4, vector<double> lamb5, vector<vector<Hyperplane>> CorridorP, vector<vector<Ellipsoid>> CorridorE, int K = 250) {
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
    // 状态变量(px, py, pz, vx, vy, vz)^T
    // 控制输入(mu1, mu2, mu3)^T
    // std::array<Eigen::Matrix<double, SizeYx - SizeX, 1>, HorizonNum + 1> new_center;
    VectorX x_init;
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
    //res.print("RESULT");
    return;
}


int solveMpc() {
    double q=0.35, st=0, wei=10, weig=50; int K=250;
    cin >> K;
    struct timeval T1,T2;
    double timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<double> dt, Px, Vx, Amax, Amin, Vcon, Vlaw, Xsafe1, Xsafe2, B1, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, lamb7, lamb8, lamb9, lamb10;
    vector<vector<double>> center;
    vector<vector<Hyperplane>> CorridorP;//分段切平面约束
    vector<vector<Ellipsoid>> CorridorE;//分段

    
    
    
    



    

    for(int i = 0;i <= HorizonNum; ++i) {
        // dt.push_back(0.2);
        // Px.push_back(0);
        // Vx.push_back(0);
        // Amax.push_back(10);
        // Amin.push_back(-10);

        // 4.将采样路点的位置作为Px、Py、Pz

        // 5.将采样路点的速度作为Vx、Vy、Vz

        lamb1.push_back(1);
        lamb2.push_back(1);
        lamb3.push_back(1);
        lamb4.push_back(1);
        lamb5.push_back(1);
    }
    solveunit3D(center, dt, Px, Vx, Amin, Amax, Xsafe1, Xsafe2, Vcon, Vlaw, B1, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5, CorridorP, CorridorE, K);
    gettimeofday(&T2,NULL);
    timeuse = (T2.tv_sec - T1.tv_sec) + (double)(T2.tv_usec - T1.tv_usec)/1000000.0;
    cout<<"time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
    return 0;
}

int main(int argc, char ** argv){
  ros::init(argc, argv, "test");
  ros::NodeHandle nh("~");

  ros::Publisher cloud_pub = nh.advertise<sensor_msgs::PointCloud>("cloud", 1, true);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("path", 1, true);
  ros::Publisher es_pub = nh.advertise<decomp_ros_msgs::EllipsoidArray>("ellipsoid_array", 1, true);
  ros::Publisher poly_pub = nh.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedron_array", 1, true);

  std::string file_name, topic_name, path_file;

  nh.param("path_file", path_file, std::string("path.txt"));
  nh.param("bag_file", file_name, std::string("voxel_map"));
  nh.param("bag_topic", topic_name, std::string("voxel_map"));
  //Read the point cloud from bag
  sensor_msgs::PointCloud cloud = read_bag<sensor_msgs::PointCloud>(file_name, topic_name);
  cloud.header.frame_id = "map";
  cloud_pub.publish(cloud);
  std::cout << " cloud size: " << cloud.points.size() << std::endl;
  vec_Vec3f obs = DecompROS::cloud_to_vec(cloud);
  std::cout << " obs size: " << obs.size() << std::endl;
  //Read path from txt
  vec_Vec3f path;
  if(!read_path<3>(path_file, path))
    ROS_ERROR("Fail to read a path!");

  nav_msgs::Path path_msg = DecompROS::vec_to_path(path);
  path_msg.header.frame_id = "map";
  path_pub.publish(path_msg);

  auto start_time = std::chrono::high_resolution_clock::now();
  //Using ellipsoid decomposition
  EllipsoidDecomp3D decomp_util;
  decomp_util.set_obs(obs);
  decomp_util.set_local_bbox(Vec3f(1, 2, 1));
  decomp_util.dilate(path); //Set max iteration number of 10, do fix the path
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  std::cout << "Time taken by code1: " << duration.count() << " microseconds " << std::endl;
  
  //1.基于曲率采样HorizonNum段弧长，得到HorizonNum+1个采样路点


  //2.根据弧长vector与对应的V求解非均衡dt


  //3.获取凸走廊结果，一个二维vector，分段的表示全部切平面约束
  vec_E<Polyhedron<3>> polyhedrons = decomp_util.get_polyhedrons();
  for (const auto& polyhedron : polyhedrons) {
    vec_E<Hyperplane<3>> hyperplanes = polyhedron.hyperplanes();
  }
    

  //4.将采样路点与分段凸走廊匹配，horizon num + 1个点对应horizon num + 1个凸走廊，在第i段有m个凸走廊约束,让Hxi正确匹配到对应的凸走廊约束
  vec_E<Polyhedron<3>> mpc_polyhedrons;
  Polyhedron<3> correspond_polyhedron;
  for (/*HorizonNum+1个采样路点*/) {
    //判断路点处于那段，并范回对应的段数id
     //对各段segment上的两个端点与采样点的向量计算点积，结果为0则返回i

    mpc_polyhedrons.push_back(/*根据id取出polyhedrons对应的polyhedron*/);
  }

  std::cout << "start solve MPC" << std::endl;
  solveMpc();
  std::cout << "end solveMpc" << std::endl;

  //Publish visualization msgs
  decomp_ros_msgs::EllipsoidArray es_msg = DecompROS::ellipsoid_array_to_ros(decomp_util.get_ellipsoids());
  es_msg.header.frame_id = "map";
  es_pub.publish(es_msg);

  decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(decomp_util.get_polyhedrons());
  poly_msg.header.frame_id = "map";
  poly_pub.publish(poly_msg);

    //Convert to inequality constraints Ax < b
  auto polys = decomp_util.get_polyhedrons();
  for(size_t i = 0; i < path.size() - 1; i++) {
    const auto pt_inside = (path[i] + path[i+1]) / 2;
    LinearConstraint3D cs(pt_inside, polys[i].hyperplanes());
    printf("i: %zu\n", i);
    std::cout << "A: " << cs.A() << std::endl;
    std::cout << "b: " << cs.b() << std::endl;
    std::cout << "point: " << path[i].transpose();
    if(cs.inside(path[i]))
      std::cout << " is inside!" << std::endl;
    else
      std::cout << " is outside!" << std::endl;

    std::cout << "point: " << path[i+1].transpose();
    if(cs.inside(path[i+1]))
      std::cout << " is inside!" << std::endl;
    else
      std::cout << " is outside!" << std::endl;
  }


  ros::spin();

  return 0;
}
