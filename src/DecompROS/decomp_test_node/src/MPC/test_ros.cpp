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


void solveunit(vector<double> dt, vector<double> Xref, vector<double> Vref, vector<double> Amin, vector<double> Amax, vector<double> Xsafe1, vector<double> Xsafe2, vector<double> Vlim, vector<double> Vlaw, vector<double> B1, vector<double> lamb1, vector<double> lamb2, vector<double> lamb3, vector<double> lamb4, vector<double> lamb5, vector<double> lamb6, vector<double> lamb7, vector<double> lamb8, int K = 250) {
    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    std::array<MatrixX, HorizonNum + 1> Hx;
    std::array<MatrixU, HorizonNum + 1> Hu;
    VectorX x_init;
    std::array<std::array<FunctionG<double>, SizeX + SizeU>, HorizonNum + 1> g;
    x_init << 0, 0, 0;
    Hxi << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;
    Hui << 1;
    for(int i = 0; i <= HorizonNum; ++i) {
        Hx[i] = Hxi; Hu[i] = Hui;
    }
    for(int i = 0; i < HorizonNum; ++i) {
        Ai << 1, dt[i], dt[i] * dt[i] * 0.5, 
              0, 1 ,    dt[i],
              0, 0 ,    1;
        Bi << 0,
              dt[i] * dt[i] * 0.5,
              dt[i];
        ci << 0, 0, 0;
        A[i] = Ai; B[i] = Bi; c[i] = ci;
    }
    for(int i = 0; i <= HorizonNum; ++i) {
        Qi << lamb1[i], 0, 0,
              0, lamb2[i], 0,
              0, 0, lamb3[i];
        Ri << lamb4[i];
        Li << -Xref[i] * lamb1[i],
              -Vref[i] * lamb2[i],
              0;
        Wi << 0;
        Q[i] = Qi; R[i] = Ri; L[i] = Li; W[i] = Wi;
    }
    for(int i = 0;i <= 3; ++i) {
        aaa.push_back(0);
        bbb.push_back(0);
        ccc.push_back(0);
        ppp.push_back(0);
    }
    for(int i = 0; i <= HorizonNum; ++i) {
        g[i][1].AddIndicator(-inf, inf);
        g[i][2].AddIndicator(Amin[i], Amax[i]);
        
        double lamb = lamb5[i] / 100 * 0;
        double Xsafem = (Xsafe1[i] + Xsafe2[i]) / 2;
        ppp[0] = -inf;
        aaa[0] = (lamb5[i] + lamb); bbb[0] = (lamb5[i] + lamb) * (-2) * Xsafe1[i]; ccc[0] = (lamb5[i] + lamb) * Xsafe1[i] * Xsafe1[i] + lamb * (Xsafem - Xsafe1[i]) * (Xsafem - Xsafe1[i]);
        ppp[1] = Xsafe1[i];
        aaa[1] = lamb; bbb[1] = lamb * (-2) * Xsafem; ccc[1] = lamb * Xsafem * Xsafem; 
        ppp[2] = Xsafe2[i];
        aaa[2] = (lamb5[i] + lamb); bbb[2] = (lamb5[i] + lamb) * (-2) * Xsafe2[i]; ccc[2] = (lamb5[i] + lamb) * Xsafe2[i] * Xsafe2[i] + lamb * (Xsafem - Xsafe2[i]) * (Xsafem - Xsafe2[i]);
        ppp[3] = inf;
        g[i][0].AddQuadratic(3, aaa, bbb, ccc, ppp);

        double lim1, lim2, weig1, weig2;
        if(Vlim[i] < Vlaw[i]) {
            lim1 = Vlim[i]; lim2 = Vlaw[i];
            weig1 = lamb6[i]; weig2 = lamb8[i];
        } else {
            lim2 = Vlim[i]; lim1 = Vlaw[i];
            weig2 = lamb6[i]; weig1 = lamb8[i];
        }
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0; ppp[0] = -inf;
        aaa[1] = weig1; bbb[1] = weig1 * (-2) * lim1; ccc[1] = weig1 * lim1 * lim1; ppp[1] = lim1;
        aaa[2] = weig1 + weig2;
        bbb[2] = weig1 * (-2) * lim1 + weig2 * (-2) * lim2;
        ccc[2] = weig1 * lim1 * lim1 + weig2 * lim2 * lim2;
        ppp[2] = lim2;
        ppp[3] = inf;
        g[i][1].AddQuadratic(3, aaa, bbb, ccc, ppp);

        ppp[0] = -inf;
        aaa[0] = lamb7[i]; bbb[0] = lamb7[i] * (-2) * B1[i]; ccc[0] = lamb7[i] * B1[i] * B1[i];
        aaa[1] = 0; bbb[1] = 0; ccc[1] = 0; ppp[1] = B1[i]; ppp[2] = inf;
        g[i][2].AddQuadratic(2, aaa, bbb, ccc, ppp);
    }


    MPC_ADMMSolver<double, HorizonNum, SizeX, SizeU, SizeX, SizeU> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, g, 100, K);
    BlockVector<double, HorizonNum + 1, SizeX + SizeU> res = mpc.solve();
    

    using json = nlohmann::json;
    json j;
    json jx = json::array();
    json jv = json::array();
    json ja = json::array();
    json ju = json::array();
    // vector<double> t(HorizonNum), x(HorizonNum), v(HorizonNum), a(HorizonNum), u(HorizonNum);
    for(int i = 0; i <= HorizonNum; ++i) {
        jx.push_back(res.v[i](0, 0));
        jv.push_back(res.v[i](1, 0));
        ja.push_back(res.v[i](2, 0));
        ju.push_back(res.v[i](3, 0));
    }
    j["x"] = jx;
    j["v"] = jv;
    j["a"] = ja;
    j["u"] = ju;
    j["Xref"] = Xref;
    ofstream out("test.out");
    if(out.is_open()) {
        out << j.dump(4) << endl;
        out.close();
    }
    //res.print("RESULT");
    return ;
}


int solveMpc() {
    double q=0.35, st=0, wei=10, weig=50; int K=250;
    cin >> K;
    struct timeval T1,T2;
    double timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<double> dt, Xref, Vref, Amax, Amin, Vcon, Vlaw, Xsafe1, Xsafe2, B1, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, lamb7, lamb8;
    /*for (int i = 0; i < 15; ++i) {
        Xref.push_back(i * 0.2 * 0.33);
    }
    for (int i = 15; i < 35; ++i) {
        Xref.push_back(1 - (i - 15) * 0.2 * 0.5);
    }
    for (int i = 35; i <= HorizonNum; ++i) {
        Xref.push_back(-1 + (i - 35) * 0.2 * 0.33);
    }*/

    for(int i = 0;i <= HorizonNum; ++i) {
        dt.push_back(0.2);
        Xref.push_back(0);
        Vref.push_back(0);
        Amax.push_back(10);
        Amin.push_back(-10);
        Vcon.push_back(25);
        Vlaw.push_back(30);
        B1.push_back(-inf);
        lamb1.push_back(1);
        lamb2.push_back(1);
        // lamb3.push_back(3);
        // lamb4.push_back(2);
        // lamb5.push_back(3);
        // lamb6.push_back(100);
        // lamb7.push_back(0.6);
        // lamb8.push_back(3000);
        lamb3.push_back(1);
        lamb4.push_back(2);
        lamb5.push_back(3000);
        lamb6.push_back(1);
        lamb7.push_back(0.6);
        lamb8.push_back(3000);
        if (i < 20 && i >= 11 ){
            Xsafe1.push_back(0.7);
            Xsafe2.push_back(1.2);
        } else if (i < 40 && i >= 31) {
            Xsafe1.push_back(-1.2);
            Xsafe2.push_back(-0.5);
        } else {
            Xsafe1.push_back(-inf);
            Xsafe2.push_back(inf);
        }
    }
    std::cout << "begin solve MPC unit" << std::endl;
    solveunit(dt, Xref, Vref, Amin, Amax, Xsafe1, Xsafe2, Vcon, Vlaw, B1, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, lamb7, lamb8, K);
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
