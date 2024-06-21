#include "/home/alan/catkin_ws/src/DecompROS/decomp_test_node/src/bag_reader.hpp"
#include "/home/alan/catkin_ws/src/DecompROS/decomp_test_node/src/txt_reader.hpp"
#include <decomp_ros_utils/data_ros_utils.h>
#include <ros/ros.h>
#include <decomp_util/ellipsoid_decomp.h>
#include <sensor_msgs/point_cloud_conversion.h>
#include <nav_msgs/Path.h>
#include <chrono>
#include "presteps.h"

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

  // check path length
  float min_length = HorizonNum * 0.5;
  float length = 0;
  for (int i = 0; i < path.size()-1; i++) {
    length += distance(path[i].cast<float>(), path[i+1].cast<float>());
  }
  if (length < min_length) {
    // 补一个点令path足够长
    int back_id = path.size() - 1; 
    float delta = 1.0f;
    Eigen::Vector3f direction = (path[back_id].cast<float>() - path[back_id - 1].cast<float>()).normalized();
    Eigen::Vector3d added_point = path[back_id] + (direction * (min_length - length + delta)).cast<double>();
    // std::cout << "back id: " << back_id << "dir: " << direction << std::endl;
    // 延长point
    path.push_back(added_point);
  }

  nav_msgs::Path path_msg = DecompROS::vec_to_path(path);
  path_msg.header.frame_id = "map";
  path_pub.publish(path_msg);

  auto start_time = std::chrono::high_resolution_clock::now();
  //Using ellipsoid decomposition
  EllipsoidDecomp3D decomp_util;
  decomp_util.set_obs(obs);
  decomp_util.set_local_bbox(Vec3f(1, 2, 1));
  decomp_util.dilate_fast(path); //Set max iteration number of 10, do fix the path
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  // auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count()*1000;
  std::cout << "Time taken by code: " << duration.count() << " microseconds" << std::endl;
  // std::cout << "tusen Time taken by code: " << duration << std::endl;

  // 计算路径曲率并平滑,基于曲率采样任意段弧长，经过标准化得到1 + HorizonNum个采样路点
  std::vector<Eigen::Vector3f> path_f;
  path_f.reserve(path.size());
  for (const auto& point : path) {
      path_f.emplace_back(point.cast<float>());
  }
  vector<Eigen::Vector3f> ref_points = get_sample_point(path_f);

  // 基于弧长长度与最大单位弧长长度插值求解v，根据弧长vector与对应的V求解非均匀dt，根据v的方向计算rk
  vector<float> v_norm, dt;
  // vector<Mat3f> Rk;
  vector<Eigen::Matrix<float, 3, 3>> Rk;
  get_param(ref_points, v_norm, Rk, dt);

  // 将采样路点与分段凸走廊匹配，horizon num + 1个点对应horizon num + 1个凸走廊，在第i段有m个凸走廊约束,让Hxi正确匹配到对应的凸走廊约束
  vec_E<Polyhedron<3>> polyhedrons = decomp_util.get_polyhedrons();
  vec_E<Ellipsoid<3>> ellipsoids = decomp_util.get_ellipsoids();
  vec_E<Polyhedron<3>> mpc_polyhedrons;
  
  // vector<Mat3f> E_inverse(ellipsoids.size());
  vector<Eigen::Matrix<float, 3, 3>> E_inverse(ellipsoids.size());
  vector<Eigen::Matrix<float, 3, 1>> E_inverse_d(ellipsoids.size());
  
  // 椭球线性变换为单位球
  for (int i = 0; i < ellipsoids.size(); i++) {
    E_inverse[i] = ellipsoids[i].C().inverse().cast<float>();
    E_inverse_d[i] = E_inverse[i]*(ellipsoids[i].d().cast<float>());
  }

  std::array<Eigen::Matrix<float, SizeYx - SizeEqx, 1>, HorizonNum + 1> new_centerX;
  std::array<Eigen::Matrix<float, SizeYu - SizeEqu, 1>, HorizonNum + 1> new_centerU;
  // std::array<Mat3f, HorizonNum + 1> elliE;
  std::array<Eigen::Matrix<float, 3, 3>, HorizonNum + 1> elliE;
  
  for (int i = 0; i <= HorizonNum; ++i) {
    std::cout << "step: " << i << std::endl;
    //判断路点处于哪段，并范回对应的段数id
    for (int j = 0; j < path_f.size()-1;j++){
      if(isPointOnSegment(path_f[j], path_f[j+1], ref_points[i])) {
        mpc_polyhedrons.push_back(polyhedrons[j]);
        elliE[i] = E_inverse[j];
        new_centerX[i] = E_inverse_d[j];
        break;
      }
    }
  }

  std::cout << "start solve MPC" << std::endl;
  solveMpc(mpc_polyhedrons, new_centerX, new_centerU, elliE, dt, ref_points, v_norm, Rk);
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
