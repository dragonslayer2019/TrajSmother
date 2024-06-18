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
  std::cout << "Time taken by get corridor: " << duration.count() << " microseconds " << std::endl;
  
  //0.计算路径曲率并平滑,基于曲率采样任意段弧长，经过裁剪路点数量得到不多于HorizonNum+1个采样路点
  vec_Vec3f sampled_points = get_sample_point(path);
  constexpr int HorizonNum = 49;
  vec_Vec3f ref_points;
  ref_points.resize(HorizonNum + 1);
  if (sampled_points.size() - 1 < HorizonNum) {
    std::copy(sampled_points.begin(), sampled_points.end(), ref_points.begin());
    std::fill(ref_points.begin() + sampled_points.size(), ref_points.end(), sampled_points.back());
  } else {
    std::copy(sampled_points.begin(), sampled_points.begin() + HorizonNum + 1, ref_points.begin());
  }

  //1.基于弧长长度与最大单位弧长长度插值求解v，根据弧长vector与对应的V求解非均衡dt
  vector<double> v_norm;
  vector<double> dt;
  double v_max = 5;
  double s_max = 0.5;
  for (int i = 0; i < HorizonNum; ++i) {
    double dis = distance(ref_points[i], ref_points[i+1]);
    v_norm.push_back(dis/s_max*v_max);
    dt.push_back(dis/v_norm[i]);
  }
  v_norm.push_back(v_norm.back());
  

  //2.将采样路点与分段凸走廊匹配，horizon num + 1个点对应horizon num + 1个凸走廊，在第i段有m个凸走廊约束,让Hxi正确匹配到对应的凸走廊约束
  vec_E<Polyhedron<3>> polyhedrons = decomp_util.get_polyhedrons();
  vec_E<Ellipsoid<3>> ellipsoids = decomp_util.get_ellipsoids();
  vec_E<Polyhedron<3>> mpc_polyhedrons;
  vec_E<Ellipsoid<3>> mpc_ellipsoids;
  for (int i = 0; i <= HorizonNum; ++i) {
    //判断路点处于哪段，并范回对应的段数id
    for (int j = 0; j < path.size()-1;j++){
      if(isPointOnSegment(path[j], path[j+1], ref_points[i])) {
        mpc_polyhedrons.push_back(polyhedrons[j]);
        mpc_ellipsoids.push_back(ellipsoids[j]);
        break;
      }
    }
  }

  std::cout << "start solve MPC" << std::endl;
  solveMpc(mpc_polyhedrons, mpc_ellipsoids, dt, ref_points, v_norm);
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
