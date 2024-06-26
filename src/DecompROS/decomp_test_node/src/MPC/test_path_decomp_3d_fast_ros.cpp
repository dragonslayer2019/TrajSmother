// #include "/home/alan/TrajSmother/src/DecompROS/decomp_test_node/src/bag_reader.hpp"
// #include "/home/alan/TrajSmother/src/DecompROS/decomp_test_node/src/txt_reader.hpp"
#include "../bag_reader.hpp"
#include "../txt_reader.hpp"
#include <decomp_ros_utils/data_ros_utils.h>
#include <ros/ros.h>
#include <decomp_util/ellipsoid_decomp.h>
#include <sensor_msgs/point_cloud_conversion.h>
#include <nav_msgs/Path.h>
#include <chrono>
#include <visualization_msgs/Marker.h>
#include "presteps.h"

int main(int argc, char ** argv){
  ros::init(argc, argv, "test");
  ros::NodeHandle nh("~");

  ros::Publisher cloud_pub = nh.advertise<sensor_msgs::PointCloud>("cloud", 1, true);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("path", 1, true);
  ros::Publisher ref_path_pub = nh.advertise<nav_msgs::Path>("ref_path_path", 1, true);
  ros::Publisher smoothed_path_pub = nh.advertise<nav_msgs::Path>("smooth_path", 1, true);
  ros::Publisher es_pub = nh.advertise<decomp_ros_msgs::EllipsoidArray>("ellipsoid_array", 1, true);
  ros::Publisher poly_pub = nh.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedron_array", 1, true);
  ros::Publisher marker_pub = nh.advertise<visualization_msgs::Marker>("visualization_marker", 1);
  ros::Publisher marker_pub1 = nh.advertise<visualization_msgs::Marker>("visualization_marker1", 1);

  visualization_msgs::Marker v_points, line_strip, v_points1, line_strip1;
  v_points.header.frame_id = line_strip.header.frame_id = v_points1.header.frame_id = line_strip1.header.frame_id = "map";  // 坐标系
  v_points.header.stamp = line_strip.header.stamp = v_points1.header.stamp = line_strip1.header.stamp = ros::Time::now();
  v_points.ns = line_strip.ns = v_points1.ns = line_strip1.ns = "points_and_lines";
  v_points.action = line_strip.action = v_points1.action = line_strip1.action = visualization_msgs::Marker::ADD;
  v_points.pose.orientation.w = line_strip.pose.orientation.w = v_points1.pose.orientation.w = line_strip1.pose.orientation.w = 1.0;
    // 设置类型
  v_points.id = 0;
  line_strip.id = 1;
  v_points1.id = 3;
  line_strip1.id = 4;
  v_points.type = v_points1.type = visualization_msgs::Marker::POINTS;
  line_strip.type = line_strip1.type = visualization_msgs::Marker::LINE_STRIP;
  // 设置尺寸
  v_points.scale.x = v_points1.scale.x = 0.2;
  v_points.scale.y = v_points1.scale.y = 0.2;
  line_strip.scale.x = line_strip1.scale.x = 0.1;
  // 设置颜色
  v_points.color.g = 1.0f;
  v_points1.color.b = 1.0f;
  v_points1.color.g = 0.5f;
  v_points.color.a = v_points1.color.a = 1.0;
  line_strip.color.b = 1.0;
  line_strip1.color.r = 1.0;
  line_strip.color.a = line_strip1.color.a = 1.0;

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
  std::chrono::duration<double> duration = end_time - start_time;
  // auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count()*1000;
  std::cout << "Time taken by corridor: " << duration.count() << " seconds" << std::endl;

  // 计算路径曲率并平滑,基于曲率采样任意段弧长，经过标准化得到1 + HorizonNum个采样路点
  std::vector<Eigen::Vector3f> path_f;
  path_f.reserve(path.size());
  for (const auto& point : path) {
      path_f.emplace_back(point.cast<float>());
  }
  vector<Eigen::Vector3f> ref_points = get_sample_point(path_f);

  // 参考路径可视化
  vec_Vec3f ref_path(HorizonNum+1);
  Eigen::Vector3f temp1;
  for (int i = 0; i <= HorizonNum; i++) {
    ref_path[i] = ref_points[i].cast<double>();
  }
  nav_msgs::Path ref_path_msg = DecompROS::vec_to_path(ref_path);
  ref_path_msg.header.frame_id = "map";
  ref_path_pub.publish(ref_path_msg);

  // 基于弧长长度与最大单位弧长长度插值求解v，根据弧长vector与对应的V求解非均匀dt，根据v的方向计算rk
  vector<float> v_norm, dt;
  // vector<Mat3f> Rk;
  vector<Eigen::Matrix<float, 3, 3>> Rk;
  get_param(ref_points, v_norm, Rk, dt);

  // 将采样路点与分段凸走廊匹配，horizon num + 1个点对应horizon num + 1个凸走廊，在第i段有m个凸走廊约束,让Hxi正确匹配到对应的凸走廊约束
  vec_E<Polyhedron<3>> polyhedrons = decomp_util.get_polyhedrons();
  vec_E<Ellipsoid<3>> ellipsoids = decomp_util.get_limit_ellipsoids();
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
  for (auto& mat : new_centerU) {
    mat.setZero();
  }
  // std::array<Mat3f, HorizonNum + 1> elliE;
  std::array<Eigen::Matrix<float, 3, 3>, HorizonNum + 1> elliE;
  
  for (int i = 0; i <= HorizonNum; ++i) {
    //判断路点处于哪段，并范回对应的段数id
    for (int j = 0; j < path_f.size()-1;j++){
      if(isPointOnSegment(path_f[j], path_f[j+1], ref_points[i])) {
        mpc_polyhedrons.push_back(polyhedrons[j]);
        elliE[i] = E_inverse[j];
        new_centerX[i].block(0, 0, 3, 1) = E_inverse_d[j];
        Eigen::Vector3f vk = (ref_path[i+1] - ref_path[i]).cast<float>().normalized()*v_norm[i];
        Eigen::Matrix<float, 3, 1> v_center;
        v_center << vk[0], vk[1], vk[2]; 
        new_centerX[i].block(3, 0, 3, 1) = v_center;
        break;
      }
    }
  }

  std::cout << "start solve MPC" << std::endl;
  BlockVector<float, HorizonNum + 1, SizeX + SizeU> res;
  res.setZero();
  solveMpc(mpc_polyhedrons, new_centerX, new_centerU, elliE, dt, ref_points, v_norm, Rk, res);
  std::cout << "end solveMpc" << std::endl;


  //Publish visualization msgs
  // decomp_ros_msgs::EllipsoidArray es_msg = DecompROS::ellipsoid_array_to_ros(decomp_util.get_ellipsoids());
  decomp_ros_msgs::EllipsoidArray es_msg = DecompROS::ellipsoid_array_to_ros(decomp_util.get_limit_ellipsoids());
  es_msg.header.frame_id = "map";
  es_pub.publish(es_msg);

  decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(decomp_util.get_polyhedrons());
  poly_msg.header.frame_id = "map";
  poly_pub.publish(poly_msg);


  // 优化后的路径可视化
  vec_Vec3f smooth_path(HorizonNum+1);
  Eigen::Vector3f temp;
  for (int i = 0; i <= HorizonNum; i++) {
    temp << res.v[i](0, 0), res.v[i](1, 0), res.v[i](2, 0);
    smooth_path[i] = temp.cast<double>();
  }
  nav_msgs::Path smooth_path_msg = DecompROS::vec_to_path(smooth_path);
  smooth_path_msg.header.frame_id = "map";
  smoothed_path_pub.publish(smooth_path_msg);

  std::vector<geometry_msgs::Point> points_list;
  std::vector<geometry_msgs::Point> points1_list;
  for (int i = 0; i <= HorizonNum; i++) {
      geometry_msgs::Point p;
      p.x = res.v[i](0, 0);
      p.y = res.v[i](1, 0);
      p.z = res.v[i](2, 0);
      points_list.push_back(p);
      geometry_msgs::Point p1;
      p1.x = ref_path[i].x();
      p1.y = ref_path[i].y();
      p1.z = ref_path[i].z();
      points1_list.push_back(p1);
  }
  for (const auto& p : points_list) {
      v_points.points.push_back(p);
      line_strip.points.push_back(p);
  }
  for (const auto& p1 : points1_list) {
      v_points1.points.push_back(p1);
      line_strip1.points.push_back(p1);
  }

  std::vector<visualization_msgs::Marker> text_markers;
  for (size_t i = 0; i < points_list.size(); ++i) {
      visualization_msgs::Marker text_marker;
      text_marker.header.frame_id = "map";
      text_marker.header.stamp = ros::Time::now();
      text_marker.ns = "text_markers";
      text_marker.id = i + 4; // 从4开始，避免与点和线的ID冲突
      text_marker.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
      text_marker.action = visualization_msgs::Marker::ADD;
      text_marker.pose.position = points_list[i];
      text_marker.pose.position.z += 0.2; // 使文本稍微高于点
      text_marker.pose.orientation.w = 1.0;
      text_marker.scale.z = 0.2; // 文本的高度
      text_marker.color.r = 0.0f;
      text_marker.color.g = 1.0f;
      text_marker.color.b = 0.0f;
      text_marker.color.a = 1.0;
      text_marker.text = std::to_string(i); // 设置文本为点的ID
      text_markers.push_back(text_marker);
  }

  std::vector<visualization_msgs::Marker> text_markers1;
  for (size_t i = 0; i < points1_list.size(); ++i) {
      visualization_msgs::Marker text_marker1;
      text_marker1.header.frame_id = "map";
      text_marker1.header.stamp = ros::Time::now();
      text_marker1.ns = "text_marker1";
      text_marker1.id = i + 100; // 从100开始，避免与点和线的ID冲突
      text_marker1.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
      text_marker1.action = visualization_msgs::Marker::ADD;
      text_marker1.pose.position = points1_list[i];
      text_marker1.pose.position.z += 0.2; // 使文本稍微高于点
      text_marker1.pose.orientation.w = 1.0;
      text_marker1.scale.z = 0.2; // 文本的高度
      text_marker1.color.r = 0.0f;
      text_marker1.color.g = 0.5f;
      text_marker1.color.b = 1.0f;
      text_marker1.color.a = 1.0;
      text_marker1.text = std::to_string(i); // 设置文本为点的ID
      text_markers1.push_back(text_marker1);
  }

  // //Convert to inequality constraints Ax < b
  // auto polys = decomp_util.get_polyhedrons();
  // for(size_t i = 0; i < path.size() - 1; i++) {
  //   const auto pt_inside = (path[i] + path[i+1]) / 2;
  //   LinearConstraint3D cs(pt_inside, polys[i].hyperplanes());
  //   printf("i: %zu\n", i);
  //   std::cout << "A: " << cs.A() << std::endl;
  //   std::cout << "b: " << cs.b() << std::endl;
  //   std::cout << "point: " << path[i].transpose();
  //   if(cs.inside(path[i]))
  //     std::cout << " is inside!" << std::endl;
  //   else
  //     std::cout << " is outside!" << std::endl;

  //   std::cout << "point: " << path[i+1].transpose();
  //   if(cs.inside(path[i+1]))
  //     std::cout << " is inside!" << std::endl;
  //   else
  //     std::cout << " is outside!" << std::endl;
  // }

  ros::Rate r(30);
  while (ros::ok()) {
      marker_pub.publish(v_points);
      marker_pub.publish(line_strip);
      for (const auto& text_marker : text_markers) {
          marker_pub.publish(text_marker);
      }
      marker_pub1.publish(v_points1);
      marker_pub1.publish(line_strip1);
      for (const auto& text_marker1 : text_markers1) {
          marker_pub1.publish(text_marker1);
      }
      r.sleep();
  }
  ros::spin();

  return 0;
}
