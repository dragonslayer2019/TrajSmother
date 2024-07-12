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
#include <geometry_msgs/PointStamped.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/JointState.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_ros/transform_broadcaster.h>
// #include "controller.h"
#include "BVP.h"

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
  ros::Publisher marker_pub2 = nh.advertise<visualization_msgs::Marker>("visualization_marker2", 1);
  ros::Publisher marker_pub3 = nh.advertise<visualization_msgs::Marker>("visualization_marker3", 1);

  visualization_msgs::Marker v_points, line_strip, v_points1, line_strip1, v_points2, line_strip2, v_points3, line_strip3;
  v_points.header.frame_id = line_strip.header.frame_id = v_points1.header.frame_id = line_strip1.header.frame_id = v_points2.header.frame_id = line_strip2.header.frame_id = v_points3.header.frame_id = line_strip3.header.frame_id = "map";  // 坐标系
  v_points.header.stamp = line_strip.header.stamp = v_points1.header.stamp = line_strip1.header.stamp = v_points2.header.stamp = line_strip2.header.stamp = v_points3.header.stamp = line_strip3.header.stamp = ros::Time::now();
  v_points.ns = line_strip.ns = v_points1.ns = line_strip1.ns = v_points2.ns = line_strip2.ns = v_points3.ns = line_strip3.ns = "points_and_lines";
  v_points.action = line_strip.action = v_points1.action = line_strip1.action = v_points2.action = line_strip2.action = v_points3.action = line_strip3.action = visualization_msgs::Marker::ADD;
  v_points.pose.orientation.w = line_strip.pose.orientation.w = v_points1.pose.orientation.w = line_strip1.pose.orientation.w = v_points2.pose.orientation.w = line_strip2.pose.orientation.w = v_points3.pose.orientation.w = line_strip3.pose.orientation.w = 1.0;

  v_points.id = 0;
  line_strip.id = 1;
  v_points1.id = 3;
  line_strip1.id = 4;
  v_points2.id = 3000;
  line_strip2.id = 4000;
  v_points3.id = 3001;
  line_strip3.id = 4001;
  v_points.type = v_points1.type = v_points2.type = v_points3.type = visualization_msgs::Marker::POINTS;
  line_strip.type = line_strip1.type = line_strip2.type = line_strip3.type = visualization_msgs::Marker::LINE_STRIP;
  v_points.scale.x = v_points1.scale.x = v_points2.scale.x = v_points3.scale.x = 0.05;
  v_points.scale.y = v_points1.scale.y = v_points2.scale.y = v_points3.scale.y = 0.05;
  line_strip.scale.x = line_strip1.scale.x = line_strip2.scale.x = line_strip3.scale.x = 0.05;
  v_points.color.g = 1.0f;

  v_points1.color.b = 1.0f;
  v_points1.color.g = 0.5f;

  v_points2.color.b = 1.0f;
  v_points2.color.r = 0.5f;

  v_points3.color.g = 1.0f;

  v_points.color.a = v_points1.color.a = v_points2.color.a = v_points3.color.a = 1.0;

  line_strip.color.b = 1.0;
  line_strip1.color.r = 1.0;
  line_strip2.color.b = 1.0;
  line_strip3.color.r = 1.0;
  line_strip.color.a = line_strip1.color.a = line_strip2.color.a = line_strip3.color.a = 1.0;

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
  float min_length = HorizonNum * max_length;
  float length = 0;
  for (int i = 0; i < path.size()-1; i++) {
    length += distance(path[i].cast<float>(), path[i+1].cast<float>());
  }
  if (length < min_length) {
    // Add a point to make the path long enough
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
  std::cout << "Time taken by corridor: " << duration.count() << " seconds" << std::endl;

  // Calculate the path curvature and smooth it, sample any arc length based on the curvature,
  // and get 1 + HorizonNum sampling waypoints after standardization
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

  // Solve v based on arc length and maximum unit arc length interpolation,
  // solve non-uniform dt based on arc length vector and corresponding V,
  // and calculate rk based on the direction of v
  vector<float> v_norm, dt;
  vector<Eigen::Matrix<float, 3, 3>> Rk;
  get_param(ref_points, v_norm, Rk, dt);

  vec_E<Polyhedron<3>> polyhedrons = decomp_util.get_polyhedrons();
  vec_E<Ellipsoid<3>> ellipsoids = decomp_util.get_limit_ellipsoids();
  vec_E<Polyhedron<3>> mpc_polyhedrons;
  
  // vector<Mat3f> E_inverse(ellipsoids.size());
  vector<Eigen::Matrix<float, 3, 3>> E_inverse(ellipsoids.size());
  vector<Eigen::Matrix<float, 3, 1>> E_inverse_d(ellipsoids.size());
  
  // Linear transformation of ellipsoid into unit sphere
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
    for (int j = 0; j < path_f.size()-1;j++){
      if(isPointOnSegment(path_f[j], path_f[j+1], ref_points[i])) {
        // std::cout << "point " << i << " is on " << j << " lane" << std::endl;
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
  
  vector<Eigen::Vector3f> res_points;
  res_points.resize(HorizonNum + 1);
  for (int i = 0; i < res_points.size(); i++) {
    res_points[i].x() = res.v[i](0, 0);
    res_points[i].y() = res.v[i](1, 0);
    res_points[i].z() = res.v[i](2, 0);
  }

  vector<float> path_points;
  float temp_dis = 0.0f;
  path_points.push_back(temp_dis);
  for (int i = 0; i < res_points.size() - 1; i++) {
      temp_dis += distance(res_points[i], res_points[i+1]);
      path_points.push_back(temp_dis);
  }

  vector<float> ref_cur = computeResCurvature(ref_points, ref_points);
  vector<float> res_cur = computeResCurvature(res_points, ref_points);

  /*************Calculate the upper reference speed limit***************/
  // Get delta_x
  vector<float> delta_x = getDeltaX(res_points);

  // Generate initial reference velocity based on delta_x and curvature
  vector<float> init_refv = generateInitialSpeeds(res_cur, delta_x);

  // Process the lower envelope
  vector<float> temp_refv = smoothLowerEnvelope(init_refv, 5);

  // Correction of reference speed based on maximum acceleration and deceleration
  std::vector<float> smooth_refv = smoothSpeeds(temp_refv, delta_x);

  std::vector<float> fit_refv = smoothCurvature(smooth_refv, 8);

  // Pull each point below the lower envelope curve
  for (int i = 0; i < fit_refv.size(); i++) {
    fit_refv[i] = std::min(fit_refv[i], smooth_refv[i]);
  }

  // The reference speed is corrected again based on the maximum acceleration and deceleration
  std::vector<float> final_refv = smoothSpeeds(fit_refv, delta_x);



  /**************Trajectory Optimization****************/
  // Generate initial reference speed profile
  // 1.Find the waypoint segment closest to the initial position (note the out-of-range situation),
  // and then project the initial velocity onto it
  Eigen::Vector3f init_pos = res_points.front();
  Eigen::Vector3f init_vel(0.1, 0, 0);
  int closest_idx = findClosestSegment(init_pos, res_points);
  Eigen::Vector3f proj_v = Rk[closest_idx] * init_vel;
  // Todo: Interpolation for more accurate positions
  float start_pos = path_points[closest_idx];
  float start_vel = proj_v.norm();
  std::cout << "start_pos: " << start_pos << ", start_vel: " << start_vel << std::endl;
  // 2.Take this v as the starting point and smoothly chase the upper limit of the reference speed.
  //  During the chasing process, the acceleration and deceleration limits should be considered to try to fit the upper
  //  limit as closely as possible.
  float max_acc = 3.0;
  float min_acc = -3.0;
  
  auto start_time2 = std::chrono::high_resolution_clock::now();
  std::vector<State> smooth_ref_traj = generateSpeedProfile(path_points, final_refv, max_acc, min_acc, start_pos, start_vel);
  auto end_time2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration2 = end_time2 - start_time2;
  std::cout << "Time taken by speed generator: " << duration2.count() << " seconds" << std::endl;
 
  // Get a smooth path starting from the initial position of the drone
  std::vector<Eigen::Vector3f> sub_res_points(res_points.begin() + closest_idx, res_points.begin() + closest_idx + HorizonNum+1);
  // Get trajectory optimization reference path points
  std::vector<Eigen::Vector3f> traj_ref_points = get_traj_ref_points(smooth_ref_traj, sub_res_points, HorizonNum+1);
  std::vector<Eigen::Vector3f> traj_ref_points1 = get_traj_ref_points(smooth_ref_traj, sub_res_points, HorizonNum+2);

  vector<float> traj_v_norm;
  for (int i = 0; i <= HorizonNum; i++) {
    traj_v_norm.push_back(smooth_ref_traj[i].velocity);
  }
  std::vector<Eigen::Vector3f> traj_ref_speed;
  for (int i = 0; i <= HorizonNum; ++i) {
    // traj_ref_speed.push_back((sub_res_points[i+1] - sub_res_points[i]).normalized()*traj_v_norm[i]);
    traj_ref_speed.push_back((traj_ref_points1[i+1] - traj_ref_points1[i]).normalized()*traj_v_norm[i]);
  }
  // traj_ref_speed.push_back((sub_res_points[HorizonNum] - sub_res_points[HorizonNum-1]).normalized()*traj_v_norm[HorizonNum]);

  std::ofstream out1("/home/alan/桌面/plot/traj_refv_x.txt");
  std::ofstream out2("/home/alan/桌面/plot/traj_refv_y.txt");
  std::ofstream out3("/home/alan/桌面/plot/traj_refv_z.txt");
  for (int i = 0; i < HorizonNum + 1; i++) {
    out1 << traj_ref_speed[i].x() << "\n";
    out2 << traj_ref_speed[i].y() << "\n";
    out3 << traj_ref_speed[i].z() << "\n";
  }
  out1.close();
  out2.close(); 
  out3.close(); 



  // get dt
  vector <float> traj_dt(HorizonNum, 0.1f);

  // match corridor
  vec_E<Polyhedron<3>> traj_mpc_polyhedrons;
  std::array<Eigen::Matrix<float, 3, 3>, HorizonNum + 1> traj_elliE;
  std::array<Eigen::Matrix<float, TrjSizeYx - TrjSizeEqx, 1>, HorizonNum + 1> traj_new_centerX;
  std::array<Eigen::Matrix<float, TrjSizeYu - TrjSizeEqu, 1>, HorizonNum + 1> traj_new_centerU;
  for (auto& mat : traj_new_centerU) {
    mat.setZero();
  }

  for (int i = 0; i <= HorizonNum; ++i) {
    int idx = findClosestSegment(traj_ref_points[i], path_f);
    traj_mpc_polyhedrons.push_back(polyhedrons[idx]);
    std::cout << "point " << i << " is on " << idx << " lane" << std::endl;
    traj_elliE[i] = E_inverse[idx];
    traj_new_centerX[i].block(0, 0, 3, 1) = E_inverse_d[idx];
  }

  // compute traj_Rk
  vector<Eigen::Matrix<float, 3, 3>> traj_Rk;
  traj_Rk.resize(HorizonNum);
  Eigen::Vector3f p0;
  Eigen::Vector3f p1;
  Eigen::Vector3f p2;
  for (int i = 0; i < HorizonNum; ++i) {
    if (i == 0) {
        p0 = traj_ref_points[0];
        p1 = traj_ref_points[1];
        p2 = traj_ref_points[2];
        traj_Rk[i] = calculateTransformationMatrix(p0, p1, p2);
    } else if (i == HorizonNum - 1) {
        p0 = traj_ref_points[HorizonNum - 2];
        p1 = traj_ref_points[HorizonNum - 1];
        p2 = traj_ref_points[HorizonNum];
        traj_Rk[i] = calculateTransformationMatrix(p0, p1, p2);
    } else {
        p0 = traj_ref_points[i-1];
        p1 = traj_ref_points[i];
        p2 = traj_ref_points[i+1];
        traj_Rk[i] = calculateTransformationMatrix(p0, p1, p2);
    }
  }

  // fomulate traj planning
  std::cout << "start solve Traj MPC" << std::endl;
  BlockVector<float, HorizonNum + 1, TrjSizeX + TrjSizeU> res_traj;
  res_traj.setZero();
  solveMpcTraj(traj_mpc_polyhedrons, traj_new_centerX, traj_new_centerU, traj_elliE, traj_dt, traj_ref_points, traj_ref_speed, traj_v_norm, traj_Rk, res_traj);
  std::cout << "end solve Traj Mpc" << std::endl;

  // save traj
  std::ofstream outfile4("/home/alan/桌面/plot/traj_time.txt");
  for (int i = 0; i < HorizonNum + 1; i++) {
    outfile4 << i*0.1f << "\n";
  }
  outfile4.close();

  std::ofstream outfile2("/home/alan/桌面/plot/traj_pos.txt");
  float start_position = sqrtf(res_traj.v[0](1, 0)*res_traj.v[0](1, 0) + res_traj.v[0](2, 0)*res_traj.v[0](2, 0) + res_traj.v[0](0, 0)*res_traj.v[0](0, 0));
  for (int i = 0; i < HorizonNum + 1; i++) {
    outfile2 << sqrtf(res_traj.v[i](1, 0)*res_traj.v[i](1, 0) + res_traj.v[i](2, 0)*res_traj.v[i](2, 0) + res_traj.v[i](0, 0)*res_traj.v[i](0, 0)) - start_position << "\n";
  }
  outfile2.close();

  std::ofstream outfile3("/home/alan/桌面/plot/traj_speed.txt");
  for (int i = 0; i < HorizonNum + 1; i++) {
    outfile3 << sqrtf(res_traj.v[i](3, 0)*res_traj.v[i](3, 0) + res_traj.v[i](4, 0)*res_traj.v[i](4, 0) + res_traj.v[i](5, 0)*res_traj.v[i](5, 0)) << "\n";
  }
  outfile3.close();

  //Publish visualization msgs
  decomp_ros_msgs::EllipsoidArray es_msg = DecompROS::ellipsoid_array_to_ros(decomp_util.get_limit_ellipsoids());
  es_msg.header.frame_id = "map";
  es_pub.publish(es_msg);

  decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(decomp_util.get_polyhedrons());
  poly_msg.header.frame_id = "map";
  poly_pub.publish(poly_msg);


  // define marker
  std::vector<geometry_msgs::Point> points_list; // smooth path
  std::vector<geometry_msgs::Point> points1_list; // ref path
  std::vector<geometry_msgs::Point> points2_list; // smooth traj
  std::vector<geometry_msgs::Point> points3_list; // ref traj
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
      geometry_msgs::Point p2;
      p2.x = res_traj.v[i](0, 0);
      p2.y = res_traj.v[i](1, 0);
      p2.z = res_traj.v[i](2, 0);
      points2_list.push_back(p2);
      geometry_msgs::Point p3;
      p3.x = traj_ref_points[i].x();
      p3.y = traj_ref_points[i].y();
      p3.z = traj_ref_points[i].z();
      points3_list.push_back(p3);
  }
  for (const auto& p : points_list) {
      v_points.points.push_back(p);
      line_strip.points.push_back(p);
  }
  for (const auto& p1 : points1_list) {
      v_points1.points.push_back(p1);
      line_strip1.points.push_back(p1);
  }
  for (const auto& p2 : points2_list) {
      v_points2.points.push_back(p2);
      line_strip2.points.push_back(p2);
  }
  for (const auto& p3 : points3_list) {
      v_points3.points.push_back(p3);
      line_strip3.points.push_back(p3);
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

  std::vector<visualization_msgs::Marker> text_markers2;
  for (size_t i = 0; i < points2_list.size(); ++i) {
      visualization_msgs::Marker text_marker2;
      text_marker2.header.frame_id = "map";
      text_marker2.header.stamp = ros::Time::now();
      text_marker2.ns = "text_marker2";
      text_marker2.id = i + 200; // 从200开始,避免ID冲突
      text_marker2.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
      text_marker2.action = visualization_msgs::Marker::ADD;
      text_marker2.pose.position = points2_list[i];
      text_marker2.pose.position.z += 0.2; // 使文本稍微高于点
      text_marker2.pose.orientation.w = 1.0;
      text_marker2.scale.z = 0.2; // 文本的高度
      text_marker2.color.r = 0.0f;
      text_marker2.color.g = 0.5f;
      text_marker2.color.b = 1.0f;
      text_marker2.color.a = 1.0;
      text_marker2.text = std::to_string(i); // 设置文本为点的ID
      text_markers2.push_back(text_marker2);
  }

  std::vector<visualization_msgs::Marker> text_markers3;
  for (size_t i = 0; i < points3_list.size(); ++i) {
      visualization_msgs::Marker text_marker3;
      text_marker3.header.frame_id = "map";
      text_marker3.header.stamp = ros::Time::now();
      text_marker3.ns = "text_marker3";
      text_marker3.id = i + 300; // 从300开始,避免ID冲突
      text_marker3.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
      text_marker3.action = visualization_msgs::Marker::ADD;
      text_marker3.pose.position = points3_list[i];
      text_marker3.pose.position.z += 0.2; // 使文本稍微高于点
      text_marker3.pose.orientation.w = 1.0;
      text_marker3.scale.z = 0.2; // 文本的高度
      text_marker3.color.r = 0.0f;
      text_marker3.color.g = 0.5f;
      text_marker3.color.b = 1.0f;
      text_marker3.color.a = 1.0;
      text_marker3.text = std::to_string(i); // 设置文本为点的ID
      text_markers3.push_back(text_marker3);
  }





  // 凸走廊结果
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

  ros::Rate r(10);
  ros::Publisher traj_pub = nh.advertise<geometry_msgs::PointStamped>("point_to_visualize", 1);
  int count = 0;

  ros::Publisher joint_pub = nh.advertise<sensor_msgs::JointState>("joint_states", 1);
  tf2_ros::TransformBroadcaster broadcaster;
  const double degree = M_PI / 180;
  // robot state
  double angle = 0;
  double radius = 2.0; // Define the radius of the circular path
  // message declarations
  sensor_msgs::JointState joint_state;
  geometry_msgs::TransformStamped odom_trans;
  odom_trans.header.frame_id = "map";
  odom_trans.child_frame_id = "base_link";

  // publish visulation msg
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

      marker_pub2.publish(v_points2);
      marker_pub2.publish(line_strip2);
      for (const auto& text_marker2 : text_markers2) {
          marker_pub2.publish(text_marker2);
      }

      marker_pub3.publish(v_points3);
      marker_pub3.publish(line_strip3);
      for (const auto& text_marker3 : text_markers3) {
          marker_pub3.publish(text_marker3);
      }

      geometry_msgs::PointStamped trj_point;
      trj_point.header.stamp = ros::Time::now();
      trj_point.header.frame_id = "map"; 
      // 设置点的位置
      trj_point.point.x = res_traj.v[count](0, 0);
      trj_point.point.y = res_traj.v[count](1, 0);
      trj_point.point.z = res_traj.v[count](2, 0);
      traj_pub.publish(trj_point);

      // urdf模型移动
      // update joint_state
      joint_state.header.stamp = ros::Time::now();
      joint_state.name.resize(1);
      joint_state.position.resize(1);
      joint_state.name[0] = "base_to_propeller";
      joint_state.position[0] = angle;

      // update transform
      odom_trans.header.stamp = ros::Time::now();
      odom_trans.transform.translation.x = res_traj.v[count](0, 0);
      odom_trans.transform.translation.y = res_traj.v[count](1, 0);
      odom_trans.transform.translation.z = res_traj.v[count](2, 0);

      tf2::Quaternion q;
      q.setRPY(0, 0, angle);  // If you want to rotate the drone as well
      odom_trans.transform.rotation.x = q.x();
      odom_trans.transform.rotation.y = q.y();
      odom_trans.transform.rotation.z = q.z();
      odom_trans.transform.rotation.w = q.w();

      // send the joint state and transform
      joint_pub.publish(joint_state);
      broadcaster.sendTransform(odom_trans);

      if (count > HorizonNum - 1) {
        count = 0;
      }
      else {
        count++;
      }
      r.sleep();
  }
  ros::spin();

  return 0;
}
