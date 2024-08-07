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
#include "utils.h"
#include <geometry_msgs/PointStamped.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/JointState.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_ros/transform_broadcaster.h>
#include "BVP.h"
#include <geometry_msgs/PoseArray.h>


int main(int argc, char ** argv){
  ros::init(argc, argv, "test");
  ros::NodeHandle nh("~");

  ros::Publisher cloud_pub = nh.advertise<sensor_msgs::PointCloud>("cloud", 1, true);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("path", 1, true);
  ros::Publisher ref_path_pub = nh.advertise<nav_msgs::Path>("ref_path_path", 1, true);
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

  vec_Vec3f obs = DecompROS::cloud_to_vec(cloud);
  vec_Vec2f obs2d;
  for(const auto& it: obs)
    obs2d.push_back(it.topRows<2>());

  //Read path from txt
  vec_Vec2f path;
  if(!read_path<2>(path_file, path))
    ROS_ERROR("Fail to read a path!");

  // check path length
  float min_length = HorizonNum * max_length;
  float length = 0;
  for (int i = 0; i < path.size()-1; i++) {
    length += distance_2d(path[i].cast<float>(), path[i+1].cast<float>());
  }
  if (length < min_length) {
    // Add a point to make the path long enough
    int back_id = path.size() - 1; 
    float delta = 1.0f;
    Eigen::Vector2f direction = (path[back_id].cast<float>() - path[back_id - 1].cast<float>()).normalized();
    Eigen::Vector2d added_point = path[back_id] + (direction * (min_length - length + delta)).cast<double>();
    std::cout << "back id: " << back_id << "dir: " << direction << std::endl;
    // 延长point
    path.push_back(added_point);
  }

  nav_msgs::Path path_msg = DecompROS::vec_to_path(path);
  path_msg.header.frame_id = "map";
  path_pub.publish(path_msg);

  //Using ellipsoid decomposition
  EllipsoidDecomp2D decomp_util;
  decomp_util.set_obs(obs2d);
  decomp_util.set_local_bbox(Vec2f(1, 2));
  decomp_util.dilate(path); //Set max iteration number of 10, do fix the path

  std::vector<Eigen::Vector2f> path_f;
  path_f.reserve(path.size());
  for (const auto& point : path) {
      path_f.emplace_back(point.cast<float>());
  }

  vector<Eigen::Vector2f> ref_points = get_sample_point_2d(path_f);
  for (const auto& pt : ref_points) {
    std::cout << pt.x() << " " << pt.y() << std::endl;
  }

  // 参考路径可视化
  vec_Vec2f ref_path(HorizonNum+1);
  for (int i = 0; i <= HorizonNum; i++) {
    ref_path[i] = ref_points[i].cast<double>();
  }
  nav_msgs::Path ref_path_msg = DecompROS::vec_to_path(ref_path);
  ref_path_msg.header.frame_id = "map";
  ref_path_pub.publish(ref_path_msg);

  vector<float> v_norm, dt;
  vector<Eigen::Matrix<float, 2, 2>> Rk;
  get_param_2d(ref_points, v_norm, Rk, dt);

  vec_E<Polyhedron<3>> polyhedrons = decomp_util.get_polyhedrons();
  vec_E<Ellipsoid<3>> ellipsoids = decomp_util.get_limit_ellipsoids();
  vec_E<Polyhedron<3>> mpc_polyhedrons;  


















  //Publish visualization msgs
  decomp_ros_msgs::EllipsoidArray es_msg = DecompROS::ellipsoid_array_to_ros(decomp_util.get_ellipsoids());
  es_msg.header.frame_id = "map";
  es_pub.publish(es_msg);

  decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(decomp_util.get_polyhedrons());
  poly_msg.header.frame_id = "map";
  poly_pub.publish(poly_msg);


  //Convert to inequality constraints Ax < b
  auto polys = decomp_util.get_polyhedrons();
  // for(size_t i = 0; i < path.size() - 1; i++) {
  //   const auto pt_inside = (path[i] + path[i+1]) / 2;
  //   LinearConstraint2D cs(pt_inside, polys[i].hyperplanes());
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


  ros::spin();

  return 0;
}
