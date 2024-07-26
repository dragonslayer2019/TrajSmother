#include <ros/ros.h>
#include <std_msgs/String.h>
#include <geometry_msgs/PoseArray.h>
#include "/home/alan/catkin_ws/src/global_path_search/global_path_search.h"
#include "/home/alan/catkin_ws/src/global_path_search/src/zyb_octree.cpp"
#include "/home/alan/catkin_ws/src/global_path_search/src/oct_cube.cpp"
#include "/home/alan/catkin_ws/src/global_path_search/include/a_star.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <chrono>
#include <visualization_msgs/Marker.h>

#include <iostream>
#include <memory>
#include <vector>
#include <queue>

namespace global_path_search {
class GlobalPathSearchTest {
 public:
    GlobalPathSearchTest() {
    }
    ~GlobalPathSearchTest() {}

 public:
    std::shared_ptr<Searcher> searcher_;

 public:
    void AddOctNode(OcTree::OctNode* oct_node, size_t cur_level = 0, size_t rela_to_node = 0, size_t cur_index = 0) {
        searcher_->AddOctNode(oct_node, cur_level, rela_to_node, cur_index);
        return;
    }

    NodeId GetIndex(Point3D pt, int level) {
        std::cout << "test get index" << std::endl;
        return searcher_->GetIndex(pt, level);
    }

    void HalfSolve() {
        std::cout << "wh ??" << std::endl;
        searcher_->HalfSolve(Point3D(1., 1., 1.), Point3D(127., 127., 127.));
        return;
    }

    Point3D BivariateGaussian(double x, double y, double sigma_x, double sigma_y, double mu_x, double mu_y, double norm = 62.) {
        // double norm = 62.;
        double gauss = norm * std::exp(-(std::pow(x - mu_x, 2.) / (2. * std::pow(sigma_x, 2.)) + std::pow(y - mu_y, 2.) / (2. * std::pow(sigma_y, 2.))));
        return Point3D(x, y, std::max(gauss - 12., 0.));
    }

    std::vector<Point3D> DrawMountain(double x, double y, double len, double dense, double x_ra = 1., double y_ra = 1., double norm = 62.) {
        std::vector<Point3D> res;
        double eps = 1e-5;
        double x_min = std::max(x - len / 2., eps);
        double x_max = x + len / 2.;
        double y_min = std::max(y - len / 2., eps);
        double y_max = y + len / 2.;
        for (double cur_x = x_min; cur_x < x_max + eps; cur_x += dense) {
            for (double cur_y = y_min; cur_y < y_max + eps; cur_y += dense) {
                res.push_back(BivariateGaussian(cur_x, cur_y, 32., 32., cur_x + x_ra * (cur_x - x), cur_y + y_ra * (cur_y - y), norm));
            }
        }
        for (auto pt: res) {
            std::cout << "mountain pt: " << pt.ToString() << std::endl;
        }
        return res;
    }

    std::vector<Point3D> DrawObsCube(Point3D centre_pt, double len_x, double len_y, double len_z, double dense) {
        double eps = 1e-5;
        std::vector<Point3D> res;
        double x_min = centre_pt.x - len_x / 2.;
        double x_max = centre_pt.x + len_x / 2.;
        double y_min = centre_pt.y - len_y / 2.;
        double y_max = centre_pt.y + len_y / 2.;
        double z_min = centre_pt.z - len_z / 2.;
        double z_max = centre_pt.z + len_z / 2.;
        for (double x = x_min; x < x_max + eps; x += dense) {
            for (double y = y_min; y < y_max + eps; y += dense) {
                res.push_back(Point3D(x, y, z_min));
                res.push_back(Point3D(x, y, z_max));
            }
        }
        for (double x = x_min; x < x_max + eps; x += dense) {
            for (double z = z_min; z < z_max + eps; z += dense) {
                res.push_back(Point3D(x, y_min, z));
                res.push_back(Point3D(x, y_max, z));
            }
        }
        for (double y = y_min; y < y_max + eps; y += dense) {
            for (double z = z_min; z < z_max + eps; z += dense) {
                res.push_back(Point3D(x_min, y, z));
                res.push_back(Point3D(x_max, y, z));
            }
        }
        return res;
    }

    std::vector<Point3D> DrawWall(Point3D centre_pt, double len_x, double len_y, double len_z, size_t dir, double dense) {
        double eps = 1e-5;
        std::vector<Point3D> res;
        double x_min = centre_pt.x - len_x / 2.;
        double x_max = centre_pt.x + len_x / 2.;
        double y_min = centre_pt.y - len_y / 2.;
        double y_max = centre_pt.y + len_y / 2.;
        double z_min = centre_pt.z - len_z / 2.;
        double z_max = centre_pt.z + len_z / 2.;
        if(dir == 0) {
            for (double y = y_min; y < y_max + eps; y += dense) {
                for (double z = z_min; z < z_max + eps; z += dense) {
                    res.push_back(Point3D(x_min, y, z));
                    res.push_back(Point3D(x_max, y, z));
                }
            }
            return res;
        } else if (dir == 1) {
            for (double x = x_min; x < x_max + eps; x += dense) {
                for (double z = z_min; z < z_max + eps; z += dense) {
                    res.push_back(Point3D(x, y_min, z));
                    res.push_back(Point3D(x, y_max, z));
                }
            }
            return res;
        } else {
            for (double x = x_min; x < x_max + eps; x += dense) {
                for (double y = y_min; y < y_max + eps; y += dense) {
                    res.push_back(Point3D(x, y, z_min));
                    res.push_back(Point3D(x, y, z_max));
                }
            }
            return res;
        }
        return res;
    }
};
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "global_path_search_node");
  ros::NodeHandle nh;

  ros::Publisher path_pub = nh.advertise<geometry_msgs::PoseArray>("global_path", 10);
  ros::Publisher obs_pub = nh.advertise<geometry_msgs::PoseArray>("obs_vec", 10);
  ros::Publisher mountain_pub = nh.advertise<visualization_msgs::Marker>("visualization_marker", 1);



  auto FILE = freopen("sce1_stdout", "w", stdout);
  nlohmann::json obs_vizer;
  std::ofstream os;
  os.open("obs.json");
  std::shared_ptr<global_path_search::GlobalPathSearchTest> p_test = std::make_shared<global_path_search::GlobalPathSearchTest>();
  std::shared_ptr<OcTree> A = std::make_shared<OcTree>(10);
  std::cout << "start construct obs" << std::endl;
  std::vector<Point3D> pt_vec;
  auto m_pts0 = p_test->DrawMountain(140., 140., 160., 1.35, 0.65, 1, 62.);
  auto m_pts1 = p_test->DrawMountain(80., 248., 100., 1.6, 1., 1., 62.);
  auto m_pts2 = p_test->DrawMountain(100., 200., 100., 1.35, 1., 1., 50.);
  auto pts0 = p_test->DrawObsCube({166., 216., 12.5}, 15., 10., 25., 1.6);
  auto pts1 = p_test->DrawObsCube({195., 226., 12.5}, 15., 10., 25., 1.35);
  auto pts2 = p_test->DrawObsCube({219., 226., 12.5}, 15., 10., 25., 1.35);
  auto pts3 = p_test->DrawWall({207., 226., 18.}, 10., 3.6, 2.8, 2, 1.35);
  auto pts4 = p_test->DrawObsCube({45., 76., 12.5}, 15., 10., 25., 1.35);
  auto pts5 = p_test->DrawObsCube({69., 76., 12.5}, 15., 10., 25., 1.35);
  auto pts6 = p_test->DrawWall({57., 76., 18.}, 10., 3.6, 2.8, 2, 1.35);

  // auto pts4 = p_test->DrawObsCube({80., 40., 20.01}, 6., 2., 40., 1.);
  pt_vec.insert(pt_vec.end(), m_pts0.begin(), m_pts0.end());
  pt_vec.insert(pt_vec.end(), m_pts1.begin(), m_pts1.end());
  pt_vec.insert(pt_vec.end(), m_pts2.begin(), m_pts2.end());
  pt_vec.insert(pt_vec.end(), pts0.begin(), pts0.end());
  pt_vec.insert(pt_vec.end(), pts1.begin(), pts1.end());
  pt_vec.insert(pt_vec.end(), pts2.begin(), pts2.end());
  pt_vec.insert(pt_vec.end(), pts3.begin(), pts3.end());
  pt_vec.insert(pt_vec.end(), pts4.begin(), pts4.end());
  pt_vec.insert(pt_vec.end(), pts5.begin(), pts5.end());
  pt_vec.insert(pt_vec.end(), pts6.begin(), pts6.end());
  obs_vizer[0]["pt"] = {166., 216., 12.5};
  obs_vizer[0]["len"] = {15., 10., 25.};
  // obs_vizer[1]["pt"] = std::make_tuple(40., 40., 20.);
  // obs_vizer[1]["len"] = std::make_tuple(6., 2., 40.);
  obs_vizer[1]["pt"] = {196., 226., 12.5};
  obs_vizer[1]["len"] = {15., 10., 25.};
  obs_vizer[2]["pt"] = {218., 226., 12.5};
  obs_vizer[2]["len"] = {15., 10., 25.};
  obs_vizer[3]["pt"] = {207., 226., 18.};
  obs_vizer[3]["len"] = {10., 3.6, 2.8};

  obs_vizer[4]["pt"] = {46., 76., 12.5};
  obs_vizer[4]["len"] = {15., 10., 25.};
  obs_vizer[5]["pt"] = {68., 76., 12.5};
  obs_vizer[5]["len"] = {15., 10., 25.};
  obs_vizer[6]["pt"] = {57., 76., 18.};
  obs_vizer[6]["len"] = {10., 3.6, 2.8};

  os << obs_vizer << std::endl;
  // std::cout << "check wa ??" << std::endl;

  A->MainUpdate(A->GetRoot(), {}, pt_vec);
  std::cout << "main update finish" << std::endl;
  std::vector<size_t> level_leaf_node_cnt(13);
  std::vector<size_t> level_node_cnt(13);
  auto all_nodes = A->GetAllLeafNodes();
  for(auto node_pair: all_nodes) {
      auto oct_node = node_pair.first;
      if (oct_node->level == 1) {
          std::cout << "LEVEL 1 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
          std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
      }
      if (oct_node->level == 2) {
          std::cout << "LEVEL 2 oct node: index = " << oct_node->index << " , Point = " << oct_node->CentralPoint().ToString() << std::endl;
          std::cout << "oct node cost: " << oct_node->cost_min << " ," <<oct_node->cost_max << std::endl;
          std::cout << "is leaf = " << node_pair.second << std::endl;
      }
      if(node_pair.second == true) {
          level_leaf_node_cnt[node_pair.first->level] ++;
      } else {

      }
      level_node_cnt[node_pair.first->level] ++;
  }
  for (size_t i = 0; i < 12; i++) {
      std::cout << "LEVEL : " << i << " leaf cnt = " << level_leaf_node_cnt[i] << std::endl;
      std::cout << "LEVEL : " << i << " cnt = " << level_node_cnt[i] << std::endl;
  }

  // auto FILE = freopen("sce1_stdout", "w", stdout);
  std::cout << "start running" << std::endl;
  auto start_ = std::chrono::high_resolution_clock::now();
  p_test->searcher_ = std::make_shared<global_path_search::Searcher> (A);
  p_test->searcher_->HalfSolve(Point3D(57., 76., 12.), Point3D(207., 226., 12.));
  p_test->searcher_->AstarTest();
  p_test->searcher_->PostProcessing();
  auto end_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_cost = end_ - start_;
  std::cout << "searching time: " << time_cost.count() << std::endl;
  std::cout << "running finish, any key to continue" << std::endl;



  // 画出障碍物点
  visualization_msgs::Marker v_points;
  v_points.header.frame_id = "map";
  v_points.header.stamp = ros::Time::now();
  v_points.ns = "mountain";
  v_points.action = visualization_msgs::Marker::ADD;
  v_points.pose.orientation.w = 1.0;
  v_points.id = 2;
  v_points.type = visualization_msgs::Marker::POINTS;
  v_points.scale.x = v_points.scale.y = v_points.scale.z = 0.2; // 点的大小
  v_points.color.r = 0.5;
  v_points.color.g = 0.5;
  v_points.color.b = 0.5;
  v_points.color.a = 1.0;
  for (const auto& pt : pt_vec) {
    geometry_msgs::Point p;
    p.x = pt.x;
    p.y = pt.y;
    p.z = pt.z;
    v_points.points.push_back(p);
  }


  std::vector<Point3D> final_global_path = p_test->searcher_->GetFinalWaypts();
  // 获取障碍物点
  std::vector<Point3D> obs_vec;
  for (int i = 0; i <final_global_path.size() - 1; i++) {
    std::vector<Point3D> temp_obs = A->GetObs(1.0, 2.0, 1.0, final_global_path[i], final_global_path[i+1]);
    obs_vec.insert(obs_vec.end(), temp_obs.begin(), temp_obs.end());
  }

  std::vector<geometry_msgs::Pose> path;
  std::vector<geometry_msgs::Pose> obstacles;

  /***********补密*************/
    std::vector<Point3D> densePolyline;
    double stepSize = 10.0;
    for (size_t i = 0; i < final_global_path.size() - 1; ++i) {
        Point3D p1 = final_global_path[i];
        Point3D p2 = final_global_path[i + 1];
        
        densePolyline.push_back(p1);
        
        double dist = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2) + std::pow(p2.z - p1.z, 2));
        int numSteps = std::max(1, static_cast<int>(std::round(dist / stepSize)));

        for (int j = 1; j < numSteps; ++j) {
            double t = static_cast<double>(j) / numSteps;
            double x = p1.x + t * (p2.x - p1.x);
            double y = p1.y + t * (p2.y - p1.y);
            double z = p1.z + t * (p2.z - p1.z);
            densePolyline.push_back(Point3D(x, y, z));
        }
    }

    // Add the last point
    densePolyline.push_back(final_global_path.back());
  /************************/
  /***********画全局点补密*************/
    std::vector<Point3D> denseDrawPolyline;
    double stepDrawSize = 2.0;
    for (size_t i = 0; i < final_global_path.size() - 1; ++i) {
        Point3D p1 = final_global_path[i];
        Point3D p2 = final_global_path[i + 1];
        
        denseDrawPolyline.push_back(p1);
        
        double dist = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2) + std::pow(p2.z - p1.z, 2));
        int numSteps = std::max(1, static_cast<int>(std::round(dist / stepDrawSize)));

        for (int j = 1; j < numSteps; ++j) {
            double t = static_cast<double>(j) / numSteps;
            double x = p1.x + t * (p2.x - p1.x);
            double y = p1.y + t * (p2.y - p1.y);
            double z = p1.z + t * (p2.z - p1.z);
            denseDrawPolyline.push_back(Point3D(x, y, z));
        }
    }

    // Add the last point
    densePolyline.push_back(final_global_path.back());
    denseDrawPolyline.push_back(final_global_path.back());
  /************************/

  ros::Publisher marker_pub = nh.advertise<visualization_msgs::Marker>("vis_marker", 10);
  visualization_msgs::Marker glob_pts;
  glob_pts.header.frame_id = "map";         // 设置坐标系
  glob_pts.header.stamp = ros::Time::now(); // 设置时间戳
  glob_pts.ns = "global_path";
  glob_pts.id = 0;
  glob_pts.type = visualization_msgs::Marker::POINTS;
  glob_pts.action = visualization_msgs::Marker::ADD;
  glob_pts.scale.x = glob_pts.scale.y = glob_pts.scale.z = 0.2; // 点的大小
  glob_pts.color.r = 1.0;
  glob_pts.color.g = 0.2;
  glob_pts.color.b = 0.2;
  glob_pts.color.a = 1.0;
  glob_pts.pose.orientation.w = 1.0;

  for (const auto& gp : densePolyline) {
    geometry_msgs::Pose pose;
    pose.position.x = gp.x;
    pose.position.y = gp.y;
    pose.position.z = gp.z;
    path.push_back(pose);
    // geometry_msgs::Point p;
    // p.x = gp.x;
    // p.y = gp.y;
    // p.z = gp.z;
    // glob_pts.points.push_back(p);
  }
  for (const auto& gp : denseDrawPolyline) {
    geometry_msgs::Point p;
    p.x = gp.x;
    p.y = gp.y;
    p.z = gp.z;
    glob_pts.points.push_back(p);
  }  
  for (const auto& obs : obs_vec) {
    geometry_msgs::Pose pose;
    pose.position.x = obs.x;
    pose.position.y = obs.y;
    pose.position.z = obs.z;
    obstacles.push_back(pose);   
  }
  // 将路径发布为PoseArray消息
  geometry_msgs::PoseArray path_msg;
  path_msg.header.frame_id = "map";  // 设置坐标系
  path_msg.header.stamp = ros::Time::now();  // 设置时间戳
  path_msg.poses = path;

  // 将障碍物发布为PoseArray消息
  geometry_msgs::PoseArray obs_msg;
  obs_msg.header.frame_id = "map";  // 设置坐标系
  obs_msg.header.stamp = ros::Time::now();  // 设置时间戳
  obs_msg.poses = obstacles;
  ROS_INFO("obs size: %d", int(obstacles.size()));

  ros::Rate loop_rate(3);

  while (ros::ok()) {
    mountain_pub.publish(v_points);
    path_pub.publish(path_msg);
    marker_pub.publish(glob_pts);
    obs_pub.publish(obs_msg);
    loop_rate.sleep();
  }
  ros::spin();
  return 0;
}
