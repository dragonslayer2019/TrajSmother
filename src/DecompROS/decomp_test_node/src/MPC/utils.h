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
using namespace std;

float distance_2d(const Eigen::Vector2f& p1, const Eigen::Vector2f& p2) {
    return sqrt(pow(p1(0) - p2(0), 2) + pow(p1(1) - p2(1), 2));
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