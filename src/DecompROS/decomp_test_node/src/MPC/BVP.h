#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include "presteps.h"
#include "spline.h" // Spline库

struct State {
    float t;
    float position;
    float velocity;
    float acceleration;
    float jerk;
};

State interplotTraj(float cur_t, std::vector<State> traj) {
    State inter_state;
    for (int i = 1; i < traj.size(); i++) {
        if (cur_t <= traj[i].t) {
            inter_state.t = cur_t;
            float par = (cur_t - traj[i-1].t) / (traj[i].t - traj[i-1].t);
            inter_state.position = traj[i-1].position + (traj[i].position - traj[i-1].position) * par;
            inter_state.velocity = traj[i-1].velocity + (traj[i].velocity - traj[i-1].velocity) * par;
            inter_state.acceleration = traj[i-1].acceleration + (traj[i].acceleration - traj[i-1].acceleration) * par;
            inter_state.jerk = traj[i-1].jerk + (traj[i].jerk - traj[i-1].jerk) * par;
            break;
        }
    }
    return inter_state;
}

// 生成平滑速度曲线的函数
std::vector<State> generateSpeedProfile(std::vector<float>& x_data, std::vector<float>& y_data, float max_acc,
                                        float min_acc, float init_x, float init_y) {
    for (auto& y : y_data) {
        y *= 0.64f;
    }

    int start_id = 0;
    if (init_x < x_data.front()) {
        start_id = 0;
    } else if (init_x > x_data.back()) {
        start_id = x_data.size() - 1;
    } else {
        for (int i = 0; i < x_data.size() - 1; i++)
        if (init_x >= x_data[i] && init_x < x_data[i+1]) {
            start_id = i;
        }
    }

    std::vector<float> new_x;
    new_x.push_back(init_x);
    std::vector<float> new_y;
    new_y.push_back(init_y);

    for (int i = start_id + 1; i < start_id+7 && i < y_data.size(); i++) {
        if (new_y.back() - y_data[i] < -1.0f) {
            new_x.push_back(x_data[i]);
            new_y.push_back(new_y.back() + (y_data[i] - new_y.back()) * (0.1f + 0.2f*(i-start_id-1)));
        }
    }

    for (int i = start_id+7; i < x_data.size() - 1; i+=3) {
        new_x.push_back(x_data[i]);
    }

    for (int i = start_id+7; i < y_data.size() - 1; i+=3) {
        new_y.push_back(y_data[i]);
    }
    // 拟合样条曲线
    tk::spline s;
    s.set_points(new_x, new_y);

    // 生成速度曲线
    std::vector<float> smooth_speeds;
    float dx = 0.1; // 生成点的步长
    std::ofstream outfile("/home/alan/桌面/plot/final_t.txt");
    std::ofstream outfile0("/home/alan/桌面/plot/speed.txt");
    // std::cout << "temp_x: " << std::endl; 
    std::vector<float> pos;
    for (float x = new_x.front(); x <= new_x.back(); x += dx) {
        // std::cout << x << std::endl;
        outfile << x << "\n";
        outfile0 << s(x) << "\n";
        pos.push_back(x);
        smooth_speeds.push_back(s(x));
    }
    outfile.close();
    outfile0.close();

    // 保存轨迹
    std::vector<State> traj(smooth_speeds.size());
    traj[0] = {0.0, init_x, init_y, 0.0, 0.0};


    // 修正加速度限制
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float acc = (smooth_speeds[i]*smooth_speeds[i] - smooth_speeds[i - 1]*smooth_speeds[i - 1]) / 2.0f / dx;
        if (acc > max_acc) {
            smooth_speeds[i] = sqrtf(std::max(smooth_speeds[i - 1] * smooth_speeds[i - 1] + 2 * max_acc * dx, 0.1f));
        } else if (acc < min_acc) {
            smooth_speeds[i] = sqrtf(std::max(smooth_speeds[i - 1] * smooth_speeds[i - 1] + 2 * min_acc * dx, 0.1f));
        }
        traj[i].velocity = smooth_speeds[i]; // save speed
    }

    for (int i = 1; i < pos.size(); i++) {
        traj[i].position = pos[i]; // save position
    }

    // 计算时间
    std::vector<float> time = {0.0f};
    std::ofstream outfile1("/home/alan/桌面/plot/tx.txt");
    outfile1 << 0.0f << "\n";
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float average_speed = (smooth_speeds[i] + smooth_speeds[i - 1]) / 2.0f;
        float real_t = dx / average_speed;
        outfile1 << time.back() + real_t << "\n";
        time.push_back(time.back() + real_t);
        traj[i].t = time[i]; // save time
    }
    outfile1.close();
    
    // 计算加速度
    std::vector<float> accelerations;
    accelerations.push_back(0.0f); // 初始点加速度为0
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float acc = (smooth_speeds[i] - smooth_speeds[i - 1]) / (time[i] - time[i - 1]);
        accelerations.push_back(acc);
        traj[i].acceleration = accelerations[i]; // save acc
    }

    // 将加速度保存到文件
    std::ofstream outfile2("/home/alan/桌面/plot/acceleration.txt");
    for (const auto& acc : accelerations) {
        outfile2 << acc << "\n";
    }
    outfile2.close();

    // 计算jerk
    std::vector<float> jerks;
    jerks.push_back(0.0f); // 初始点jerk为0
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float jerk = (accelerations[i] - accelerations[i - 1]) / (time[i] - time[i - 1]);
        jerks.push_back(jerk);
        traj[i].jerk = jerks[i]; // save jerk
    }

    // 将jerk保存到文件
    std::ofstream outfile3("/home/alan/桌面/plot/jerks.txt");
    for (const auto& jerk : jerks) {
        outfile3 << jerk << "\n";
    }
    outfile3.close();

    float dt = 0.1;
    for (size_t i = 1; i < traj.size(); ++i) {
        float cur_t = i * dt;
        traj[i] = interplotTraj(cur_t, traj);
    }

    return traj;
}

// 获取参考路径点
std::vector<Eigen::Vector3f> get_traj_ref_points(std::vector<State>& smooth_ref_traj, vector<Eigen::Vector3f> res_points, int step) {
    std::vector<Eigen::Vector3f> traj_ref_points;
    traj_ref_points.push_back(res_points.front());
    for (size_t i = 1; i < step; ++i) {
        float length = smooth_ref_traj[i].position;
        Eigen::Vector3f unit_point;
        float total_seg_length = 0.0f;
        for (int j = 1; j < res_points.size(); j++) {
            float segment_length = distance(res_points[j], res_points[j - 1]);
            total_seg_length += segment_length;
            if (length <= total_seg_length) {
                float par = (length - total_seg_length + segment_length) / segment_length;
                unit_point.x() = res_points[j-1].x() + (res_points[j].x() - res_points[j-1].x()) * par;
                unit_point.y() = res_points[j-1].y() + (res_points[j].y() - res_points[j-1].y()) * par;
                unit_point.z() = res_points[j-1].z() + (res_points[j].z() - res_points[j-1].z()) * par;
                break;
            }
        }
        traj_ref_points.push_back(unit_point);
    }
    return traj_ref_points;
}