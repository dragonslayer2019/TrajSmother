#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "spline.h" // Spline库

// 生成平滑速度曲线的函数
std::vector<float> generateSpeedProfile(std::vector<float>& x_data, std::vector<float>& y_data, float max_acc,
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

    for (int i = start_id+5; i < x_data.size() - 1; i+=3) {
        new_x.push_back(x_data[i]);
    }

    for (int i = start_id+5; i < y_data.size() - 1; i+=3) {
        new_y.push_back(y_data[i]);
    }
    // 拟合样条曲线
    tk::spline s;
    s.set_points(new_x, new_y);

    // 生成速度曲线
    std::vector<float> smooth_speeds;
    float dx = 0.1; // 生成点的步长
    std::ofstream outfile("/home/alan/桌面/plot/final_t.txt");
    // std::cout << "temp_x: " << std::endl;    
    for (float x = new_x.front(); x <= new_x.back(); x += dx) {
        // std::cout << x << std::endl;
        outfile << x << "\n";
        smooth_speeds.push_back(s(x));     
    }
    outfile.close();

    // 修正加速度限制
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float acc = (smooth_speeds[i]*smooth_speeds[i] - smooth_speeds[i - 1]*smooth_speeds[i - 1]) / 2.0f / dx;
        if (acc > max_acc) {
            smooth_speeds[i] = sqrtf(std::max(smooth_speeds[i - 1] * smooth_speeds[i - 1] + 2 * max_acc * dx, 0.1f));
        } else if (acc < min_acc) {
            smooth_speeds[i] = sqrtf(std::max(smooth_speeds[i - 1] * smooth_speeds[i - 1] + 2 * min_acc * dx, 0.1f));
        }
    }

    std::vector<float> time = {0.0f};
    std::ofstream outfile1("/home/alan/桌面/plot/tx.txt");
    outfile1 << 0.0f << "\n";
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float average_speed = (smooth_speeds[i] + smooth_speeds[i - 1]) / 2.0f;
        float real_t = dx / average_speed;
        outfile1 << time.back() + real_t << "\n";
        time.push_back(time.back() + real_t);
    }
    outfile1.close();
    
    // 计算加速度
    std::vector<float> accelerations;
    accelerations.push_back(0.0f); // 初始点加速度为0
    for (size_t i = 1; i < smooth_speeds.size(); ++i) {
        float acc = (smooth_speeds[i] - smooth_speeds[i - 1]) / (time[i] - time[i - 1]);
        accelerations.push_back(acc);
    }

    // 将加速度保存到文件
    std::ofstream outfile2("/home/alan/桌面/plot/acceleration.txt");
    for (const auto& acc : accelerations) {
        outfile2 << acc << "\n";
    }
    outfile2.close();

    return smooth_speeds;
}