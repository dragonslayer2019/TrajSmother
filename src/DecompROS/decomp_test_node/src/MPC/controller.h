// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <algorithm>

// struct State {
//     float position;
//     float velocity;
//     float acceleration;
//     // float jerk;
// };

// class PathPlanner {
// public:
//     PathPlanner(const std::vector<float>& path_points, const std::vector<float>& speed_limits, float v0, float a0)
//         : path_points(path_points), speed_limits(speed_limits), dt(0.1) {
//         state = {0.0f, v0, a0};
//     }

//     void calculateReferenceSpeedProfile() {
//         std::vector<State> trajectory;
//         trajectory.push_back(state);
//         std::vector<float> jerk_history;
//         size_t step = 0;

//         while (step < max_steps && state.position < path_points.back()) {
//             bool valid_jerk_found = false;
//             for (float jerk = 5.0f; jerk >= -5.0f; jerk -= 1.0f) {
//                 State next_state = calculateNextState(state, jerk);
//                 if (isValidState(next_state, step)) {
//                     state = next_state;
//                     trajectory.push_back(state);
//                     jerk_history.push_back(jerk);
//                     valid_jerk_found = true;
//                     break;
//                 }
//             }

//             if (!valid_jerk_found) {
//                 if (step > 0) {
//                     trajectory.pop_back();
//                     state = trajectory.back();
//                     jerk_history.pop_back();
//                     step--;

//                     bool previous_valid_jerk_found = false;
//                     float last_jerk = jerk_history.back();
//                     for (float jerk = last_jerk - 1.0f; jerk >= -5.0f; jerk -= 1.0f) {
//                         State prev_state = calculateNextState(trajectory.back(), jerk);
//                         if (isValidState(prev_state, step)) {
//                             state = prev_state;
//                             trajectory.back() = state;
//                             jerk_history.back() = jerk;
//                             previous_valid_jerk_found = true;
//                             break;
//                         }
//                     }

//                     if (!previous_valid_jerk_found) {
//                         std::cerr << "No valid jerk value found for previous step" << std::endl;
//                         break;
//                     }
//                 } else {
//                     std::cerr << "No valid jerk value found at step 0" << std::endl;
//                     break;
//                 }
//             } else {
//                 step++;
//             }
//         }

//         printTrajectory(trajectory);
//     }

// private:
//     std::vector<float> path_points;
//     std::vector<float> speed_limits;
//     State state;
//     const float dt;
//     const size_t max_steps = 50;

//     State calculateNextState(const State& current_state, float jerk) {
        // State next_state;
        // // next_state.jerk = current_state.jerk + jerk * dt;
        // // next_state.acceleration = current_state.acceleration + current_state.jerk * dt + 0.5f * jerk * dt * dt;
        // // next_state.velocity = current_state.velocity + current_state.acceleration * dt + 0.5f * current_state.jerk * dt * dt + (1.0f / 6.0f) * jerk * dt * dt * dt;
        // // next_state.position = current_state.position + current_state.velocity * dt + 0.5f * current_state.acceleration * dt * dt + (1.0f / 6.0f) * current_state.jerk * dt * dt * dt + (1.0f / 24.0f) * jerk * dt * dt * dt * dt;
        // next_state.acceleration = current_state.acceleration + jerk * dt;
        // next_state.velocity = current_state.velocity + current_state.acceleration * dt + 0.5f * jerk * dt * dt;
        // next_state.position = current_state.position + current_state.velocity * dt + 0.5f * current_state.acceleration * dt * dt + (1.0f / 6.0f) * jerk * dt * dt * dt;
        // return next_state;
//     }

//     bool isValidState(const State& state, size_t step) {
//         if (state.acceleration < -2.5f || state.acceleration > 2.0f) return false;
//         float speed_limit = 0;
//         for (int i = 0; i < path_points.size() - 1; i++) {
//             if (state.position >= path_points[i] && state.position < path_points[i+1]) {
//                 speed_limit = (state.position - path_points[i]) / (path_points[i+1] - path_points[i]) * (speed_limits[i+1] - speed_limits[i]) + speed_limits[i];
//             }
//         }
//         if (state.velocity < 0 || state.velocity > speed_limit) return false;
//         if (state.position > path_points.back()) return false;
//         return true;
//     }

//     void printTrajectory(const std::vector<State>& trajectory) {
//         for (const auto& state : trajectory) {
//             std::cout << "Position: " << state.position << ", Velocity: " << state.velocity
//                       << ", Acceleration: " << state.acceleration << std::endl;
//         }
//     }
// };

// // int main() {
// //     std::vector<float> path_points = {0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f};
// //     std::vector<float> speed_limits = {10.0f, 15.0f, 20.0f, 15.0f, 10.0f, 5.0f};
// //     float v0 = 0.0f;
// //     float a0 = 0.0f;
// //     float j0 = 0.0f;

// //     PathPlanner planner(path_points, speed_limits, v0, a0, j0);
// //     planner.calculateReferenceSpeedProfile();

// //     return 0;
// // }
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

struct State {
    float position;
    float velocity;
    float acceleration;
    float jerk;
};


State updateState(const State& current_state, float jerk, float dt) {
    State next_state;
    next_state.jerk = jerk;
    next_state.acceleration = current_state.acceleration + jerk * dt;
    next_state.velocity = current_state.velocity + current_state.acceleration * dt + 0.5f * jerk * dt * dt;
    next_state.position = current_state.position + current_state.velocity * dt + 0.5f * current_state.acceleration * dt * dt + (1.0f / 6.0f) * jerk * dt * dt * dt;
    return next_state;
}

float interpolateSpeedLimit(float position, const std::vector<float>& path, const std::vector<float>& speed_limits) {
    if (position <= path.front()) {
        return speed_limits.front();
    } else if (position >= path.back()) {
        return speed_limits.back();
    } else {
        for (size_t j = 1; j < path.size(); ++j) {
            if (position < path[j]) {
                float t = (position - path[j - 1]) / (path[j] - path[j - 1]);
                return speed_limits[j - 1] + t * (speed_limits[j] - speed_limits[j - 1]);
            }
        }
    }
    return speed_limits.back(); // default case, shouldn't be reached
}

// std::vector<State> generateTrajectory(
//     const std::vector<float>& path,
//     const std::vector<float>& speed_limits,
//     float p0, float v0, float a0, float j0,
//     float time_interval, int max_intervals) {
    
//     std::vector<State> trajectory;
//     State current_state = {p0, v0, a0, 0};
//     float previous_jerk = 0;
//     float max_jerk_change = 1.0f; // 限制jerk的变化速率
    
//     for (int i = 0; i < max_intervals; ++i) {
//         // 插值当前位置的速度上限
//         float speed_limit = interpolateSpeedLimit(current_state.position, path, speed_limits);

//         // 根据速度差调整加加速度，确保速度不会超过速度上限
//         float velocity_error = speed_limit - current_state.velocity;
//         float desired_jerk = std::max(-5.0f, std::min(5.0f, (velocity_error / time_interval - current_state.acceleration) / time_interval));
        
//         // 限制jerk的变化速率
//         float jerk = std::max(previous_jerk - max_jerk_change, std::min(previous_jerk + max_jerk_change, desired_jerk));
        
//         // 暂时计算下一个状态
//         State next_state = updateState(current_state, jerk, time_interval);
        
//         // 如果计算出的速度超过了速度上限，则调整jerk
//         if (next_state.velocity > speed_limit) {
//             float required_jerk = (2 * (speed_limit - current_state.velocity) - current_state.acceleration * time_interval) / (time_interval * time_interval);
//             required_jerk = std::max(-5.0f, std::min(5.0f, required_jerk));
//             jerk = std::max(previous_jerk - max_jerk_change, std::min(previous_jerk + max_jerk_change, required_jerk));
//             next_state = updateState(current_state, jerk, time_interval);
//         }
        
//         // 确保速度和加速度在限制范围内
//         next_state.acceleration = std::max(-2.5f, std::min(2.0f, next_state.acceleration));
//         next_state.velocity = std::max(0.0f, std::min(speed_limit, next_state.velocity));
        
//         // 将当前状态加入轨迹
//         trajectory.push_back(next_state);
        
//         // 更新当前状态和jerk
//         current_state = next_state;
//         previous_jerk = jerk;
        
//         // 路径超出最大长度则退出
//         if (current_state.position >= path.back()) break;
//     }
    
//     return trajectory;
// }

// struct State {
//     float t;
//     float position;
//     float velocity;
//     float acceleration;
//     float jerk;
// };

// // 采样n条轨迹
// std::vector<std::vector<state>> sampleTrajs(const std::vector<float>& path, const std::vector<float>& speed_limits,float p0, float v0, float a0, float j0, float time_interval, int max_intervals) {
//     std::vector<std::vector<state>> trajectorys;
//     std::vector<State> trajectory;
//     State current_state = {p0, v0, a0, 0};
//     for (float t = 0.0f; t < 5.0f; t+=time_interval) {
        
//         for (float j = 0.0f; j <= 5.0f; j+=1.0f) {
//             // 获取下一时刻状态
//             State next_state = updateState(current_state, j, time_interval);
//             // 检查是否超出参考速度上界或加速度上界

//             // 不超出则存入轨迹

//             // 如果超出则不再继续采样
//         }
//         for (float j = -1.0f; j >= -5.0f; j-=1.0f) {
//             // 获取下一时刻状态

//             // 检查是否超出参考速度或加速度下界

//             // 不超出则存入轨迹

//             // 如果超出则不再继续采样
//         }
//     }
// }


void getTraj(const std::vector<float>& path, const std::vector<float>& speed_limits, 
             State& init_state, State& cur_state, std::vector<State>& traj, float dt, int max_step, std::vector<std::vector<State>>& trajs, int step, int count) {
    step++;
    if (step == max_step) {
        step = 0;
        trajs.push_back(traj);
        traj.clear();
        cur_state = init_state;
        count++;
        std::cout << "get traj " << count << std::endl;
    }
    traj.push_back(cur_state);

    for (float j = -5.0f; j <= 5.0f; j+=1.0f) {
        // 计算下一时刻状态
        State next_state = updateState(cur_state, j, dt);
        float speed_limit = interpolateSpeedLimit(next_state.position, path, speed_limits);
        if (next_state.velocity > 0.0f && next_state.velocity < speed_limit && next_state.acceleration > -2.5f && next_state.acceleration < 2.0f) {
            getTraj(path , speed_limits, init_state, next_state, traj, dt, max_step, trajs, step, count);
        }
    }
}

std::vector<std::vector<State>> sampleTrajs(const std::vector<float>& path, const std::vector<float>& speed_limits,
                                            float p0, float v0, float a0, float j0, float time_interval, float total_time) {
    int count = 0;
    int step = 0;
    int max_step = int(total_time / time_interval);
    std::vector<std::vector<State>> trajectories;
    std::vector<State> trajectory;
    State current_state = {p0, v0, a0, 0};
    getTraj(path, speed_limits, current_state, current_state, trajectory, time_interval, max_step, trajectories, step, count);
    return trajectories;
}

// 删除掉加速度与速度不满足约束的轨迹
// std::vector<std::vector<state>> selectTrajs(std::vector<std::vector<state>> sample_trajs) {
//     for (auto& traj : sample_trajs) {
//         // 计算每条轨迹距离参考速度上界的最小误差

//     }
//     // 选出误差最小的轨迹

// }