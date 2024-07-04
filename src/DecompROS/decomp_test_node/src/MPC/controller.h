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
};

State updateState(const State& current_state, float jerk, float dt) {
    State next_state;
    next_state.acceleration = current_state.acceleration + jerk * dt;
    next_state.velocity = current_state.velocity + current_state.acceleration * dt + 0.5f * jerk * dt * dt;
    next_state.position = current_state.position + current_state.velocity * dt + 0.5f * current_state.acceleration * dt * dt + (1.0f / 6.0f) * jerk * dt * dt * dt;
    return next_state;
}

std::vector<State> generateTrajectory(
    const std::vector<float>& path,
    const std::vector<float>& speed_limits,
    float p0, float v0, float a0, float j0,
    float time_interval, int max_intervals) {
    
    std::vector<State> trajectory;
    State current_state = {p0, v0, a0};
    
    for (int i = 0; i < max_intervals; ++i) {
        // 插值当前位置的速度上限
        float speed_limit = 0;
        if (current_state.position <= path.front()) {
            speed_limit = speed_limits.front();
        } else if (current_state.position >= path.back()) {
            speed_limit = speed_limits.back();
        } else {
            for (size_t j = 1; j < path.size(); ++j) {
                if (current_state.position < path[j]) {
                    float t = (current_state.position - path[j - 1]) / (path[j] - path[j - 1]);
                    speed_limit = speed_limits[j - 1] + t * (speed_limits[j] - speed_limits[j - 1]);
                    break;
                }
            }
        }
        
        // 根据速度差调整加加速度
        float velocity_error = speed_limit - current_state.velocity;
        float jerk = std::max(-5.0f, std::min(5.0f, velocity_error / time_interval));

        // 更新状态
        State next_state = updateState(current_state, jerk, time_interval);
        
        // 速度和加速度的限制
        next_state.acceleration = std::max(-2.5f, std::min(2.0f, next_state.acceleration));
        next_state.velocity = std::max(0.0f, std::min(speed_limit, next_state.velocity));
        
        // 将当前状态加入轨迹
        trajectory.push_back(next_state);
        
        // 更新当前状态
        current_state = next_state;
        
        // 路径超出最大长度则退出
        if (current_state.position >= path.back()) break;
    }
    
    return trajectory;
}

// int main() {
//     std::vector<float> path = {0, 10, 20, 30, 40, 50};
//     std::vector<float> speed_limits = {5, 10, 15, 20, 25, 30};
    
//     float p0 = 0;
//     float v0 = 0;
//     float a0 = 0;
//     float j0 = 1;
//     float time_interval = 0.1f;
//     int max_intervals = 1000;
    
//     std::vector<State> trajectory = generateTrajectory(path, speed_limits, p0, v0, a0, j0, time_interval, max_intervals);
    
//     for (const auto& state : trajectory) {
//         std::cout << "Position: " << state.position << ", Velocity: " << state.velocity << ", Acceleration: " << state.acceleration << std::endl;
//     }
    
//     return 0;
// }