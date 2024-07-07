#pragma once
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <queue>

struct State {
    float t;
    float position;
    float velocity;
    float acceleration;
    float jerk;
};

State updateState(const State& current_state, float jerk, float dt) {
    State next_state;
    next_state.t = current_state.t + dt;
    next_state.jerk = jerk;
    next_state.acceleration = current_state.acceleration + jerk * dt;
    next_state.velocity = current_state.velocity + current_state.acceleration * dt + 0.5f * jerk * dt * dt;
    next_state.position = current_state.position + current_state.velocity * dt + 0.5f * current_state.acceleration * dt * dt + (1.0f / 6.0f) * jerk * dt * dt * dt;
    return next_state;
}


float interpolateSpeedLimit(const std::vector<float>& path, const std::vector<float>& speed_limits, float position) {
    // 使用线性插值法来获取当前位置的速度上限
    if (position <= path.front()) {
        return speed_limits.front();
    } else if (position >= path.back()) {
        return speed_limits.back();
    } else {
        for (size_t i = 0; i < path.size() - 1; ++i) {
            if (position >= path[i] && position <= path[i + 1]) {
                float t = (position - path[i]) / (path[i + 1] - path[i]);
                return speed_limits[i] + t * (speed_limits[i + 1] - speed_limits[i]);
            }
        }
    }
    return speed_limits.back(); // 默认返回最后一个速度上限
}


std::vector<std::vector<State>> sampleTrajs(const std::vector<float>& path, const std::vector<float>& speed_limits,
                                            float p0, float v0, float a0, float j0, float time_interval, float total_time) {
    int max_step = int(total_time / time_interval);
    std::vector<std::vector<State>> trajectories;
    std::queue<std::pair<int, std::vector<State>>> state_queue;
    State initial_state = {0.0f, p0, v0, a0, 0.0f};
    state_queue.push({0, {initial_state}});
    std::vector<float> jerk_vec_acc = {.0f, .2f, .4f, .6f, .8f, 1.0f, 3.0f, 5.0f};
    std::vector<float> jerk_vec_dec = {-.2f, -.4f, -.6f, -.8f, -1.0f, -3.0f, -5.0f};
    while (!state_queue.empty()) {
        auto [step, trajectory] = state_queue.front();
        state_queue.pop();
        
        if (step == max_step) {
            trajectories.push_back(trajectory);
            continue;
        }

        State current_state = trajectory.back();
        
        for (auto& j : jerk_vec_acc) {
            State next_state = updateState(current_state, j, time_interval);
            float speed_limit = interpolateSpeedLimit(path, speed_limits, next_state.position);
            if (next_state.acceleration > 2.0f || next_state.velocity > speed_limit) {
                break;
            } else if (next_state.acceleration >= -2.5f && next_state.velocity >= 0.0f) {
                std::vector<State> new_trajectory = trajectory;
                new_trajectory.push_back(next_state);
                state_queue.push({step + 1, new_trajectory});
            }
        }
        for (auto& j : jerk_vec_dec) {
            State next_state = updateState(current_state, j, time_interval);
            float speed_limit = interpolateSpeedLimit(path, speed_limits, next_state.position);
            if (next_state.acceleration < -2.5f || next_state.velocity < 0.0f) {
                break;
            } else if (next_state.acceleration <= 2.0f && next_state.velocity <= speed_limit) {
                std::vector<State> new_trajectory = trajectory;
                new_trajectory.push_back(next_state);
                state_queue.push({step + 1, new_trajectory});
            }
        }
    }
    std::cout << "get traj size: " << trajectories.size() << std::endl;
    return trajectories;
}

std::vector<State> selectTrajs(const std::vector<std::vector<State>>& sample_trajs, const std::vector<float>& path, const std::vector<float>& speed_limits, float time_interval) {
    std::vector<State> best_trajectory;
    float min_error = std::numeric_limits<float>::max();

    // Todo: 在补密中检查速度是否超限
    for (const auto& traj : sample_trajs) {
        float total_error = 0.0f;
        std::vector<State> smooth_best_traj;
        State current_state;
        float ddt = time_interval / 4.0f;
        for (int j = 0; j < traj.size() - 1; j++) {
            current_state = traj[j];
            float jerk = traj[j+1].jerk;
            smooth_best_traj.push_back(current_state);
            for (float i = 0.0f; i < time_interval - ddt; i += ddt) {
                State next_state = updateState(current_state, jerk, ddt);
                smooth_best_traj.push_back(next_state);
                current_state = next_state;
            }
        }
        smooth_best_traj.push_back(traj.back());
        for (int i = 0; i < smooth_best_traj.size() - 1; i++) {
            State front_state = smooth_best_traj[i];
            float speed_limit = interpolateSpeedLimit(path, speed_limits, front_state.position);
            float err = front_state.velocity - speed_limit;
            if (err > 0.0f) {
                total_error += 10000.0f;
            } else {
                total_error += std::fabs(err);
            }
        }

        if (total_error < min_error) {
            min_error = total_error;
            best_trajectory = traj;
        }
    }

    return best_trajectory;
}