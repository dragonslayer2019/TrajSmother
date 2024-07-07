# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd

# 加载 CSV 文件
df = pd.read_csv('/home/alan/TrajSmother/src/DecompROS/decomp_test_node/src/MPC/trajectories.csv', header=None, names=['position', 'velocity'])

# 分割不同轨迹的数据
split_indices = df.index[df['position'].isnull()].tolist()
trajectories = [df[start:end] for start, end in zip([0] + split_indices, split_indices + [None])]

# 绘制曲线图
plt.figure(figsize=(10, 6))
for trajectory in trajectories:
    plt.plot(trajectory['position'], trajectory['velocity'], marker='o')

data_fit = []
with open('/home/alan/TrajSmother/src/DecompROS/decomp_test_node/src/MPC/fit_v.txt', 'r') as file:
    for line in file:
        data_fit.append(float(line.strip()))
x_axil = []
with open('/home/alan/TrajSmother/src/DecompROS/decomp_test_node/src/MPC/x_data.txt', 'r') as file:
    for line in file:
        x_axil.append(float(line.strip()))

plt.plot(x_axil, data_fit, color='b', marker='o', linestyle='-', label=unicode('参考速度上限', 'utf-8'))

plt.xlabel('Position')
plt.ylabel('Velocity')
plt.title('Velocity vs Position for Trajectories')
plt.grid(True)
plt.legend(['Trajectory 1', 'Trajectory 2'])  # 根据实际情况添加轨迹标签
plt.show()







