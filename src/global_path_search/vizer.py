import numpy as np
import cv2 as cv
import json
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

def bivariate_gaussian(x, y, sigmax = 36, sigmay = 36, mux = 0.0, muy = 0.0, norm = 62):
    # norm = 62
    gauss = norm * np.exp(-((x - mux) ** 2 / (2.0 * sigmax ** 2) + (y - muy) ** 2 / (2.0 * sigmay ** 2)))
    return max(gauss - 12, 0)

def draw_mountains(x, y):
    left_mountains = max(bivariate_gaussian(140. + 0.65 * (x - 140.), y, 32., 32., 140., 140.), bivariate_gaussian(x, y, 30, 30, 80, 248))
    # right_mountains = max(bivariate_gaussian(x, y, 32, 32, 148, 148), bivariate_gaussian(x, y, 32, 32, -20, 148))
    right_mountains = bivariate_gaussian(x, y, 28, 28, 100, 200, 50)
    return max(left_mountains, right_mountains)

def cal_z(x, y):
    return bivariate_gaussian(140. + 0.65 * (x - 140.), y, 32., 32., 140., 140.)

def sol(jpath, obs_jpath, pts_jpath):

    #X, Y = np.mgrid[0:6*np.pi:0.25, 0:4*np.pi:0.25]
    #Z = np.sqrt(np.abs(np.cos(X) + np.cos(Y)))
    #print(X)
    #print(len(X))
    #print(X[0])
    #print(Y)
    #print(Z)
    # fig = plt.figure()
    # ax = fig.add_subplot(projection = '3d')

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    
    with open(obs_jpath, 'r') as f:
        xx = [-1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1]
        yy = [1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1]
        zz = [1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1]
        ojf = json.load(f)
        obs_x = []
        obs_y = []
        obs_z = []
        for i in range(len(ojf)):
            obs = ojf[i]
            xr = obs['len'][0] / 2.
            yr = obs['len'][1] / 2.
            zr = obs['len'][2] / 2.
            for j in range(len(xx)):
                obs_x.append(obs['pt'][0] + xr * xx[j])
                obs_y.append(obs['pt'][1] + yr * yy[j])
                obs_z.append(obs['pt'][2] + zr * zz[j])
            ax.plot(obs_x, obs_y, obs_z, 'b', label = "obstacles")
            obs_x.clear()
            obs_y.clear()
            obs_z.clear()
            for j in range(len(xx)):
                obs_x.append(obs['pt'][0] - xr * xx[j])
                obs_y.append(obs['pt'][1] - yr * yy[j])
                obs_z.append(obs['pt'][2] - zr * zz[j])
            ax.plot(obs_x, obs_y, obs_z, 'b')
            obs_x.clear()
            obs_y.clear()
            obs_z.clear()
    
    with open(pts_jpath, 'r') as f:
        jf = json.load(f)
        pts = jf['waypt']
        # pts = []
        pt_x = []
        pt_y = []
        pt_z = []
        for i in range(len(pts)):
            pt_x.append(pts[i][0])
            pt_y.append(pts[i][1])
            pt_z.append(pts[i][2])
        ax.plot(pt_x, pt_y, pt_z, 'r')


    with open(jpath, 'r') as f:    
        jf = json.load(f)
        traj = jf['traj']
        # traj = []
        traj_x = []
        traj_y = []
        traj_z = []
        for i in range(len(traj)):
            traj_x.append(traj[i][0])
            traj_y.append(traj[i][1])
            traj_z.append(traj[i][2])

        obs = []
        if 'obs' in jf:
            obs = jf['obs']
        obs = []

        obs_x = []
        obs_y = []
        obs_z = []
        bias = 1
        
        xx = [-1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1]
        yy = [1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1]
        zz = [1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1]
        xr = 20
        yr = 1
        zr = 16
        for i in range(len(obs)):
            for j in range(len(xx)):
                obs_x.append(obs[i][0] + bias * xr * xx[j])
                obs_y.append(obs[i][1] + bias * yr * yy[j])
                obs_z.append(obs[i][2] + bias * zr * zz[j])
            ax.plot(obs_x, obs_y, obs_z, 'b', label = "obstacles")
            obs_x.clear()
            obs_y.clear()
            obs_z.clear()
            for j in range(len(xx)):
                obs_x.append(obs[i][0] - bias * xr * xx[j])
                obs_y.append(obs[i][1] - bias * yr * yy[j])
                obs_z.append(obs[i][2] - bias * zr * zz[j])
            ax.plot(obs_x, obs_y, obs_z, 'b')
            obs_x.clear()
            obs_y.clear()
            obs_z.clear()
        end_pt = jf['end_point']
        # end_pt = []
        end_pt_x = []
        end_pt_y = []
        end_pt_z = []
        

        for i in range(len(end_pt)):
            end_pt_x.append(end_pt[i][0])
            end_pt_y.append(end_pt[i][1])
            end_pt_z.append(end_pt[i][2])

        sp_x = [0, 128]
        sp_y = [0, 128]
        sp_z = [0, 128]
        sp = []
        # sp = jf["sp"]

        for i in range(len(sp)):
            sp_x.append(sp[i][0])
            sp_y.append(sp[i][1])
            sp_z.append(sp[i][2])

        # spp_x = [0, 128]
        # spp_y = [0, 128]
        # spp_z = [0, 128]
        spp_x = []
        spp_y = []
        spp_z = []
        spp = []
        # spp = jf["spp"]
        

        for i in range(len(spp)):
            spp_x.append(spp[i][0])
            spp_y.append(spp[i][1])
            spp_z.append(spp[i][2])


    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    X, Y, Z = axes3d.get_test_data(0.05)

    # ax.contour(X, Y, Z, cmap=cm.coolwarm)
    # ax.plot_surface(X, Y, Z, edgecolor='royalblue', lw=0.5, rstride=8, cstride=8, alpha=0.3)
    
    # ax.plot(obs_x, obs_y, obs_z, 'b', label = "obstacles")
    ax.plot(end_pt_x, end_pt_y, end_pt_z, 'ko', markersize = 4, label = "end_point")
    ax.plot(sp_x, sp_y, sp_z, 'yo', markersize = 2)
    ax.plot(spp_x, spp_y, spp_z, 'go', markersize = 2)
    
    # XX = np.arange(-36., 158., 0.5)
    # YY = np.arange(-36., 158., 0.5)
    # print(XX)
    # print(YY)
    # XX, YY = np.meshgrid(XX, YY)
    # print(XX)
    # print(YY)
    # print("check")
    # print(bivariate_gaussian(0., 0.))
    # print(type(XX))
    # print(type(XX[0]))
    # print(type(XX[0][0]))
    # ZZ = draw_mountains(XX, YY)
    # print(XX.ndim)
    # print(type(XX))
    # print(type(ZZ))
    # print(ZZ)
    # ax.plot_surface(XX, YY, ZZ, vmin = ZZ.min() * 2, cmap=cm.Blues)
    # 

    x_range = []
    y_range = []
    for i in range (50):
        x_range.append(-10 + i * 6.4)
        y_range.append(-10 + i * 6.4)
    x_mesh = []
    y_mesh = []
    for i in range(len(x_range)):
        tmpx = []
        for j in range(len(y_range)):
            tmpx.append(x_range[i])
        x_mesh.append(tmpx)

    for i in range(len(y_range)):
        tmpy = []
        for j in range(len(x_range)):
            tmpy.append(y_range[i])
        y_mesh.append(tmpy)

    # print(x_mesh)

    for i in range(len(x_range)):
        zc = []
        for j in range(len(y_range)):
            zc.append(draw_mountains(x_mesh[i][j], y_range[j]))
            # zc.append(bivariate_gaussian(40. + 0.65 * (x_mesh[i][j] - 40.), y_range[j], 32., 32., 40., 40.))
        y_t = []
        for j in range(len(y_range)):
            y_t.append(30 + 1.5 * (y_range[j] - 30))
        # ax.plot(x_mesh[i], y_t, zc, 'g')
        ax.plot(x_mesh[i], y_range, zc, 'g')

    for i in range(len(y_range)):
        zc = []
        for j in range(len(x_range)):
            # zc.append(bivariate_gaussian(40. + 0.65 * (x_range[j] - 40.), y_mesh[i][j], 32., 32., 40., 40.))
            zc.append(draw_mountains(x_range[j], y_mesh[i][j]))
        y_t = []
        for j in range(len(x_range)):
            y_t.append(30 + 1.5 * (y_mesh[i][j] - 30))
        # ax.plot(x_range, y_t, zc, 'g')
        ax.plot(x_range, y_mesh[i], zc, 'g')

    ax.plot(traj_x, traj_y, traj_z, 'y', label = "search result")

    plt.legend()
    plt.show()

if __name__ == '__main__':
    # print(bivariate_gaussian(1., 1., 32., 32., 36., 36))
    # print(cal_z(42.8, 8.24))
    # print(cal_z(27.84, 6.72))
    # print(cal_z(1.6, 1.6))
    # print(cal_z(23.36, 10.56))
    print(cal_z(77.28, 116.64))
    print(cal_z(69.12, 115.2))
    jpath = "/home/alan/下载/global_path_search/build/astar_result.json"
    obs_jpath = "/home/alan/下载/global_path_search/build/obs.json"
    pts_jpath = "/home/alan/下载/global_path_search/build/final_waypts.json"
    sol(jpath, obs_jpath, pts_jpath)
