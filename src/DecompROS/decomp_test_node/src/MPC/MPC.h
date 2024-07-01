#pragma once

#include "BlockMatrix.h"
#include "FunctionG.h"
#include <Eigen/Dense>
#include "json.hpp"
#include "cmath"

template <typename T, int HorizonNum, int SizeX, int SizeU, int SizeYx, int SizeYu, int SizeEqx, int SizeEqu, int NumEllx, int NumEllu, int SizeG, int SizeEllx[], int SizeEllu[]>
class MPC_ADMMSolver{
private:
    typedef Eigen::Matrix<T, SizeX, SizeX> MatrixX;
    typedef Eigen::Matrix<T, SizeU, SizeU> MatrixU;
    typedef Eigen::Matrix<T, SizeX, SizeU> MatrixB;
    typedef Eigen::Matrix<T, SizeX, 1> VectorX;
    typedef Eigen::Matrix<T, SizeU, 1> VectorU;
    typedef Eigen::Matrix<T, SizeYx, SizeX> MatrixHx;
    typedef Eigen::Matrix<T, SizeYu, SizeU> MatrixHu;
    typedef BlockVector<T, HorizonNum + 1, SizeX + SizeU> BlockVectorComb;
    typedef BlockVector<T, HorizonNum + 1, SizeX> BlockVectorX;
    typedef BlockVector<T, HorizonNum + 1, SizeU> BlockVectorU;
    typedef BlockVector<T, HorizonNum + 1, SizeYx + SizeYu> BlockVectorW;
    const float inf = 1e6;

public:
    /**
     * MPC problem coefficients
     */
    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    VectorX x_init;
    std::array<MatrixHx, HorizonNum + 1> Hx;
    std::array<MatrixHu, HorizonNum + 1> Hu;
    std::array<Eigen::Matrix<T, SizeYx - SizeEqx, 1>, HorizonNum + 1> new_centerX;
    std::array<Eigen::Matrix<T, SizeYu - SizeEqu, 1>, HorizonNum + 1> new_centerU;
    std::array<std::array<FunctionG<T>, SizeG>, HorizonNum + 1> g;
    //const std::vector<int> SizeEllx;
    //const std::vector<int> SizeEllu;
    int K;
    int KAcc;
    T costweight;

    /**
     * ADMM coefficients
     */
    BlockLowerBidiagonalMatrix<T, HorizonNum + 1, SizeX, SizeX + SizeU> bA;
    BlockLowerBidiagonalMatrix<T, HorizonNum + 1, SizeX, SizeX> bL;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeYx + SizeYu, SizeX + SizeU> bE;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeX + SizeU, SizeYx + SizeYu> bET;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeX + SizeU, SizeX + SizeU> bG, biG;
    BlockSymmetricTridiagonalMatrix<T, HorizonNum + 1, SizeX> bAiGAT;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeX + SizeU, SizeX + SizeU> barE;
    BlockVectorX bb;
    BlockVectorComb bF;
    T rho = 1;

    /**
     * ADMM new coefficients
     */
    std::array<MatrixX, HorizonNum> barA;
    std::array<MatrixB, HorizonNum> barB;
    std::array<Eigen::Matrix<T, SizeU, SizeX>, HorizonNum> barBT;
    std::array<VectorX, HorizonNum> barc;
    std::array<MatrixHx, HorizonNum + 1> bEx;
    std::array<MatrixHu, HorizonNum + 1> bEu;
    /**
     * LQR backward matrix
     */
    std::array<MatrixX, HorizonNum + 1> Q_lqr;
    std::array<MatrixU, HorizonNum + 1> R_lqr;
    std::array<MatrixU, HorizonNum + 1> G_lqr;
    std::array<MatrixB, HorizonNum + 1> K_lqr;
    std::array<MatrixX, HorizonNum + 1> Fx_lqr;
    std::array<MatrixX, HorizonNum + 1> P_lqr;
    std::array<VectorU, HorizonNum + 1> BTPc;
    std::array<MatrixX, HorizonNum + 1> FTP;
    std::array<MatrixB, HorizonNum + 1> KRT;
    VectorX bar_x_init;
    MPC_ADMMSolver(
        std::array<MatrixX, HorizonNum + 1> _Q,
        std::array<MatrixU, HorizonNum + 1> _R,
        std::array<VectorX, HorizonNum + 1> _L,
        std::array<VectorU, HorizonNum + 1> _W,
        std::array<MatrixX, HorizonNum> _A,
        std::array<MatrixB, HorizonNum> _B,
        std::array<VectorX, HorizonNum> _c,
        VectorX _x_init,
        std::array<MatrixHx, HorizonNum + 1> _Hx,
        std::array<MatrixHu, HorizonNum + 1> _Hu,
        std::array<Eigen::Matrix<T, SizeYx - SizeEqx, 1>, HorizonNum + 1> _new_centerX,
        std::array<Eigen::Matrix<T, SizeYu - SizeEqu, 1>, HorizonNum + 1> _new_centerU,
        std::array<std::array<FunctionG<T>, SizeG>, HorizonNum + 1> _g,
        //const std::vector<int> _SizeEllx,
        //const std::vector<int> _SizeEllu,
        T _costweight,
        int _KAcc,
        int _K = 250
    ) : Q(_Q), R(_R), L(_L), W(_W), A(_A), B(_B), c(_c), x_init(_x_init), Hx(_Hx), Hu(_Hu), new_centerX(_new_centerX), new_centerU(_new_centerU), g(_g), costweight(_costweight), KAcc(_KAcc), K(_K) {}

    void PreScaling1() {
        // P_x, P_u
        std::array<MatrixX, HorizonNum + 1> Px, iPx, PxT, iPxT;
        std::array<MatrixU, HorizonNum + 1> Pu, iPu, PuT, iPuT;
        for(int i = 0; i <= HorizonNum; ++i) {
            Px[i] = Q[i].llt().matrixL();
            Pu[i] = R[i].llt().matrixL();
            iPx[i] = Px[i].inverse();
            iPu[i] = Pu[i].inverse();
            PxT[i] = Px[i].transpose();
            PuT[i] = Pu[i].transpose();
            iPxT[i] = PxT[i].inverse();
            iPuT[i] = PuT[i].inverse();
        }
        bar_x_init = PxT[0] * x_init;
        
        std::array<Eigen::Matrix<T, SizeYx, SizeYx>, HorizonNum + 1> LambX;
        std::array<Eigen::Matrix<T, SizeYu, SizeYu>, HorizonNum + 1> LambU;
        Eigen::Matrix<T, SizeEqx, SizeEqx> LambXi;
        Eigen::Matrix<T, SizeEqu, SizeEqu> LambUi;
    
        for(int i = 0; i <= HorizonNum; ++i) { //这里原来是<=,但是程序会崩溃，取不到49
        
            // std::cout << "loop: " << i << std::endl;
        
            // \bar{A}
            if (i < HorizonNum) {
                barA[i] = PxT[i + 1] * A[i] * iPxT[i];
            }
        
            // std::cout << "here1" << std::endl;
            // \bar{b}
            if (i < HorizonNum) {
                barB[i] = PxT[i + 1] * B[i] * iPuT[i];
                barBT[i] = barB[i].transpose();
            }
            // \bar{c}
            if (i < HorizonNum) {
                barc[i] = PxT[i + 1] * c[i];
            }
            // \bar{F}
            bF.v[i] << iPx[i] * L[i], iPu[i] *  W[i];
           
            // \Lambda_x, \Lambda_u
            LambXi = (Hx[i].block(0, 0, SizeEqx, SizeX) * iPxT[i]).rowwise().norm().asDiagonal().inverse();
            for (int i = 0; i < SizeEqx; i++) {
                if (LambXi(i, i) > inf) {
                    LambXi(i, i) = 1;
                }
            }
            LambX[i] << LambXi, Eigen::Matrix<T, SizeEqx, SizeYx-SizeEqx>::Zero(SizeEqx, SizeYx - SizeEqx),
                        Eigen::Matrix<T, SizeYx - SizeEqx, SizeEqx>::Zero(SizeYx - SizeEqx, SizeEqx), Eigen::Matrix<T, SizeYx - SizeEqx, SizeYx - SizeEqx>::Identity(SizeYx - SizeEqx, SizeYx - SizeEqx);
        
            LambUi = (Hu[i].block(0, 0, SizeEqu, SizeU) * iPuT[i]).rowwise().norm().asDiagonal().inverse();
            for (int i = 0; i < SizeEqu; i++) {
                if (LambUi(i, i) > inf) {
                    LambUi(i, i) = 1;
                }
            }
            LambU[i] << LambUi, Eigen::Matrix<T, SizeEqu, SizeYu - SizeEqu>::Zero(SizeEqu, SizeYu - SizeEqu),
                        Eigen::Matrix<T, SizeYu - SizeEqu, SizeEqu>::Zero(SizeYu - SizeEqu, SizeEqu), Eigen::Matrix<T, SizeYu - SizeEqu, SizeYu - SizeEqu>::Identity(SizeYu - SizeEqu, SizeYu - SizeEqu);
            
            for(int j = 0; j < SizeEqx; ++j) g[i][j].Prescaling(LambX[i].diagonal()(j));
            for(int j = 0; j < SizeEqu; ++j) g[i][j + SizeEqx + NumEllx].Prescaling(LambU[i].diagonal()(j));
            // (x,u) = barE * barz; w = bE * barz;
            
            barE.d[i] << iPxT[i], Eigen::Matrix<T, SizeX, SizeU>::Zero(SizeX, SizeU),
                      Eigen::Matrix<T, SizeU, SizeX>::Zero(SizeU, SizeX), iPuT[i];
            bE.d[i] << LambX[i] * Hx[i] * iPxT[i], Eigen::Matrix<T, SizeYx, SizeU>::Zero(SizeYx, SizeU),
                      Eigen::Matrix<T, SizeYu, SizeX>::Zero(SizeYu, SizeX), LambU[i] * Hu[i] * iPuT[i];
            bET.d[i] = bE.d[i].transpose();
            // std::cout << "here7" << std::endl;
            bEx[i] = LambX[i] * Hx[i] * iPxT[i];
            bEu[i] = LambU[i] * Hu[i] * iPuT[i];
            /**/
            
            // std::cout << "loop: " << i << std::endl;
        }
        
        // std::cout << "func finish: " << std::endl;
    }


    void ADMMPrework1() {
        Eigen::Matrix<T, SizeX, SizeX> Ix;
        Ix.setIdentity(SizeX, SizeX);
        Eigen::Matrix<T, SizeU, SizeU> Iu;
        Iu.setIdentity(SizeU, SizeU);

        for (int i = HorizonNum; i >= 0; --i) {
            Q_lqr[i] = (rho * bEx[i].transpose() * bEx[i] + Ix) * 0.5;
            R_lqr[i] = (rho * bEu[i].transpose() * bEu[i] + Iu) * 0.5;
        }
        P_lqr[HorizonNum] = Q_lqr[HorizonNum];
        K_lqr[HorizonNum] = MatrixB::Zero(SizeX, SizeU);
        for (int i = HorizonNum - 1; i >= 0; --i) {
            G_lqr[i] = (barBT[i] * P_lqr[i + 1] * barB[i] + R_lqr[i]).inverse();
            K_lqr[i] = -barA[i].transpose() * P_lqr[i + 1] * barB[i] * G_lqr[i].transpose();
            Fx_lqr[i] = barA[i] + barB[i] * K_lqr[i].transpose();
            P_lqr[i] = Fx_lqr[i].transpose() * P_lqr[i + 1] * Fx_lqr[i] + K_lqr[i] * R_lqr[i] * K_lqr[i].transpose() + Q_lqr[i];
            BTPc[i] = -2 * barBT[i] * P_lqr[i + 1] * barc[i];
            FTP[i] = Fx_lqr[i].transpose() * P_lqr[i + 1];
            KRT[i] = K_lqr[i] * R_lqr[i].transpose();
        }
    }

    BlockVectorComb LQR_Solver1(BlockVectorComb f) {
        std::array<VectorU, HorizonNum + 1> r_lqr;
        std::array<VectorX, HorizonNum + 1> E_lqr;
        std::array<VectorX, HorizonNum + 1> Fr_lqr;
        for (auto &mat : r_lqr) {
            mat.setZero();
        }
        for (auto &mat : E_lqr) {
            mat.setZero();
        }
        for (auto &mat : Fr_lqr) {
            mat.setZero();
        }

        E_lqr[HorizonNum] = f.v[HorizonNum].block(0, 0, SizeX, 1);
        for (int i = HorizonNum - 1; i >= 0; --i) {
            // i = 40的时候E_lqr里有nan
            r_lqr[i] = 0.5 * G_lqr[i] * (BTPc[i] - barBT[i] * E_lqr[i + 1] - f.v[i].block(SizeX, 0, SizeU, 1));
            Fr_lqr[i] = barB[i] * r_lqr[i];
            E_lqr[i] = 2 * FTP[i] * (Fr_lqr[i] + barc[i]) + Fx_lqr[i].transpose() * E_lqr[i + 1] + 2 * KRT[i] * r_lqr[i] + f.v[i].block(0, 0, SizeX, 1) + K_lqr[i] * f.v[i].block(SizeX, 0, SizeU, 1);
        }
        BlockVectorComb result;
        result.setZero();
        result.v[0].block(0, 0, SizeX, 1) = bar_x_init;
        result.v[HorizonNum].block(SizeX, 0, SizeU, 1) = -0.5 * R_lqr[HorizonNum].inverse() * f.v[HorizonNum].block(SizeX, 0, SizeU, 1);
        for (int i = 0; i < HorizonNum; ++i) {
            result.v[i + 1].block(0, 0, SizeX, 1) = Fx_lqr[i] * result.v[i].block(0, 0, SizeX, 1) + Fr_lqr[i] + barc[i]; // 第二轮Fr_lqr里全是nan
            result.v[i].block(SizeX, 0, SizeU, 1) = K_lqr[i].transpose() * result.v[i].block(0, 0, SizeX, 1) + r_lqr[i]; // 第二轮r_lqr里全是nan
        }
        return result;
    }
    
    BlockVectorW G_Solver(BlockVectorComb barz, BlockVectorW nu, T rho) {
        BlockVectorW u = bE * barz - nu;
        for(int i = 0; i <= HorizonNum; ++i) {
            for(int j = 0; j < SizeEqx; ++j) {
                T ui = u.v[i](j); 
                u.v[i](j) = g[i][j].Minimizer(ui, rho);
            }
            int cntx = 0;
            for (int k = 0; k < NumEllx; k++) {
                Eigen::Matrix<T, Eigen::Dynamic, 1> d = u.v[i].block(SizeEqx + cntx, 0, SizeEllx[k], 1);
                Eigen::Matrix<T, Eigen::Dynamic, 1> new_center = new_centerX[i].block(cntx, 0, SizeEllx[k], 1);
                T L = std::sqrt((new_center - d).transpose() * (new_center - d));
                T res = g[i][SizeEqx + k].Minimizer(L, rho);
                u.v[i].block(SizeEqx + cntx, 0, SizeEllx[k], 1) = res / std::max(L, 1e-5f) * (d - new_center) + new_center;
                cntx += SizeEllx[k];
            }
            for (int j = 0; j < SizeEqu; j++) {
                T ui = u.v[i](SizeYx + j); 
                u.v[i](SizeYx + j) = g[i][SizeEqx + NumEllx + j].Minimizer(ui, rho);
            }
            int cntu = 0;
            for (int k = 0; k < NumEllu; k++) {
                Eigen::Matrix<T, Eigen::Dynamic, 1> d = u.v[i].block(SizeYx + SizeEqu + cntu, 0, SizeEllu[k], 1);// (0,0)^T
                Eigen::Matrix<T, Eigen::Dynamic, 1> new_center = new_centerU[i].block(cntu, 0, SizeEllu[k], 1);// (0,0)^T
                T L = std::sqrt((new_center - d).transpose() * (new_center - d));// 0
                T res = g[i][SizeEqx + NumEllx + SizeEqu + k].Minimizer(L, rho);
                u.v[i].block(SizeYx + SizeEqu + cntu, 0, SizeEllu[k], 1) = res / std::max(L, 1e-5f) * (d - new_center) + new_center;// 除0了
                cntx += SizeEllu[k];
            }
        }
        return u;
    }

    T CalculateCost(BlockVector<T, HorizonNum + 1, SizeX + SizeU> barz) {
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> z = barE * barz;
        // zhenzhengde pvmu
        BlockVectorW w = bE * barz;
        T res1 = 0, res2 = 0, res3 = 0;
        std::array<std::array<T, SizeG>, HorizonNum + 1> d;
        for(int i = 0; i <= HorizonNum; ++i) {
            VectorX x = z.v[i].block(0, 0, SizeX, 1);
            VectorU u1 = z.v[i].block(SizeX, 0, SizeU, 1);
            res1 += 0.5 * x.transpose() * (Q[i] * x) + L[i].dot(x);
            VectorU v = R[i] * u1;
            res1 += 0.5 * u1.transpose() * v + W[i].dot(u1);
            VectorX Li = Q[i].inverse() * L[i];
            VectorU Wi = R[i].inverse() * W[i];
            res1 += 0.5 * L[i].transpose() * Li;
            res1 += 0.5 * W[i].transpose() * Wi;
            //cout<< i<< ' '<< u.transpose() << ' ' << 0.5 * u.transpose() * v + W[i].dot(u) << ' ' << 0.5 * x.transpose() * (Q[i] * x) <<' '<< L[i].dot(x)<<' '<<0.5 * (L[i][0] * L[i][0] / 1 + L[i][1] * L[i][1] / 1)<<' '<<x[0]<<' '<<x[1]<<' '<<x[2]<<endl;
            std::array<T, SizeG> wi;
            for(int j = 0; j < SizeEqx; ++j) {
                wi[j] = w.v[i](j);
            }
            int cntx = 0;
            for (int k = 0; k < NumEllx; k++) {
                Eigen::Matrix<T, Eigen::Dynamic, 1> d = w.v[i].block(SizeEqx + cntx, 0, SizeEllx[k], 1);
                Eigen::Matrix<T, Eigen::Dynamic, 1> new_center = new_centerX[i].block(cntx, 0, SizeEllx[k], 1);
                wi[k + SizeEqx] = std::sqrt((new_center - d).transpose() * (new_center - d));
                cntx += SizeEllx[k];
            }
            for (int j = 0; j < SizeEqu; j++) {
                wi[j + SizeEqx + NumEllx] = w.v[i](j + SizeYx);
            }
            int cntu = 0;
            for (int k = 0; k < NumEllu; k++) {
                Eigen::Matrix<T, Eigen::Dynamic, 1> d = w.v[i].block(SizeYx + SizeEqu + cntu, 0, SizeEllu[k], 1);
                Eigen::Matrix<T, Eigen::Dynamic, 1> new_center = new_centerU[i].block(cntu, 0, SizeEllu[k], 1);
                wi[k + SizeEqx + NumEllx + SizeEqu] = std::sqrt((new_center - d).transpose() * (new_center - d));
                cntx += SizeEllu[k];
            }
            for (int j = 0; j < SizeG; ++j) {
                res2 += g[i][j].CostOfQuadraticPart(wi[j]);
                d[j][i] = g[i][j].DistanceOfIndicatorPart(wi[j]);
                // std::cout << " step:" << i << "g{" << j <<"}" << "input:" << wi[j]<< " ,g value: " << g[i][j].CostOfQuadraticPart(wi[j]) << std::endl;
            }

        }
        for(int i = 0; i < HorizonNum; ++i) {
            for (int j = 0; j < SizeG; ++j) {
                res3 += costweight * (d[j][i + 1] - d[j][i] > 0 ? d[j][i + 1] - d[j][i] : 0);
            }
        }
        // std::cout << "funcG: " << g[0][11].CostOfQuadraticPart(0.5);
        return res1 + res2 + res3;
    }

    BlockVectorComb ADMMIteration() {
        int k = 1; // iteration num
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> barz, res;
        BlockVector<T, HorizonNum + 1, SizeYx + SizeYu> w, nu;
        barz.setZero(); w.setZero(); nu.setZero(); 
        res = barz;
        T min_cost = inf;
        std::vector<T> cvec;
        while(k <= K) {
            // barz.print();
            BlockVectorComb f = bF - (bET * (w + nu)) * (rho);
            // f.print();
            barz = LQR_Solver1(f);
            // barz.print();
            w = G_Solver(barz, nu, rho);
            nu = nu + w - bE * barz;
            k++;
            T cost = CalculateCost(barz);
            cvec.push_back(cost);
            if(!(k % 10)) {
                if(cost < min_cost) {
                    res = barz;
                    min_cost = cost;
                }
                // cout << "Episodes: " << k << ' ' << cost << endl;
            }
        }

        using json = nlohmann::json;
        json j;
        json jcost = json::array();
        for(int i = 0; i < K; ++i) {
            jcost.push_back(cvec[i]);
        }
        j["cost"] = jcost;
        ofstream out("testori.out");
        if(out.is_open()) {
            out << j.dump(4) << endl;
            out.close();
        }

        cout << "Final Cost: " << min_cost << endl;
        return barE * res;
    }

    BlockVectorComb ADMMIterationQuick() {
        int k = 1; // iteration num
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> barz, res;
        BlockVector<T, HorizonNum + 1, SizeYx + SizeYu> w, nu, nu_minus, lambda, nu_best, w_best;
        barz.setZero(); w.setZero(); nu.setZero(); lambda.setZero();
        res = barz;
        T min_cost = inf;
        std::vector<T> cvec;
        int k_acc = KAcc < 0.6 * K ? KAcc : 0.6* K;
        // barz = LQR_Solver1(bF - (bET * (w + nu_minus)) * (rho));
        while(k <= K) {
            // if (k == K)
                // barz.print();
            // Nestrov Acceleration
            if (k < k_acc) {
                float step = 1;
                float gamma = 2. / (k + 1);
                //cout<<"test: "<<step / gamma<<' '<<step<<' '<<gamma<<endl;
                nu_minus = nu * (1 - gamma) + lambda * gamma;
                // auto start_time1 = std::chrono::high_resolution_clock::now();

                barz = LQR_Solver1(bF - (bET * (w + nu_minus)) * (rho));
                // auto end_time1 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration1 = end_time1 - start_time1;
                // std::cout << "Time taken by LQR_Solver1: " << duration1.count() << " seconds" << std::endl;

                // auto start_time2 = std::chrono::high_resolution_clock::now();
                w = G_Solver(barz, nu_minus, rho);
                // auto end_time2 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration2 = end_time2 - start_time2;
                // std::cout << "Time taken by G_Solver: " << duration2.count() << " seconds" << std::endl;
                lambda = lambda + (w - bE * barz) * (step / gamma);
                nu = nu * (1 - gamma) + lambda * gamma;
                k++;
                
                // cvec.push_back(cost);
                if(!(k % 10)) {
                    T cost = CalculateCost(barz);
                    if(cost < min_cost) {
                        res = barz;
                        nu_best = nu;
                        w_best = w;
                        min_cost = cost;
                    }
                    // cout << "Episodes: " << k << ' ' << cost << endl;
                }
            } else { // Normal Iteration
                if (k == k_acc) {
                    barz = res;
                    w = w_best;
                    nu = nu_best;
                } else {
                    barz = LQR_Solver1(bF - (bET * (w + nu)) * (rho));
                    w = G_Solver(barz, nu, rho);
                    nu = nu + w - bE * barz;
                }
                k++;
                
                // cvec.push_back(cost);
                if(!(k % 10)) {
                    T cost = CalculateCost(barz);
                    if(cost < min_cost) {
                        res = barz;
                        min_cost = cost;
                    }
                    // cout << "Episodes: " << k << ' ' << cost << endl;
                }
            }
        }
        // using json = nlohmann::json;
        // json j;
        // json jcost = json::array();
        // for(int i = 0; i < K; ++i) {
        //     jcost.push_back(cvec[i]);
        // }
        // j["cost"] = jcost;
        // ofstream out("testacc.out");
        // if(out.is_open()) {
        //     out << j.dump(4) << endl;
        //     out.close();
        // }

        cout << "Final Cost: " << min_cost << endl;
        return barE * res;
    }
    
    BlockVector<T, HorizonNum + 1, SizeX + SizeU> solve() {
        // std::cout << "start PreScaling1" << std::endl;
        auto start_time1 = std::chrono::high_resolution_clock::now();
        PreScaling1();
        auto end_time1 = std::chrono::high_resolution_clock::now();
        // std::cout << "finish PreScaling1" << std::endl;
        auto start_time2 = std::chrono::high_resolution_clock::now();
        ADMMPrework1();
        auto end_time2 = std::chrono::high_resolution_clock::now();
        // std::cout << "start PreScaling2" << std::endl;
        // BlockVector<T, HorizonNum + 1, SizeX + SizeU> res;
        // PreScaling2(res);
        // std::cout << "finish ADMMPrework1" << std::endl;
        auto start_time3 = std::chrono::high_resolution_clock::now();
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> res = ADMMIterationQuick();
        auto end_time3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration1 = end_time1 - start_time1;
        std::chrono::duration<double> duration2 = end_time2 - start_time2;
        std::chrono::duration<double> duration3 = end_time3 - start_time3;
        // std::cout << "finish ADMMIterationQuick" << std::endl;
        std::cout << "Time taken by PreScaling1: " << duration1.count() << " seconds" << std::endl;
        std::cout << "Time taken by ADMMPrework1: " << duration2.count() << " seconds" << std::endl;
        std::cout << "Time taken by ADMMIterationQuick: " << duration3.count() << " seconds" << std::endl;
        return res;
    }
    

    // void solve(BlockVector<T, HorizonNum + 1, SizeX + SizeU>& res) {
    //     // std::cout << "start PreScaling1" << std::endl;
    //     // PreScaling1();
    //     // std::cout << "finish PreScaling1" << std::endl;
    //     // ADMMPrework1();
    //     std::cout << "start PreScaling2" << std::endl;
    //     PreScaling2(res);
    //     std::cout << "finish ADMMPrework1" << std::endl;
    //     // BlockVector<T, HorizonNum + 1, SizeX + SizeU> res = ADMMIterationQuick();
    //     std::cout << "finish ADMMIterationQuick" << std::endl;
    //     return;
    // }


};
