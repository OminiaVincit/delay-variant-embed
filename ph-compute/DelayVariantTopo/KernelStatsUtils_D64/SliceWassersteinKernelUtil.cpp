#include "stdafx.h"
#include "SliceWassersteinKernelUtil.h"
#include <numeric>
#include <ppl.h>
#include <ppltasks.h>

namespace NSliceWassersteinKernelUtils {
    template <class Dtype>
    std::vector<Dtype> PersistentVecThetaProducts(const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes,
        const Dtype theta) {
        std::vector<Dtype> prod_vec;
        auto stheta = std::sin(theta);
        auto ctheta = std::cos(theta);

        for (auto bar : barcodes) {
            const Dtype tmp = std::get<0>(bar) * ctheta + std::get<1>(bar) * stheta;
            prod_vec.push_back(tmp);
        }
        std::sort(prod_vec.begin(), prod_vec.end());
        return prod_vec;
    }

    template <class Dtype>
    std::vector<Dtype> PersistentVecThetaProducts(const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes,
        const Dtype theta, const Dtype phi) {
        std::vector<Dtype> prod_vec;
        auto stheta = std::sin(theta);
        auto x = stheta * std::cos(phi);
        auto y = stheta * std::sin(phi);
        auto z = std::cos(theta);

        for (auto bar : barcodes) {
            const Dtype tmp = std::get<0>(bar) * x + std::get<1>(bar) * y + std::get<2>(bar) * z;
            prod_vec.push_back(tmp);
        }
        std::sort(prod_vec.begin(), prod_vec.end());
        return prod_vec;
    }

    template <class Dtype>
    bool ProjectToDiagonalAndAppend(std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2) {
        std::vector<std::tuple<Dtype, Dtype, Dtype>> proj_1to2;
        for (auto bar : barcodes_1) {
            auto tmp = (std::get<0>(bar) + std::get<1>(bar)) / 2.0;
            proj_1to2.push_back(std::make_tuple(tmp, tmp, std::get<2>(bar)));
        }

        for (auto bar : barcodes_2) {
            auto tmp = (std::get<0>(bar) + std::get<1>(bar)) / 2.0;
            barcodes_1.push_back(std::make_tuple(tmp, tmp, std::get<2>(bar)));
        }
        std::copy(proj_1to2.begin(), proj_1to2.end(), std::back_inserter(barcodes_2));
        return true;
    }

    template <class Dtype>
    KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, size_t M /*number of direction*/) {
        if (M <= 0)
            return false;
        std::vector<std::tuple<Dtype, Dtype, Dtype>> tbars_1 = barcodes_1;
        std::vector<std::tuple<Dtype, Dtype, Dtype>> tbars_2 = barcodes_2;
        
        ProjectToDiagonalAndAppend(tbars_1, tbars_2);
        
        Dtype sw = 0;
        Dtype theta = -M_PI_2;
        Dtype s = M_PI / M;
        for (int i = 0; i <= M; ++i) {
            auto p1vec = PersistentVecThetaProducts(tbars_1, theta);
            auto p2vec = PersistentVecThetaProducts(tbars_2, theta);
            if (p1vec.size() != p2vec.size()) {
                return false;
            }
            size_t n = p1vec.size();
            std::vector<Dtype> values(n, (Dtype)0.0);
            concurrency::parallel_for((size_t)0, n, [&](size_t idx){
                values[idx] = abs(p1vec[idx] - p2vec[idx]);
            });
            const Dtype local_sum = std::accumulate(values.begin(), values.end(), 0.0, std::plus<Dtype>());
            sw += s * local_sum;
            theta += s;
        }
        result = sw * M_1_PI;

        return true;
    }

    template <class Dtype>
    KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2,
        size_t theta_M /*number of direction for theta*/,
        size_t phi_M /*number of direction for phi*/) {
        if (theta_M <= 0 || phi_M <= 0)
            return false;

        std::vector<std::tuple<Dtype, Dtype, Dtype>> tbars_1 = barcodes_1;
        std::vector<std::tuple<Dtype, Dtype, Dtype>> tbars_2 = barcodes_2;

        ProjectToDiagonalAndAppend(tbars_1, tbars_2);
        
        Dtype sw = 0;
        Dtype phi_0 = 0, theta_0 = -M_PI_2;
        Dtype s_phi = M_PI_2 / phi_M;
        Dtype s_theta = M_PI / theta_M;
        std::vector<Dtype> phi_svalues(phi_M + 1, (Dtype)0.0);

        concurrency::parallel_for((size_t)0, phi_M + 1, [&](size_t phi_idx) {
            Dtype phi = phi_0 + phi_idx * s_phi;
            std::vector<Dtype> theta_svalues(theta_M + 1, (Dtype)0.0);
            concurrency::parallel_for((size_t)0, theta_M + 1, [&](size_t theta_idx) {
                Dtype theta = theta_0 + theta_idx * s_theta;
                auto p1vec = PersistentVecThetaProducts(tbars_1, theta, phi);
                auto p2vec = PersistentVecThetaProducts(tbars_2, theta, phi);
                if (p1vec.size() == p2vec.size()) {
                    size_t n = p1vec.size();
                    std::vector<Dtype> values(n, (Dtype)0.0);
                    concurrency::parallel_for((size_t)0, n, [&](size_t idx) {
                        values[idx] = abs(p1vec[idx] - p2vec[idx]);
                    });
                    theta_svalues[theta_idx] = std::accumulate(values.begin(), values.end(), 0.0, std::plus<Dtype>());
                }
            });
            phi_svalues[phi_idx] = M_1_PI * std::accumulate(theta_svalues.begin(), theta_svalues.end(), 0.0, std::plus<Dtype>());
        });
        result = M_2_PI * std::accumulate(phi_svalues.begin(), phi_svalues.end(), 0.0, std::plus<Dtype>());;

        return true;
    }
}