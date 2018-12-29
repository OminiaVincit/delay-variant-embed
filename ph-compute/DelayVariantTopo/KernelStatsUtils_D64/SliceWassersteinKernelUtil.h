#pragma once
#include "stdafx.h"
#include "KernelStatsUtils_D64.h"

namespace NSliceWassersteinKernelUtils {
    template <class Dtype>
    KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(Dtype& result,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, size_t M /*number of direction*/);

    template KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(float& result,
        const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2, size_t M /*number of direction*/);

    template KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(double& result,
        const std::vector<std::tuple<double, double, double>>& barcodes_1,
        const std::vector<std::tuple<double, double, double>>& barcodes_2, size_t M /*number of direction*/);

    template <class Dtype>
    KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(Dtype& result,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2,
        size_t theta_M /*number of direction for theta*/,
        size_t phi_M /*number of direction for phi*/);

    template KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(float& result,
        const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2,
        size_t theta_M /*number of direction for theta*/,
        size_t phi_M /*number of direction for phi*/);

    template KernelStatsUtils_D64_API bool SliceCircleWassersteinDistance(double& result,
        const std::vector<std::tuple<double, double, double>>& barcodes_1,
        const std::vector<std::tuple<double, double, double>>& barcodes_2,
        size_t theta_M /*number of direction for theta*/,
        size_t phi_M /*number of direction for phi*/);
}