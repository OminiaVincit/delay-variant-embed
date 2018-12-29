#pragma once

#include "stdafx.h"
#include "KernelStatsUtils_D64.h"

namespace NRiemannianManifoldKernelUtils {
    template <class Dtype>
    KernelStatsUtils_D64_API bool RiemannGeodesicMetric(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype T1, Dtype T2);

    KernelStatsUtils_D64_API bool RiemannGeodesicMetric(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2, float T1, float T2);
}