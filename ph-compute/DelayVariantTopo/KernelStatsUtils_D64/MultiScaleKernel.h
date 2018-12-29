#pragma once
#include "stdafx.h"
#include "KernelStatsUtils_D64.h"

namespace NMultiScaleKernelUtils {
    template <class Dtype>
    KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype time);

    template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2, float time);

    template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(double& result, const std::vector<std::tuple<double, double, double>>& barcodes_1,
        const std::vector<std::tuple<double, double, double>>& barcodes_2, double time);

	template <class Dtype>
	KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
		const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype T1, Dtype T2);

	template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
		const std::vector<std::tuple<float, float, float>>& barcodes_2, float T1, float T2);

	template KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(double& result, const std::vector<std::tuple<double, double, double>>& barcodes_1,
		const std::vector<std::tuple<double, double, double>>& barcodes_2, double T1, double T2);

    template <class Dtype>
    KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype time);

    template KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2, float time);

    template KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(double& result, const std::vector<std::tuple<double, double, double>>& barcodes_1,
        const std::vector<std::tuple<double, double, double>>& barcodes_2, double time);

	KernelStatsUtils_D64_API bool L2DiagramMskInnerProductSSE(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
		const std::vector<std::tuple<float, float, float>>& barcodes_2, float time);

    KernelStatsUtils_D64_API bool L2DiagramMskInnerProductSSE(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2, float T1, float T2);
}