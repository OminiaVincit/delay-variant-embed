#include "stdafx.h"
#include <numeric>
#include "MultiScaleKernel.h"

#include <ppl.h>
#include <ppltasks.h>


namespace NMultiScaleKernelUtils {
    template <class Dtype>
    KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype time)
    {
        std::vector<Dtype> a;
        std::vector<Dtype> mu_1;
        std::vector<Dtype> mu_2;
        std::vector<Dtype> mu_3;

        std::vector<Dtype> b;
        std::vector<Dtype> nu_1;
        std::vector<Dtype> nu_2;
        std::vector<Dtype> nu_3;

        int64_t num_points_1 = (int64_t)barcodes_1.size();
        for (auto tmp : barcodes_1) {

            a.push_back(1.0);
            mu_1.push_back(std::get<0>(tmp));
            mu_2.push_back(std::get<1>(tmp));
            mu_3.push_back(std::get<2>(tmp));

            // mirrored points
            a.push_back(-1.0);
            mu_1.push_back(std::get<1>(tmp));
            mu_2.push_back(std::get<0>(tmp));
            mu_3.push_back(std::get<2>(tmp));
        }

        int64_t num_points_2 = (int64_t)barcodes_2.size();
        for (auto tmp : barcodes_2) {

            b.push_back(1.0);
            nu_1.push_back(std::get<0>(tmp));
            nu_2.push_back(std::get<1>(tmp));
            nu_3.push_back(std::get<2>(tmp));

            // mirrored points
            b.push_back(-1.0);
            nu_1.push_back(std::get<1>(tmp));
            nu_2.push_back(std::get<0>(tmp));
            nu_3.push_back(std::get<2>(tmp));
        }

        size_t local_end = 2 * num_points_1;
        std::vector<Dtype> gathered_values(local_end, 0);
		Dtype thres_zero = 70, thres_one = 1e-10;
        concurrency::parallel_for((size_t)0, local_end, [&](size_t idx_1) {
            Dtype local_sum = 0;
            for (size_t idx_2 = 0; idx_2 < 2 * num_points_2; idx_2++) {
                Dtype diff_3 = mu_3[idx_1] - nu_3[idx_2];
				Dtype diff3s = (diff_3 * diff_3) / (2.0*time);
				if (diff3s >= thres_zero) continue;

				Dtype diff_1 = mu_1[idx_1] - nu_1[idx_2];
                Dtype diff_2 = mu_2[idx_1] - nu_2[idx_2];

                Dtype mu_nu_squared = diff_1 * diff_1 + diff_2 * diff_2;

				auto val = mu_nu_squared / (2.0 * time) + diff3s;
                Dtype tmp = a[idx_1] * b[idx_2];
				
				if (val >= thres_zero) continue;
				if (val > thres_one) {
					tmp *= std::exp(-val);
				}
                local_sum +=  tmp;
            }
            gathered_values[idx_1] = local_sum;
        });
        const Dtype sum = std::accumulate(gathered_values.begin(), gathered_values.end(), 0.0, std::plus<Dtype>());
        result = (sum * M_1_PI)/ (4.0 * time);
        return true;
    }

	template <class Dtype>
	KernelStatsUtils_D64_API bool L2DiagramMskInnerProduct(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
		const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype T1, Dtype T2)
	{
		float rate = std::sqrt(T1 / T2);
		std::vector<std::tuple<Dtype, Dtype, Dtype>> barvec1;
		std::vector<std::tuple<Dtype, Dtype, Dtype>> barvec2;
		for (auto tmp : barcodes_1) {
			barvec1.push_back(std::make_tuple(std::get<0>(tmp), std::get<1>(tmp), rate * std::get<2>(tmp)));
		}

		for (auto tmp : barcodes_2) {
			barvec2.push_back(std::make_tuple(std::get<0>(tmp), std::get<1>(tmp), rate * std::get<2>(tmp)));
		}
		L2DiagramMskInnerProduct(result, barvec1, barvec2, T1);
		return true;
	}

    template <class Dtype>
    KernelStatsUtils_D64_API bool L2DiagramMskDistanceSquare(Dtype& result, const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_1,
        const std::vector<std::tuple<Dtype, Dtype, Dtype>>& barcodes_2, Dtype time)
    {
        std::vector< std::tuple<Dtype, Dtype, Dtype> > barcodes = barcodes_1;
        for (auto b : barcodes_2) {
            auto birth = std::get<0>(b);
            auto death = std::get<1>(b);
            auto tau = std::get<2>(b);
            // reverse (birth, death) to be (death, birth)
            barcodes.push_back(std::make_tuple(death, birth, tau));
        }
        L2DiagramMskInnerProduct(result, barcodes, barcodes, time);

        return true;
    }

	KernelStatsUtils_D64_API bool L2DiagramMskInnerProductSSE(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
		const std::vector<std::tuple<float, float, float>>& barcodes_2, float time)
	{
        if (time <= 0)
            return false;
		std::vector<float> a;
		std::vector<__m128> mu;

		std::vector<float> b;
		std::vector<__m128> nu;

		int64_t num_points_1 = (int64_t)barcodes_1.size();
		int64_t num_points_2 = (int64_t)barcodes_2.size();

		for (auto tmp : barcodes_1) {
			a.push_back(1.0);
			mu.push_back(_mm_set_ps(0.0, std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp)));

			// mirrored points
			a.push_back(-1.0);
			mu.push_back(_mm_set_ps(0.0, std::get<1>(tmp), std::get<0>(tmp), std::get<2>(tmp)));
		}

		for (auto tmp : barcodes_2) {
			b.push_back(1.0);
			nu.push_back(_mm_set_ps(0.0, std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp)));

			// mirrored points
			b.push_back(-1.0);
			nu.push_back(_mm_set_ps(0.0, std::get<1>(tmp), std::get<0>(tmp), std::get<2>(tmp)));
		}

		float local_sum = 0.0;
		size_t local_end = 2 * num_points_1;
		std::vector<float> gathered_values(local_end, 0);
		float thres = 70;
		float thres2 = thres * 2 * time;
		concurrency::parallel_for((size_t)0, local_end, [&](size_t idx_1) {
			float local_sum = 0;
			for (size_t idx_2 = 0; idx_2 < 2 * num_points_2; idx_2++) {
				__m128 diff = _mm_sub_ps(mu[idx_1], nu[idx_2]);
				diff = _mm_mul_ps(diff, diff);
				if (diff.m128_f32[0] >= thres2) continue;
				diff = _mm_hadd_ps(diff, diff);
				diff = _mm_hadd_ps(diff, diff);
				float tmp = a[idx_1] * b[idx_2];
				auto val = diff.m128_f32[0] / (2.0 * time);
				if (val >= thres) continue;
				local_sum += std::exp(-val) * tmp;
			}
			gathered_values[idx_1] = local_sum;
		});
		const float sum = std::accumulate(gathered_values.begin(), gathered_values.end(), 0.0, std::plus<float>());
		result = (sum * M_1_PI) / (4.0 * time);
		return true;
	}
    
	KernelStatsUtils_D64_API bool L2DiagramMskInnerProductSSE(float& result, const std::vector<std::tuple<float, float, float>>& barcodes_1,
        const std::vector<std::tuple<float, float, float>>& barcodes_2, float T1, float T2)
    {
        if (T1 <= 0 || T2 <= 0)
            return false;
        float rate = std::sqrt(T1 / T2);
		std::vector<std::tuple<float, float, float>> barvec1;
		std::vector<std::tuple<float, float, float>> barvec2;
		for (auto tmp : barcodes_1) {
			barvec1.push_back(std::make_tuple(std::get<0>(tmp), std::get<1>(tmp), rate * std::get<2>(tmp)));
		}

		for (auto tmp : barcodes_2) {
			barvec2.push_back(std::make_tuple(std::get<0>(tmp), std::get<1>(tmp), rate * std::get<2>(tmp)));
		}
		L2DiagramMskInnerProductSSE(result, barvec1, barvec2, T1);
        return true;
    }
	
	
}
