#pragma once

#include "stdafx.h"
#include "KernelStatsUtils_D64.h"
#include <ppl.h>
#include <Eigen/Dense>

namespace NKernelStatsUtils {
    using namespace Eigen;
    typedef std::vector<double> KFDRVec;

	// compute KFDR(n,k), 1 <= k < n
    template <class Dtype>
	KernelStatsUtils_D64_API double ComputeKFDR(int n, int k, const Matrix<Dtype, Dynamic, Dynamic> &Gk, const double lambda);

	// return arrays that contains n-1 kfdr values
    template <class Dtype>
	KernelStatsUtils_D64_API KFDRVec ComputeKFDRs(int n /*number of points*/,
        const Matrix<Dtype, Dynamic, Dynamic>&Gk, const double lambda);
}
