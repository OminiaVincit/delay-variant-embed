#pragma once

#include "stdafx.h"
#include "KernelFisherDA.h"

namespace NKernelStatsUtils {
    using namespace Eigen;

	using Eigen::MatrixXd;
    using Eigen::MatrixXf;

	// compute KFDR(n,k), 1 <= k < n
    template KernelStatsUtils_D64_API double ComputeKFDR(int n, int k, const MatrixXd &Gk, const double lambda);
    template KernelStatsUtils_D64_API double ComputeKFDR(int n, int k, const MatrixXf &Gk, const double lambda);

    template <class Dtype>
	KernelStatsUtils_D64_API double ComputeKFDR(int n, int k, const Matrix<Dtype, Dynamic, Dynamic> &Gk, const double lambda) {
        typedef Matrix<Dtype, Dynamic, Dynamic> MatrixDtype;
        typedef Matrix<Dtype, Dynamic, 1> VectorDtype;

        int m_plus = k;
        int m_minus = n - k;
        double inv_m_plus = 1.0 / m_plus;
        double inv_m_minus = 1.0 / m_minus;
        
        // construct vector a
        VectorDtype a(n);
        for (int i = 0; i < n; ++i) {
            if (i < m_plus) a[i] = inv_m_plus;
            else a[i] = -inv_m_minus;
        }

        // construct matrix J
        MatrixDtype J(n, n);
        J.setZero();
        J.block(0, 0, m_plus, m_plus) = sqrt(inv_m_plus) * (MatrixDtype::Identity(m_plus, m_plus) - MatrixDtype::Constant(m_plus, m_plus, inv_m_plus));
        J.block(m_plus, m_plus, m_minus, m_minus) = sqrt(inv_m_minus) * (MatrixDtype::Identity(m_minus, m_minus) - MatrixDtype::Constant(m_minus, m_minus, inv_m_minus));

        double kfdr_val = -1;
        auto I = MatrixDtype::Identity(n, n);

		// compute KFDR
#if 1
        auto sigma = lambda * I + J * Gk * J;
        auto alpha = (I - J * sigma.inverse() * J * Gk) * a / lambda;
        kfdr_val = alpha.transpose() * Gk * a;
        kfdr_val *= m_plus * 1.0 * m_minus / n;
#else
        auto A = I + J * J * Gk / lambda;
        auto B = I / lambda - (J / lambda) * Gk * A.inverse() * (J / lambda);
        auto alpha = (I - J * B * J * Gk) * a / lambda;
        kfdr_val = alpha.transpose() * Gk * a;
        kfdr_val *= m_plus * 1.0 * m_minus / n;
#endif
        return kfdr_val;
	}

	// return arrays that contains n-1 kfdr values
    template KernelStatsUtils_D64_API KFDRVec ComputeKFDRs(int n, const MatrixXd &Gk, const double lambda);
    template KernelStatsUtils_D64_API KFDRVec ComputeKFDRs(int n, const MatrixXf &Gk, const double lambda);

    template <class Dtype>
	KernelStatsUtils_D64_API KFDRVec ComputeKFDRs(int n /*number of points*/, const Matrix<Dtype, Dynamic, Dynamic> &Gk, const double lambda) {
        KFDRVec kfdrs(n - 1);
#if 1
		concurrency::parallel_for(1, n, [&](int k) {
			kfdrs[k - 1] = ComputeKFDR(n, k, Gk, lambda);
		});
#else
		for (int k = 1; k < n; ++k) {
			kfdrs[k - 1] = ComputeKFDR(n, k, Gk, lambda);
        }
#endif
		return kfdrs;
	}
}
