#pragma once
#include "PersistenceUtils_D64.h"
#include <boost/tuple/tuple.hpp>
#include <boost/optional.hpp>
#include <vector>

namespace NPersistenceUtils {
    template <typename Dtype>
    class PersistenceUtils_D64_API CPersistenceBarcodes {
    public:
        CPersistenceBarcodes() {};
        CPersistenceBarcodes(const std::wstring fileName, const boost::tuple<int, Dtype, bool, Dtype> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
            const boost::tuple<Dtype, Dtype, bool> kprms /*kernel parameter opt_delay, tau_max, tau_flag*/);
        ~CPersistenceBarcodes();
        size_t dim() const { return m_dim; };
        std::wstring GetFile() { return m_filesrc; };

        std::vector<std::tuple<Dtype, Dtype, Dtype>> barcodes() { return m_barcodes; };
        const size_t numbars() { return m_barcodes.size(); };
        const bool IsEmpty() { return m_barcodes.empty(); };

        std::vector<Dtype> births();
        std::vector<Dtype> deaths();

        bool PushBar(const Dtype birth, Dtype death) {
            m_barcodes.push_back(std::make_tuple(birth, death, (Dtype)0.0)); //dummy tau
            return true;
        }

        bool SetDim(size_t _dim) {
            m_dim = _dim;
            return true;
        }

        std::wstring ToFile(const std::wstring outpath, const std::wstring fileInitial);
        Dtype GetOptimalTimeHole();
        Dtype GetOptimalTimeTau();
    private:
        std::vector<std::tuple<Dtype, Dtype, Dtype>> m_barcodes = {}; // vector of tuple(birth, death, tau) // tau for timeseries case
        size_t m_dim = 0;
        std::wstring m_filesrc = L"";
    };
}

namespace NPersistenceUtils {
    PersistenceUtils_D64_API bool dummy_func();

    template <typename Dtype>
    PersistenceUtils_D64_API bool MakeDiagramVec(std::vector<boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<Dtype>>> &diagram_vec,
        const std::wstring barcodes_path,
        const boost::tuple<int, Dtype, bool, Dtype> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
        const boost::tuple<Dtype, Dtype, bool> kpm /*opt_tau, tau_max, all_tau*/);

    template PersistenceUtils_D64_API bool MakeDiagramVec(std::vector<boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<float>>> &diagram_vec,
        const std::wstring barcodes_path,
        const boost::tuple<int, float, bool, float> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
        const boost::tuple<float, float, bool> kpm /*opt_tau, tau_max, all_tau*/);

    template PersistenceUtils_D64_API bool MakeDiagramVec(std::vector<boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<double>>> &diagram_vec,
        const std::wstring barcodes_path,
        const boost::tuple<int, double, bool, double> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
        const boost::tuple<double, double, bool> kpm /*opt_tau, tau_max, all_tau*/);
}