#include "stdafx.h"
#include <sstream>
#include <fstream>
#include <ostream>
#include "PersistenceBarcodes.h"
#include "../TopoUtils_D64/StringUtils.h"
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/tuple/tuple.hpp>

namespace NPersistenceUtils {
    namespace fs = boost::filesystem;

    template <typename Dtype>
    CPersistenceBarcodes<Dtype>::~CPersistenceBarcodes() {
    };

    template <typename Dtype>
    CPersistenceBarcodes<Dtype>::CPersistenceBarcodes(const std::wstring filename,
        const boost::tuple<int, Dtype, bool, Dtype> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
        const boost::tuple<Dtype, Dtype, bool> kprms /*opt_tau, taumax, all_tau*/) {
        auto dim = rprms.get<0>();
        auto valInfty = rprms.get<1>();
        auto skipInfty = rprms.get<2>();
        auto threshold = rprms.get<3>();

        auto opt_tau = kprms.get<0>();
        auto tau_max = kprms.get<1>();
        auto all_tau = kprms.get<2>();

        std::wifstream read_op;
        read_op.open(NStringUtil::_w2s(filename));
        if (!read_op.good()) {
            return;
        }
        else {
            while (read_op.good()) {
                std::wstring line;
                std::getline(read_op, line);

                if (line.size()) {
                    std::wstringstream s;
                    std::wstring wbirth, wdeath, wtau;
                    s << line;
                    s >> wbirth;
                    s >> wdeath;
                    s >> wtau;

                    Dtype birth, death, tau;
                    if (typeid(Dtype) == typeid(double)) {
                        birth = _wtof(wbirth.c_str());
                        death = _wtof(wdeath.c_str());
                        tau   = _wtof(wtau.c_str());
                    }
                    else {
                        birth = (float) _wtof(wbirth.c_str());
                        death = (float) _wtof(wdeath.c_str());
                        tau   = (float) _wtof(wtau.c_str());
                    }
                    // read with specified tau
                    if (all_tau == true && tau_max > 0 && tau > tau_max)
                        break;
                    if (all_tau == false && opt_tau != tau)
                        continue;
                    if (skipInfty && death == std::numeric_limits<Dtype>::infinity()) {
                        continue;
                    }
                    if (death == std::numeric_limits<Dtype>::infinity()) death = std::max(birth, valInfty);
                    if (abs(death - birth) <= threshold) continue;
                    if (all_tau == false) tau = 0.0;
                    m_barcodes.push_back(std::make_tuple(birth, death, tau));
                }
                else {
                    break;
                }
            }
        }
        read_op.close();

        m_dim = static_cast<size_t>(dim);
        m_filesrc = filename;
    }

    template <typename Dtype>
    std::vector<Dtype> CPersistenceBarcodes<Dtype>::births() {
        std::vector<Dtype> birth_vecs;
        for (auto bar : m_barcodes) birth_vecs.push_back(std::get<0>(bar));
        return birth_vecs;
    }

    template <typename Dtype>
    std::vector<Dtype> CPersistenceBarcodes<Dtype>::deaths() {
        std::vector<Dtype> death_vecs;
        for (auto bar : m_barcodes) death_vecs.push_back(std::get<1>(bar));
        return death_vecs;
    }

    template<typename Dtype>
    std::wstring CPersistenceBarcodes<Dtype>::ToFile(const std::wstring outpath, const std::wstring fileInitial)
    {
        namespace io = boost::iostreams;

        std::wstring outfile = outpath + L"/barcode_" + fileInitial + L"_dim_" + std::to_wstring(dim()) + L".txt";
        io::stream_buffer<io::file_sink> buf(NStringUtil::_w2s(outfile));
        std::ostream out(&buf);
        for (auto bar : m_barcodes) {
            out << std::get<0>(bar) << " " << std::get<1>(bar) << std::endl;
        }
        return outfile;
    }

    // Get the median of an unordered set of numbers of arbitrary 
    // type without modifying the underlying dataset.
    template <typename It>
    auto GetMedian(It begin, It end)
    {
        using T = typename std::iterator_traits<It>::value_type;
        std::vector<T> data(begin, end);
        std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
        return data[data.size() / 2];
    }

    template<typename Dtype>
    Dtype CPersistenceBarcodes<Dtype>::GetOptimalTimeHole()
    {
        std::vector<Dtype> dis_diffs;
        Dtype tHole = Dtype(0.0);
        auto n = m_barcodes.size();
        for (size_t i = 0; i < n; ++i) {
            auto bar1 = m_barcodes[i];
            for (size_t j = i + 1; j < n; ++j) {
                auto bar2 = m_barcodes[j];
                auto xdiff = std::get<0>(bar1) - std::get<0>(bar2);
                auto ydiff = std::get<1>(bar1) - std::get<1>(bar2);
                dis_diffs.push_back(xdiff * xdiff + ydiff * ydiff);
            }
        }
        if (dis_diffs.empty() == false) {
            tHole = GetMedian(dis_diffs.begin(), dis_diffs.end());
        }
        return tHole;
    }

    template<typename Dtype>
    Dtype CPersistenceBarcodes<Dtype>::GetOptimalTimeTau()
    {
        std::vector<Dtype> tau_diffs;
        Dtype tTau = Dtype(0.0);
        auto n = m_barcodes.size();
        for (size_t i = 0; i < n; ++i) {
            auto bar1 = m_barcodes[i];
            for (size_t j = i + 1; j < n; ++j) {
                auto bar2 = m_barcodes[j];
                auto tdiff = std::get<2>(bar1) - std::get<2>(bar2);
                tau_diffs.push_back(tdiff * tdiff);
            }
        }
        if (tau_diffs.empty() == false) {
            tTau = GetMedian(tau_diffs.begin(), tau_diffs.end());
        }
        return tTau;
    }

    template PersistenceUtils_D64_API class CPersistenceBarcodes<float>;
    template PersistenceUtils_D64_API class CPersistenceBarcodes<double>;
}

namespace NPersistenceUtils {
    PersistenceUtils_D64_API bool dummy_func()
    {
        return true;
    }

    bool ReadAttributesFromString(std::wstring &filename, int &tau, int &dim, std::wstring src) {
        if (src.empty())
            return true;
        std::vector<std::wstring> lpaths;
        boost::algorithm::split(lpaths, src, boost::is_any_of(";"));
        if (lpaths.empty())
            return true;
        filename = lpaths[0];

        boost::iterator_range<std::wstring::iterator> rtau = boost::algorithm::find_regex(src, boost::regex{ "tau=\\d+" });
        boost::iterator_range<std::wstring::iterator> rdim = boost::algorithm::find_regex(src, boost::regex{ "dim=\\d+" });

        if (rtau.empty() == false) {
            auto sstr = std::wstring(rtau.begin() + 4, rtau.end());
            tau = static_cast<int>(_wtof(sstr.c_str()));
        }

        if (rdim.empty() == false) {
            auto sstr = std::wstring(rdim.begin() + 4, rdim.end());
            dim = static_cast<int>(_wtof(sstr.c_str()));
        }
        return true;
    }

    template <typename Dtype>
    bool MakeDiagramVecFromBarcodesFile(std::vector<boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<Dtype>>> &diagram_vec,
        const std::wstring barlist_filename, const boost::tuple<int, Dtype, bool, Dtype> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
        const boost::tuple<Dtype, Dtype, bool> kpm /*opt_tau, tau_max, all_tau*/) {
        const int dim = rprms.get<0>();
        const Dtype valInfty = rprms.get<1>();
        const bool  skipInfty = rprms.get<2>();
        const Dtype threshold = rprms.get<3>();

        Dtype opt_tau = kpm.get<0>();
        const Dtype tau_max = kpm.get<1>();
        const bool  all_tau = kpm.get<2>();

        // Read parameters from file
        std::wifstream read_op;
        read_op.open(barlist_filename);
        if (!read_op.good()) {
            return false;
        }
        else {
            while (read_op.good()) {
                std::wstring line;
                std::getline(read_op, line);

                if (!line.empty()) {
                    std::wstring bar_path = L"";
                    int bDim = dim;
                    boost::optional<int> classIdx = boost::none;
                    int single_tau = 1;
                    ReadAttributesFromString(bar_path, single_tau, bDim, line);
                    if (all_tau == false && opt_tau <= 0)
                        opt_tau = (Dtype)single_tau;
                    {
                        fs::path bpath(bar_path);
                        if (fs::is_regular_file(bpath) == false) {
                            fs::path barlist_path(barlist_filename);
                            bpath = barlist_path.remove_filename() / bpath;
                            if (fs::is_regular_file(bpath)) bar_path = bpath.c_str();
                            else continue;
                        }
                    }
                    const boost::tuple<Dtype, Dtype, bool> kprms = boost::make_tuple(opt_tau, tau_max, all_tau);
                    const boost::tuple<int, Dtype, bool, Dtype> rprms2 = boost::make_tuple(bDim, valInfty, skipInfty, threshold);

                    boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<Dtype>> bar_ptr(new CPersistenceBarcodes<Dtype>(bar_path, rprms2, kprms));
                    diagram_vec.push_back(bar_ptr);
                }
            }
        }
        read_op.close();
        return true;
    }

    template <typename Dtype>
    PersistenceUtils_D64_API bool MakeDiagramVec(std::vector<boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<Dtype>>>& diagram_vec /*ouput vector of diagrams*/,
        const std::wstring barcodes_path,
        const boost::tuple<int, Dtype, bool, Dtype> rprms /*hole parameter dim, valInfty, skipInfty, threshold*/,
        const boost::tuple<Dtype, Dtype, bool> kpm /*opt_tau, tau_max, all_tau*/) {
        fs::path input_path(barcodes_path);

        if (fs::is_directory(input_path)) {
            // find barcode files in folder
            fs::directory_iterator it(input_path), eod;
            BOOST_FOREACH(fs::path const &p, std::make_pair(it, eod)) {
                if (fs::is_regular_file(p)) {
                    boost::shared_ptr<NPersistenceUtils::CPersistenceBarcodes<Dtype>> bar_ptr(new CPersistenceBarcodes<Dtype>(p.c_str(), rprms, boost::make_tuple((Dtype)0.0, (Dtype)1.0, false)));
                    if (bar_ptr->IsEmpty() == false) {
                        diagram_vec.push_back(bar_ptr);
                    }
                }
            }
        }
        else if (fs::is_regular_file(input_path)) {
            MakeDiagramVecFromBarcodesFile(diagram_vec, barcodes_path, rprms, kpm);
        }
        else {
            return false;
        }
        return true;
    }
}