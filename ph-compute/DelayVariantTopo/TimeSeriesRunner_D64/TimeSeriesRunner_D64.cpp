// TimeSeriesRunner_D64.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "TimeSeriesRunner_D64.h"

#include "../TopoUtils_D64/StringUtils.h"
#include "../RipserLib_D64/RipserComputeUtils.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/foreach.hpp>
#include "boost/program_options.hpp"
#include <fstream>


using namespace NRipserComputeUtils;

typedef float TimeType;
typedef float DataType;
typedef CPointCloudPtr<DataType>::Type CPointCloudPtrType;

// TODO: deal with unequal sampling data
// Firstly, we'll do with timeseries data with equal sample sampling time
class CTimeSeries {
public:
    CTimeSeries() : m_series() {};
	CTimeSeries(size_t _scale) : m_series(), m_scale(_scale) {};
    size_t GetSize() const { return m_series.size(); };
    bool MakeSimulation() { return true; };
    bool ReadFromSingleCsvFile(const std::wstring& filepath) { 
        std::ifstream read_op;
        read_op.open(filepath);
        if (!read_op.good()) {
            return false;
        }
        else {
            size_t index = 0;
            while (read_op.good()) {
                std::string line;
                std::getline(read_op, line);

                if (line.size()) {
                    std::stringstream s;
                    std::string sval;
                    s << line;
                    s >> sval;
                    const std::wstring wval(NStringUtil::_s2w(sval));

                    DataType val;
                    if (typeid(DataType) == typeid(double)) {
                        val = _wtof(wval.c_str());
                    }
                    else {
                        val = (float) _wtof(wval.c_str());
                    }
                    m_series.push_back(std::make_pair((TimeType)index, val));
                    index++;
                }
                else {
                    break;
                }
            }
        }
        read_op.close();
        return true;
    };
    bool MakeFromSingleLine(const std::wstring& line, const bool &has_class) {
        std::vector<std::wstring> liststr;
        boost::algorithm::split(liststr, line, boost::is_any_of(","));
        for (size_t i = 0; i < liststr.size(); ++i) {
            const std::wstring wval = liststr[i];
            DataType val;
            if (typeid(DataType) == typeid(double)) {
                val = _wtof(wval.c_str());
            }
            else {
                val = (float)_wtof(wval.c_str());
            }
            if (has_class && i == 0) {
                // read class number from the first
                m_class = (int)val;
                continue;
            }
            m_series.push_back(std::make_pair((TimeType)(i-has_class), val));
        }
        return true;
    }
    CPointCloudPtrType MakePointCloudDelay(size_t m, size_t tau, size_t pointmax = 0) const {
        CPointCloudPtrType pcl_pp(new CPointCloud<DataType>(m));
        if (m == 0) return pcl_pp;
        // TODO: rewrite this function to deal with TimeType
        size_t len = GetSize() / m_scale;
        if (pointmax > 0) len = std::min(len, pointmax);
        for (size_t i = 0; i + (m-1)*tau < len; ++i) {
            std::vector<DataType> elem;
            for (size_t j = 0; j < m; ++j) elem.push_back(m_series[m_scale*(i + j*tau)].second);
            pcl_pp->AddElement(elem);
        }
        return pcl_pp;
    }

    // This function shouldn't be called due to extral memory
    std::vector<CPointCloudPtrType> AllDelayPointCloud(size_t m, size_t tau_max) const {
        std::vector<CPointCloudPtrType> pcl_vec = {};
        if (m <= 0 || m < GetSize() / 2) 
            return pcl_vec;
        
        // TODO: rewrite this function to deal with TimeType
        for (size_t tau = 1; tau < tau_max; ++tau) {
            pcl_vec.push_back(MakePointCloudDelay(m, tau));
        }
        return pcl_vec;
    }
    boost::optional<int> GetClass() { return m_class; }
private:
    std::vector<std::pair<TimeType, DataType>> m_series;
    boost::optional<int> m_class = boost::none;
	size_t m_scale = 1;
};
typedef boost::shared_ptr<CTimeSeries> CTimeSeriesPtr;

std::vector<std::wstring> ReadMultiLinesFromFile(const std::wstring& filepath) {
    std::vector<std::wstring> rstr;
    std::ifstream read_op;
    read_op.open(filepath);
    if (!read_op.good()) {
        return rstr;
    }
    else {
        size_t index = 0;
        while (read_op.good()) {
            std::string line;
            std::getline(read_op, line);
            if (line.empty() == false)
                rstr.push_back(NStringUtil::_s2w(line));
        }
    }
    read_op.close();
    return rstr;
}
bool WriteDelayBarcodesToFile(std::vector<std::pair<size_t, RipsComputePrmPtr>> delay_barcodes,
    std::wstring out_dir, size_t maxdim, std::wstring basename) {
    namespace io = boost::iostreams;
    // Make output directory if not exist
    boost::filesystem::path outpath(out_dir);
    if (!boost::filesystem::is_directory(outpath)) {
        boost::filesystem::create_directory(outpath);
    }
    for (size_t dim = 0; dim <= maxdim; ++dim) {
        std::wstring outfile = out_dir + L"/delay_barcode_" + basename + L"_dim_" + std::to_wstring(dim) + L".txt";
        io::stream_buffer<io::file_sink> buf(NStringUtil::_w2s(outfile));
        std::ostream out(&buf);

        for (auto dbar : delay_barcodes) {
            size_t tau = dbar.first;
            RipsComputePrmPtr prm = dbar.second;
            if (!prm) continue;
            OutputPrmPtr output_prm = prm->output_prm;
            if (!output_prm) continue;
            auto bvec = output_prm->Barcodes();
            if (bvec.size() <= dim) continue;
            for (auto bar : bvec[dim]->barcodes()) {
                out << std::get<0>(bar) << " " << std::get<1>(bar) << " " << tau << " " << std::endl;
            }
        }
    }
    return true;
}
int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
    using namespace boost::program_options;
    namespace fs = boost::filesystem;

    int nRetCode = 0;
    options_description description("Timeseries");
    description.add_options()
        ("debug,x", value<bool>()->default_value(false), "Print Debug Information")
        ("nthreads,n", value<int>()->default_value(-1), "Number of threads")
        ("emd,e", value<size_t>()->default_value(3), "Embedded dimension")
        ("modulus,m", value<coefficient_t>()->default_value(2), "Compute homology with coefficients in the prime field Z/<p>Z")
        ("maxdim,d", value<index_t>()->default_value(0), "Compute persistent homology up to dimension <k>")
        ("threshold,th", value<value_t>()->default_value(std::numeric_limits<value_t>::max()), "Compute Rips complexes up to diameter <t>")
        ("output_dir,o", value<std::string>()->default_value("Output"), "Output directory")
        ("input,i", value<std::string>()->default_value("timeseries"), "Input time series file or folder")
        ("multi,s", value<bool>()->default_value(false), "Use multi files in one input file")
        ("taumax,tm", value<size_t>()->default_value(0), "Compute with delay time upt to taumax")
        ("numpoints,p", value<size_t>()->default_value(0), "Compute up to number of points")
        ("scale,S", value<size_t>()->default_value(1), "Scale sampling in time series")
        ("numskip,N", value<size_t>()->default_value(0), "Skip time series")
        ("help,H", "Help: Usage TimeSeriesRunner [options] filename")
        ("version,v", "v1.0")
        ;
    variables_map vm;
    store(parse_command_line(argc, argv, description), vm);
    notify(vm);

    if (vm.count("help")) {
        std::cout << description << std::endl;
        return nRetCode;
    }
    size_t emb_dim = vm["emd"].as<size_t>();
    if (emb_dim <= 1)
        return false;
    size_t taumax = vm["taumax"].as<size_t>();
    size_t numpoints = vm["numpoints"].as<size_t>();
    size_t scale = vm["scale"].as<size_t>();
    if(scale <= 0)
        return false;
    bool multi_in_one = vm["multi"].as<bool>();
    size_t numskip = vm["numskip"].as<size_t>();

    RipsComputePrmPtr prm(new RipsComputePrm());
    RipsPrmPtr rip_prm = prm->rip_prm;
    InputPrmPtr input_prm = prm->input_prm;
    OutputPrmPtr output_prm = prm->output_prm;
    if(rip_prm == NULL || input_prm == NULL || output_prm == NULL)
        return false;

    auto maxdim = vm["maxdim"].as<index_t>();
    auto out_dir = NStringUtil::_s2w(vm["output_dir"].as<std::string>());

    rip_prm->modulus = vm["modulus"].as<coefficient_t>();
    if(false == NRipserUtils::isPrime(rip_prm->modulus))
        return false;
    rip_prm->dim_max = maxdim;
    rip_prm->threshold = vm["threshold"].as<value_t>();

    output_prm->write_mode = INOUT_MODE::INNER_MODE;
    output_prm->out_dir = out_dir;
    input_prm->read_mode = INOUT_MODE::INNER_MODE;

    // read all input time series file
    std::vector<fs::path> input_files;
    fs::path input_path(vm["input"].as<std::string>());
    std::vector<std::wstring> rstr;

    if (fs::is_directory(input_path)) {
        // loop for all files in folder
        fs::directory_iterator it(input_path), eod;
        BOOST_FOREACH(fs::path const &p, std::make_pair(it, eod)) {
            if (fs::is_regular_file(p)) {
                input_files.push_back(p);
            }
        }
    }
    else if (fs::is_regular_file(input_path)) {
        if (multi_in_one) {
            rstr = ReadMultiLinesFromFile(input_path.c_str());
            // add dummy path
            auto p = input_path;
            auto extstr = p.extension().c_str();
            p.replace_extension();

            for (size_t i = 0; i < rstr.size(); ++i) {
                std::wstring tmpstr = p.c_str();
                input_files.push_back(fs::path(tmpstr + L"_" + std::to_wstring(i) + extstr));
            }
        }
        else {
            input_files.push_back(input_path);
        }
    }
    else {
        return nRetCode;
    }
    for (size_t i = 0; i < input_files.size(); ++i) {
        if (i < numskip) continue;
        auto p = input_files[i];
        CTimeSeriesPtr ts(new CTimeSeries(scale));

        if (multi_in_one) {
            if (i >= rstr.size())
                break;
            ts->MakeFromSingleLine(rstr[i], false);
        }
        else {
            ts->ReadFromSingleCsvFile(p.c_str());
        }

        input_prm->dim_point = emb_dim;
        p.replace_extension();
        std::wstring basename = p.filename().c_str();
        std::vector<std::pair<size_t, RipsComputePrmPtr>> delay_barcodes;
        std::vector<RipsComputePrmPtr> prm_vec;
        if (taumax == 0) taumax = ts->GetSize() / (scale * (emb_dim - 1));
        for (size_t tau = 1; tau <= taumax; ++tau) {
            RipsComputePrmPtr prm_f(new RipsComputePrm());
            prm_f->ParamCopy(prm);
            prm_f->input_prm->pcl = ts->MakePointCloudDelay(input_prm->dim_point, tau, numpoints);
            if (prm_f->input_prm->pcl->Nums() < 2) continue;

            if (output_prm->write_mode == INOUT_MODE::FILE_MODE) {
                output_prm->basename = L"tau_" + std::to_wstring(tau);
            }
            prm_vec.push_back(prm_f);
            delay_barcodes.push_back(std::make_pair(tau, prm_f));
        }
        if (!prm_vec.empty()) {
            auto nthreads = vm["nthreads"].as<int>();
            ComputeRipPHMultiFiles(nthreads, prm_vec);
        }
        WriteDelayBarcodesToFile(delay_barcodes, out_dir, maxdim, basename);
    }
	return nRetCode;
}