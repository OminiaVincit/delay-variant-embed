// -------------------------------------------------------------------
// $ProjectName     : PersistentRunner_D64$
// $FileName        : PersistentRunner_D64.cpp$
// $Programmer      : Tran Quoc Hoan$
// -------------------------------------------------------------------

#include "stdafx.h"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

#include "../TopoUtils_D64/StringUtils.h"
#include "../TopoUtils_D64/FileUtils.h"

#include "../RipserLib_D64/RipserComputeUtils.h"
#include "../RipserLib_D64/RipserInOutUtils.h"

#include <ppl.h>
#include <ppltasks.h>

using namespace NRipserComputeUtils;
using namespace NRipserInOutUtils;

int main(int argc, char** argv) {
    using namespace boost::program_options;

    int nRetCode = 0;
    options_description description("PersistentRunner");
    description.add_options()
        ("debug,x", value<bool>()->default_value(false), "Print Debug Information")
        ("nthreads,n", value<int>()->default_value(-1), "Number of threads to use, -1: use all as possible")
        ("modulus,m", value<coefficient_t>()->default_value(2), "Compute homology with coefficients in the prime field Z/<p>Z")
        ("maxdim,d", value<index_t>()->default_value(0), "Compute persistent homology up to dimension <k>")
        ("thres,th", value<value_t>()->default_value(std::numeric_limits<value_t>::max()), "Compute Rips complexes up to diameter <t>")
        ("format,f", value<std::string>()->default_value("point-cloud"), 
            "Use the specified file format for the input. Options are point-cloud, lower-distance, upper-distance, distance, dipha")
        ("outdir,o", value<std::string>()->default_value("output"), "Output directory")
        ("input,i", value<std::string>()->default_value("pointcloud"), "Input file")
        ("multi,s", value<bool>()->default_value(false), "Use multi files in one input file")
        ("help,H", "Help: Usage PersistentRunner_D64 [options]")
        ("version,v", "v1.0")
        ;
    variables_map vm;
    store(parse_command_line(argc, argv, description), vm);
    notify(vm);

    if (vm.count("help")) {
        std::cout << description << std::endl;
        return nRetCode;
    }
    auto nthreads = vm["nthreads"].as<int>();
    if (nthreads <= 0) {
        nthreads = std::thread::hardware_concurrency();
    }
    std::cout << "Number of threads (= " << nthreads << ")" << std::endl;

    RipsComputePrmPtr prm(new RipsComputePrm());
    RipsPrmPtr rip_prm = prm->rip_prm;
    InputPrmPtr input_prm = prm->input_prm;
    OutputPrmPtr output_prm = prm->output_prm;
    if(!rip_prm || !input_prm || !output_prm)
        return nRetCode;
    input_prm->read_mode = INOUT_MODE::FILE_MODE;
    output_prm->write_mode = INOUT_MODE::FILE_MODE;

    rip_prm->modulus = vm["modulus"].as<coefficient_t>();
    if (false == NRipserUtils::isPrime(rip_prm->modulus)) {
        std::cout << "Modulus must be a prime number." << std::endl;
        return nRetCode;
    }

    rip_prm->dim_max     = vm["maxdim"].as<index_t>();
    rip_prm->threshold   = vm["thres"].as<value_t>();

    output_prm->out_dir      = NStringUtil::_s2w(vm["outdir"].as<std::string>());
    std::string file_format  = vm["format"].as<std::string>();

    NRipserComputeUtils::FILEFORMAT format = NRipserComputeUtils::FILEFORMAT::POINT_CLOUD;

    if (file_format == "lower-distance")
        format = NRipserComputeUtils::FILEFORMAT::LOWER_DISTANCE_MATRIX;
    else if (file_format == "upper-distance")
        format = NRipserComputeUtils::FILEFORMAT::UPPER_DISTANCE_MATRIX;
    else if (file_format == "distance")
        format = NRipserComputeUtils::FILEFORMAT::DISTANCE_MATRIX;
    else if (file_format == "point-cloud")
        format = NRipserComputeUtils::FILEFORMAT::POINT_CLOUD;
    else if (file_format == "dipha")
        format = NRipserComputeUtils::FILEFORMAT::DIPHA;

    input_prm->format = format;
    
    auto multi = vm["multi"].as<bool>();
    auto input_file = vm["input"].as<std::string>();
    if (!boost::filesystem::is_regular_file(input_file)) {
        std::cout << "Input file: \"" << input_file << "\" not found!!!" << std::endl;
        return nRetCode;
    }
    if (!multi) {
        input_prm->input_file = NStringUtil::_s2w(input_file);
        ComputeRipPH(prm);
    }
    else {
        std::vector<std::string> file_list;
        NFileUtils::ReadMultiFilesFromInput(file_list, input_file);
        std::vector<RipsComputePrmPtr> prm_vec;
        for (auto fname : file_list) {
            RipsComputePrmPtr prm_f(new RipsComputePrm());
            prm_f->ParamCopy(prm);
            // change the input file name
            prm_f->input_prm->input_file = NStringUtil::_s2w(fname);
            prm_vec.push_back(prm_f);
        }
        if (!prm_vec.empty()) {
            std::cout << "Number of files to compute:  " << prm_vec.size() << std::endl;
            ComputeRipPHMultiFiles(nthreads, prm_vec);
        }
        else {
            std::cout << "Nothing to compute. Please check the input." << std::endl;
        }
    }

    return nRetCode;
}
