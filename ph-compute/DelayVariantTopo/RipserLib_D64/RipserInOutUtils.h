#pragma once
#include "RipserComputeUtils.h"
#include "../TopoUtils_D64/DataUtils.h"

using namespace NRipserComputeUtils;
using namespace NDataUtils;

namespace NRipserInOutUtils {
    CompressedLowerDistMatType read_point_cloud(std::istream& input_stream);
    CompressedLowerDistMatType read_lower_distance_matrix(std::istream& input_stream);
    CompressedLowerDistMatType read_upper_distance_matrix(std::istream& input_stream);
    CompressedLowerDistMatType read_distance_matrix(std::istream& input_stream);
    CompressedLowerDistMatType read_dipha(std::istream& input_stream);
    CompressedLowerDistMatType read_file(std::istream& input_stream, FILEFORMAT format);
    CompressedLowerDistMatType MakeDistFromPointCloud(CPointCloudPtr<value_t>::Type pcl);
}