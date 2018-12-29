#pragma once
#include "TopoUtils_D64.h"

namespace NFileUtils {
    TopoUtils_D64_API bool ReadMultiFilesFromInput(std::vector<std::string>& file_list, std::string input_file);
}