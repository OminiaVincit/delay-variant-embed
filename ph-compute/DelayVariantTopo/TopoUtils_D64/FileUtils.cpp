#include "stdafx.h"
#include "FileUtils.h"
#include <fstream>
#include <iostream>

namespace NFileUtils {
    TopoUtils_D64_API bool ReadMultiFilesFromInput(std::vector<std::string>& file_list, std::string input_file) {
        std::ifstream read_op;
        read_op.open(input_file);
        if (!read_op.good()) {
            return false;
        }
        else {
            while (read_op.good()) {
                std::string line;
                std::getline(read_op, line);

                if (!line.empty()) {
                    file_list.push_back(line);
                }
            }
        }
        read_op.close();
        return true;
    }
}
