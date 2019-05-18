//******************************************************************************************************************
//     The core C++ code is derived from Ripser, 
//     which is also available under an MIT license and copyright to Ulrich Baeur.
//        Bauer, Ulrich. "Ripser: a lean C++ code for the computation of Vietoris–Rips persistence barcodes." 
//        Software available at https ://github. com/Ripser/ripser (2017).
//
//******************************************************************************************************************

#include "stdafx.h"
#include "RipserInOutUtils.h"

namespace NRipserInOutUtils {
    CompressedLowerDistMatType read_point_cloud(std::istream& input_stream) {
        std::vector<std::vector<value_t>> points;

        std::string line;
        value_t value;
        while (std::getline(input_stream, line)) {
            std::vector<value_t> point;
            std::istringstream s(line);
            while (s >> value) {
                point.push_back(value);
                s.ignore();
            }
            if (!point.empty()) points.push_back(point);
            assert(point.size() == points.front().size());
        }

        CEuclideanDistanceMatrix eucl_dist(std::move(points));

        index_t n = eucl_dist.size();
        std::vector<value_t> distances;

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

        return CompressedLowerDistMatType(std::move(distances));
    }

    CompressedLowerDistMatType read_lower_distance_matrix(std::istream& input_stream) {
        std::vector<value_t> distances;
        value_t value;
        while (input_stream >> value) {
            distances.push_back(value);
            input_stream.ignore();
        }

        return CompressedLowerDistMatType(std::move(distances));
    }

    CompressedLowerDistMatType read_upper_distance_matrix(std::istream& input_stream) {
        std::vector<value_t> distances;
        value_t value;
        while (input_stream >> value) {
            distances.push_back(value);
            input_stream.ignore();
        }

        return CompressedLowerDistMatType(CompressedUpperDistMatType(std::move(distances)));
    }

    CompressedLowerDistMatType read_distance_matrix(std::istream& input_stream) {
        std::vector<value_t> distances;

        std::string line;
        value_t value;
        for (int i = 0; std::getline(input_stream, line); ++i) {
            std::istringstream s(line);
            for (int j = 0; j < i && s >> value; ++j) {
                distances.push_back(value);
                s.ignore();
            }
        }

        return CompressedLowerDistMatType(std::move(distances));
    }

    CompressedLowerDistMatType read_dipha(std::istream& input_stream) {
        if (read<int64_t>(input_stream) != 8067171840) {
            std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
            exit(-1);
        }

        if (read<int64_t>(input_stream) != 7) {
            std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
            exit(-1);
        }

        index_t n = read<int64_t>(input_stream);

        std::vector<value_t> distances;

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i > j)
                    distances.push_back(read<double>(input_stream));
                else
                    read<double>(input_stream);

        return CompressedLowerDistMatType(std::move(distances));
    }

    CompressedLowerDistMatType read_file(std::istream& input_stream, FILEFORMAT format) {
        switch (format) {
        case LOWER_DISTANCE_MATRIX:
            return read_lower_distance_matrix(input_stream);
        case UPPER_DISTANCE_MATRIX:
            return read_upper_distance_matrix(input_stream);
        case DISTANCE_MATRIX:
            return read_distance_matrix(input_stream);
        case POINT_CLOUD:
            return read_point_cloud(input_stream);
        case DIPHA:
            return read_dipha(input_stream);
        }
        return read_point_cloud(input_stream);
    }

    CompressedLowerDistMatType MakeDistFromPointCloud(CPointCloudPtr<value_t>::Type pcl) {
        CEuclideanDistanceMatrix eucl_dist(pcl->Pcl());
        index_t n = eucl_dist.size();
        std::vector<value_t> distances;

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

        return CompressedLowerDistMatType(std::move(distances));
    }
}