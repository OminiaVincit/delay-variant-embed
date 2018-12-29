#include "stdafx.h"
#include "RipserComputeUtils.h"
#include "RipserInOutUtils.h"
#include "../PersistenceUtils_D64/PersistenceDeclarations.h"
#include "../PersistenceUtils_D64/PersistenceUtils_D64.h"
#include "../TopoUtils_D64/FileUtils.h"

#include <boost/filesystem.hpp>

#include <ppl.h>
#include <ppltasks.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <chrono>

namespace NRipserComputeUtils {
    using namespace NRipserUtils;
    using namespace NPersistenceUtils;
    using namespace NRipserInOutUtils;

    template <> void CCompressedDistanceMatrix<LOWER_TRIANGULAR>::init_rows() {
        value_t* pointer = &m_distances[0];
        for (index_t i = 1; i < size(); ++i) {
            m_rows[i] = pointer;
            pointer += i;
        }
    }

    template <> void CCompressedDistanceMatrix<UPPER_TRIANGULAR>::init_rows() {
        value_t* pointer = &m_distances[0] - 1;
        for (index_t i = 0; i < size() - 1; ++i) {
            m_rows[i] = pointer;
            pointer += size() - i - 2;
        }
    }

    template <> value_t CCompressedDistanceMatrix<UPPER_TRIANGULAR>::operator()(index_t i, index_t j) const {
        if (i > j) std::swap(i, j);
        return i == j ? 0 : m_rows[i][j];
    }

    template <> value_t CCompressedDistanceMatrix<LOWER_TRIANGULAR>::operator()(index_t i, index_t j) const {
        if (i > j) std::swap(i, j);
        return i == j ? 0 : m_rows[j][i];
    }


    template <typename Heap> diameter_entry_t pop_pivot(Heap& column, coefficient_t modulus) {
        if (column.empty())
            return diameter_entry_t(-1);
        else {
            auto pivot = column.top();

#ifdef USE_COEFFICIENTS
            coefficient_t coefficient = 0;
            do {
                coefficient = (coefficient + getCoefficient(column.top())) % modulus;
                column.pop();

                if (coefficient == 0) {
                    if (column.empty())
                        return diameter_entry_t(-1);
                    else
                        pivot = column.top();
                }
            } while (!column.empty() && getIndex(column.top()) == getIndex(pivot));
            if (getIndex(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
            column.pop();
            while (!column.empty() && getIndex(column.top()) == getIndex(pivot)) {
                column.pop();
                if (column.empty())
                    return diameter_entry_t(-1);
                else {
                    pivot = column.top();
                    column.pop();
                }
            }
#endif
            return pivot;
        }
    }

    template <typename Heap> diameter_entry_t get_pivot(Heap& column, coefficient_t modulus) {
        diameter_entry_t result = pop_pivot(column, modulus);
        if (getIndex(result) != -1) column.push(result);
        return result;
    }

    template <typename Heap> void push_entry(Heap& column, index_t i, coefficient_t c, value_t diameter) {
        entry_t e = makeEntry(i, c);
        column.push(std::make_pair(diameter, e));
    }

    template <typename Comparator>
    void assemble_columns_to_reduce(std::vector<CDiameterIndexT>& columns_to_reduce,
        hash_map<index_t, index_t>& pivot_column_index, const Comparator& comp, index_t dim,
        index_t n, value_t threshold, const CBinomialCoeffTable& binomial_coeff) {
        index_t num_simplices = binomial_coeff(n, dim + 2);

        columns_to_reduce.clear();

#ifdef INDICATE_PROGRESS
        std::cout << "\033[K"
            << "assembling " << num_simplices << " columns" << std::flush << "\r";
#endif

        for (index_t index = 0; index < num_simplices; ++index) {
            if (pivot_column_index.find(index) == pivot_column_index.end()) {
                value_t diameter = comp.diameter(index);
                if (diameter <= threshold) columns_to_reduce.push_back(std::make_pair(diameter, index));
#ifdef INDICATE_PROGRESS
                if ((index + 1) % 1000 == 0)
                    std::cout << "\033[K"
                    << "assembled " << columns_to_reduce.size() << " out of " << (index + 1) << "/"
                    << num_simplices << " columns" << std::flush << "\r";
#endif
            }
        }

#ifdef INDICATE_PROGRESS
        std::cout << "\033[K"
            << "sorting " << num_simplices << " columns" << std::flush << "\r";
#endif

        std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
            greater_diameter_or_smaller_index<CDiameterIndexT>());
#ifdef INDICATE_PROGRESS
        std::cout << "\033[K";
#endif
    }

    template <typename DistanceMatrix, typename ComparatorCofaces, typename Comparator>
    void compute_pairs(CPersistenceBarcodesPtr& bar_ptr, std::vector<CDiameterIndexT>& columns_to_reduce, 
        hash_map<index_t, index_t>& pivot_column_index,
        index_t dim, index_t n, value_t threshold, coefficient_t modulus,
        const std::vector<coefficient_t>& multiplicative_inverse, const DistanceMatrix& dist,
        const ComparatorCofaces& comp, const Comparator& comp_prev,
        const CBinomialCoeffTable& binomial_coeff) {

#ifdef ASSEMBLE_REDUCTION_MATRIX
        CCompressedSparseMatrix<diameter_entry_t> reduction_coefficients;
#else
#ifdef USE_COEFFICIENTS
        std::vector<diameter_entry_t> reduction_coefficients;
#endif
#endif

        std::vector<diameter_entry_t> coface_entries;

        for (index_t i = 0; i < columns_to_reduce.size(); ++i) {
            auto column_to_reduce = columns_to_reduce[i];

#ifdef ASSEMBLE_REDUCTION_MATRIX
            std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>, smaller_index<diameter_entry_t>>
                reduction_column;
#endif

            std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
                greater_diameter_or_smaller_index<diameter_entry_t>>
                working_coboundary;

            value_t diameter = getDiameter(column_to_reduce);

            index_t j = i;

            // start with a dummy pivot entry with coefficient -1 in order to initialize
            // working_coboundary with the coboundary of the simplex with index column_to_reduce
            diameter_entry_t pivot(0, -1, -1 + modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
            // initialize reduction_coefficients as identity matrix
            reduction_coefficients.append_column();
            reduction_coefficients.push_back(diameter_entry_t(column_to_reduce, 1));
#else
#ifdef USE_COEFFICIENTS
            reduction_coefficients.push_back(diameter_entry_t(column_to_reduce, 1));
#endif
#endif

            bool might_be_apparent_pair = (i == j);
            do {
                const coefficient_t factor = modulus - getCoefficient(pivot);

#ifdef ASSEMBLE_REDUCTION_MATRIX
                auto coeffs_begin = reduction_coefficients.cbegin(j), coeffs_end = reduction_coefficients.cend(j);
#else
#ifdef USE_COEFFICIENTS
                auto coeffs_begin = &reduction_coefficients[j], coeffs_end = &reduction_coefficients[j] + 1;
#else
                auto coeffs_begin = &columns_to_reduce[j], coeffs_end = &columns_to_reduce[j] + 1;
#endif
#endif

                for (auto it = coeffs_begin; it != coeffs_end; ++it) {
                    diameter_entry_t simplex = *it;
                    setCoefficient(simplex, getCoefficient(simplex) * factor % modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
                    reduction_column.push(simplex);
#endif

                    coface_entries.clear();
                    CSimplexCoboundaryEnumerator<decltype(dist)> cofaces(simplex, dim, n, modulus, dist, binomial_coeff);
                    while (cofaces.has_next()) {
                        diameter_entry_t coface = cofaces.next();
                        if (getDiameter(coface) <= threshold) {
                            coface_entries.push_back(coface);
                            if (might_be_apparent_pair && (getDiameter(simplex) == getDiameter(coface))) {
                                if (pivot_column_index.find(getIndex(coface)) == pivot_column_index.end()) {
                                    pivot = coface;
                                    goto found_persistence_pair;
                                }
                                might_be_apparent_pair = false;
                            }
                        }
                    }
                    for (auto e : coface_entries) working_coboundary.push(e);
                }

                pivot = get_pivot(working_coboundary, modulus);

                if (getIndex(pivot) != -1) {
                    auto pair = pivot_column_index.find(getIndex(pivot));

                    if (pair != pivot_column_index.end()) {
                        j = pair->second;
                        continue;
                    }
                }
                else {
                    bar_ptr->PushBar(diameter, std::numeric_limits<value_t>::infinity());
                    break;
                }

            found_persistence_pair:
                value_t death = getDiameter(pivot);
                if (diameter != death) {
                    bar_ptr->PushBar(diameter, death);
                }

                pivot_column_index.insert(std::make_pair(getIndex(pivot), i));

#ifdef USE_COEFFICIENTS
                const coefficient_t inverse = multiplicative_inverse[getCoefficient(pivot)];
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
                // replace current column of reduction_coefficients (with a single diagonal 1 entry)
                // by reduction_column (possibly with a different entry on the diagonal)
                reduction_coefficients.pop_back();
                while (true) {
                    diameter_entry_t e = pop_pivot(reduction_column, modulus);
                    if (getIndex(e) == -1) break;
#ifdef USE_COEFFICIENTS
                    setCoefficient(e, inverse * getCoefficient(e) % modulus);
                    assert(getCoefficient(e) > 0);
#endif
                    reduction_coefficients.push_back(e);
                }
#else
#ifdef USE_COEFFICIENTS
                reduction_coefficients.pop_back();
                reduction_coefficients.push_back(diameter_entry_t(column_to_reduce, inverse));
#endif
#endif
                break;
            } while (true);
        }
    }

    bool WriteBarcodesToFile(CPersistenceBarcodesPtr bar_ptr, std::wstring out_dir, std::wstring basename) {
        if(!bar_ptr)
            return false;
        if (bar_ptr->IsEmpty() == true) {
            return false;
        }
        // Make output directory if not exist
        boost::filesystem::path outpath(out_dir);
        if (!boost::filesystem::is_directory(outpath)) {
            boost::filesystem::create_directory(outpath);
        }
        auto outfile_name = bar_ptr->ToFile(outpath.c_str(), basename);
        return true;
    }

    bool ComputeRipPH(RipsComputePrmPtr prm) {
        RipsPrmPtr rip_prm = prm->rip_prm;
        InputPrmPtr input_prm = prm->input_prm;
        OutputPrmPtr output_prm = prm->output_prm;
        if(!rip_prm || !input_prm || !output_prm)
            return false;
        output_prm->ClearBarcodes();

        if(!isPrime(rip_prm->modulus))
            return false;

        CompressedLowerDistMatType dist;
        std::wstring basename = output_prm->basename;

        if (input_prm->read_mode == INOUT_MODE::FILE_MODE) {
            // read from file
            // get basename for input file
            boost::filesystem::path p(input_prm->input_file);
            p.replace_extension();
            basename = p.filename().c_str();

            std::ifstream file_stream(input_prm->input_file);
            if(false==file_stream.good())
                return false;

            dist = read_file(file_stream, input_prm->format);
        }
        else {
            if (input_prm->format == POINT_CLOUD) {
                // read from memory (point cloud)
                dist = MakeDistFromPointCloud(input_prm->pcl);
            }
            else {
				// make lower mat from distance vec
                dist = CompressedLowerDistMatType(std::move(input_prm->distances));
            }
        }

        index_t n = dist.size();
        if (false) {
            auto distances = dist.Distances();
            auto value_range = std::minmax_element(distances.begin(), distances.end());
            auto min_idx = value_range.first - distances.begin();
            auto max_idx = value_range.second - distances.begin();
        }


        auto dim_max = std::min(rip_prm->dim_max, n - 2);

        CBinomialCoeffTable binomial_coeff(n, dim_max + 2);
        std::vector<coefficient_t> multiplicative_inverse(NRipserUtils::MultiplicativeInverseVector(rip_prm->modulus));

        std::vector<CDiameterIndexT> columns_to_reduce;

        // PERSISTENT HOMOLOGY with DIM = 0
        {
            union_find dset(n);
            std::vector<CDiameterIndexT> edges;
            CRipsFiltrationComparator<decltype(dist)> comp(dist, 1, binomial_coeff);
            for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
                value_t diameter = comp.diameter(index);
                if (diameter <= rip_prm->threshold) edges.push_back(std::make_pair(diameter, index));
            }
            std::sort(edges.rbegin(), edges.rend(), greater_diameter_or_smaller_index<CDiameterIndexT>());

            // Generate barcode for dim 0
            CPersistenceBarcodesPtr bar0_ptr(new CPersistenceBarcodes<FType>());

            std::vector<index_t> vertices_of_edge(2);
            for (auto e : edges) {
                vertices_of_edge.clear();
                getSimplexVertices(getIndex(e), 1, n, binomial_coeff, std::back_inserter(vertices_of_edge));
                index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

                if (u != v) {
                    auto diameter = getDiameter(e);
                    if (diameter > 0) {
                        bar0_ptr->PushBar(0, diameter);
                    }
                    dset.link(u, v);
                }
                else
                    columns_to_reduce.push_back(e);
            }
            std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

            for (index_t i = 0; i < n; ++i) {
                if (dset.find(i) == i) {
                    bar0_ptr->PushBar(0, std::numeric_limits<value_t>::infinity());
                }
            }
            output_prm->AddBarcodes(bar0_ptr);

            // To file
            if (output_prm->write_mode == INOUT_MODE::FILE_MODE) {
                WriteBarcodesToFile(bar0_ptr, output_prm->out_dir, basename);
            }
        }

        for (index_t dim = 1; dim <= dim_max; ++dim) {
            CRipsFiltrationComparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
            CRipsFiltrationComparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);

            hash_map<index_t, index_t> pivot_column_index;
            pivot_column_index.reserve(columns_to_reduce.size());
            CPersistenceBarcodesPtr bar_ptr(new CPersistenceBarcodes<FType>());
            bar_ptr->SetDim(dim);
            compute_pairs(bar_ptr, columns_to_reduce, pivot_column_index, dim, n, rip_prm->threshold, rip_prm->modulus, multiplicative_inverse, dist,
                comp, comp_prev, binomial_coeff);
            
            output_prm->AddBarcodes(bar_ptr);
            
            // To file
            if (output_prm->write_mode == INOUT_MODE::FILE_MODE) {
                WriteBarcodesToFile(bar_ptr, output_prm->out_dir, basename);
            }
            if (dim < dim_max) {
                assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, comp, dim, n, rip_prm->threshold, binomial_coeff);
            }
        }
        return true;
    }

    bool ComputeRipPHMultiFiles(int nprocs, std::vector<RipsComputePrmPtr> prm_vec) {
        int num_files = static_cast<int>(prm_vec.size());
        if (nprocs <= 0) {
            // should get all processors
            return false;
        }
        const int n_slice = (num_files + nprocs - 1) / nprocs;
        concurrency::parallel_for(0, nprocs, [&](int proc) {
            const int zbg = n_slice * proc;
            const int zed = std::min(num_files, zbg + n_slice);
			for (int i = zbg; i < zed; ++i) {
        		ComputeRipPH(prm_vec[i]);
			}
        });
        return true;
    }
}