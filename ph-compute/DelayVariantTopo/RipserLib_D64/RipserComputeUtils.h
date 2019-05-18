//******************************************************************************************************************
//     The core C++ code is derived from Ripser, 
//     which is also available under an MIT license and copyright to Ulrich Baeur.
//        Bauer, Ulrich. "Ripser: a lean C++ code for the computation of Vietoris–Rips persistence barcodes." 
//        Software available at https ://github. com/Ripser/ripser (2017).
//
//******************************************************************************************************************

#pragma once
#include "../PersistenceUtils_D64/PersistenceDeclarations.h"
#include "../PersistenceUtils_D64/PersistenceUtils_D64.h"
#include "../RipserLib_D64/RipserLib_D64.h"
#include "../RipserLib_D64/RipserUtils.h"
#include "../TopoUtils_D64/DataUtils.h"
#include <numeric>

using namespace NRipserUtils;
using namespace NDataUtils;

namespace NRipserComputeUtils {
    enum FILEFORMAT {
        LOWER_DISTANCE_MATRIX,
        UPPER_DISTANCE_MATRIX,
        DISTANCE_MATRIX,
        POINT_CLOUD,
        DIPHA
    };

    enum CompressedMatrixLayout {
        LOWER_TRIANGULAR,
        UPPER_TRIANGULAR
    };

    template <typename Entry> struct greater_diameter_or_smaller_index {
        bool operator()(const Entry& a, const Entry& b) {
            return (getDiameter(a) > getDiameter(b)) ||
                ((getDiameter(a) == getDiameter(b)) && (getIndex(a) < getIndex(b)));
        }
    };

    template <typename T> T read(std::istream& s) {
        T result;
        s.read(reinterpret_cast<char*>(&result), sizeof(T));
        return result; // on little endian: boost::endian::little_to_native(result);
    }

    class CEuclideanDistanceMatrix {
    public:
        std::vector<std::vector<value_t>> points;

        CEuclideanDistanceMatrix(std::vector<std::vector<value_t>>&& _points) : points(_points) {}

        value_t operator()(const index_t i, const index_t j) const {
            return std::sqrt(std::inner_product(points[i].begin(), points[i].end(), points[j].begin(), value_t(),
                std::plus<value_t>(),
                [](value_t u, value_t v) { return (u - v) * (u - v); }));
        }

        size_t size() const { return points.size(); }
    };

    template <typename DistanceMatrix>
    class CRipsFiltrationComparator {
        public:
            CRipsFiltrationComparator(const DistanceMatrix& _dist, const index_t _dim,
                const CBinomialCoeffTable& _binomial_coeff)
                : m_dist(_dist), m_dim(_dim), m_vertices(_dim + 1), m_binomial_coeff(_binomial_coeff) {};

            value_t diameter(const index_t index) const {
                value_t diam = 0;
                getSimplexVertices(index, m_dim, m_dist.size(), m_binomial_coeff, m_vertices.begin());

                for (index_t i = 0; i <= m_dim; ++i)
                    for (index_t j = 0; j < i; ++j) { diam = std::max(diam, m_dist(m_vertices[i], m_vertices[j])); }
                return diam;
            };

            bool operator()(const index_t a, const index_t b) const {
                assert(a < binomial_coeff(dist.size(), dim + 1));
                assert(b < binomial_coeff(dist.size(), dim + 1));

                return greater_diameter_or_smaller_index<CDiameterIndexT>()(CDiameterIndexT(diameter(a), a),
                    CDiameterIndexT(diameter(b), b));
            };

            template <typename Entry> bool operator()(const Entry& a, const Entry& b) const {
                return operator()(getIndex(a), getIndex(b));
            };

        public:
            const DistanceMatrix& m_dist;
            const index_t m_dim;

        private:
            mutable std::vector<index_t> m_vertices;
            const CBinomialCoeffTable& m_binomial_coeff;
    };

    template <class DistanceMatrix>
    class CSimplexCoboundaryEnumerator {
    private:
        index_t m_idx_below, m_idx_above, m_v, m_k;
        std::vector<index_t> m_vertices;
        const diameter_entry_t m_simplex;
        const coefficient_t m_modulus;
        const DistanceMatrix& m_dist;
        const CBinomialCoeffTable& m_binomial_coeff;

    public:
        CSimplexCoboundaryEnumerator(const diameter_entry_t _simplex, index_t _dim, index_t _n,
            const coefficient_t _modulus, const DistanceMatrix& _dist,
            const CBinomialCoeffTable& _binomial_coeff)
            : m_simplex(_simplex), m_idx_below(getIndex(_simplex)), m_idx_above(0), m_v(_n - 1), m_k(_dim + 1), m_modulus(_modulus),
            m_binomial_coeff(_binomial_coeff), m_dist(_dist), m_vertices(_dim + 1) {
            getSimplexVertices(getIndex(_simplex), _dim, _n, m_binomial_coeff, m_vertices.begin());
        }

        bool has_next() {
            while ((m_v != -1) && (m_binomial_coeff(m_v, m_k) <= m_idx_below)) {
                m_idx_below -= m_binomial_coeff(m_v, m_k);
                m_idx_above += m_binomial_coeff(m_v, m_k + 1);

                --m_v;
                --m_k;
                assert(m_k != -1);
            }
            return m_v != -1;
        }

        index_t next_index() { return idx_above + binomial_coeff(m_v--, m_k + 1) + m_idx_below; }

        diameter_entry_t next() {
            value_t coface_diameter = getDiameter(m_simplex);
            for (index_t w : m_vertices) coface_diameter = std::max(coface_diameter, m_dist(m_v, w));
            coefficient_t coface_coefficient = (m_k & 1 ? -1 + m_modulus : 1) * getCoefficient(m_simplex) % m_modulus;
            return diameter_entry_t(coface_diameter, m_idx_above + m_binomial_coeff(m_v--, m_k + 1) + m_idx_below,
                coface_coefficient);
        }
    };



    template <CompressedMatrixLayout Layout>
    class CCompressedDistanceMatrix {
    public:
        CCompressedDistanceMatrix() : m_distances(), m_rows() {};
        void init_rows();

        CCompressedDistanceMatrix(std::vector<value_t>&& _distances)
            : m_distances(_distances), m_rows((1 + std::sqrt(1 + 8 * m_distances.size())) / 2) {
            assert(m_distances.size() == size() * (size() - 1) / 2);
            init_rows();
        }

        template <typename DistanceMatrix>
        CCompressedDistanceMatrix(const DistanceMatrix& mat)
            : m_distances(mat.size() * (mat.size() - 1) / 2), m_rows(mat.size()) {
            init_rows();

            for (index_t i = 1; i < size(); ++i)
                for (index_t j = 0; j < i; ++j) m_rows[i][j] = mat(i, j);
        }

        value_t operator()(const index_t i, const index_t j) const;

        size_t size() const { return m_rows.size(); }
        std::vector<value_t> Distances() const { return m_distances; };

    private:
        std::vector<value_t> m_distances;
        std::vector<value_t*> m_rows;
    };

    typedef CCompressedDistanceMatrix<LOWER_TRIANGULAR> CompressedLowerDistMatType;
    typedef CCompressedDistanceMatrix<UPPER_TRIANGULAR> CompressedUpperDistMatType;

    template <typename ValueType>
    class CCompressedSparseMatrix {
    public:
        size_t size() const { return m_bounds.size(); }

        typename std::vector<ValueType>::const_iterator cbegin(size_t index) const {
            assert(index < size());
            return index == 0 ? m_entries.cbegin() : m_entries.cbegin() + m_bounds[index - 1];
        }

        typename std::vector<ValueType>::const_iterator cend(size_t index) const {
            assert(index < size());
            return m_entries.cbegin() + m_bounds[index];
        }

        template <typename Iterator> void append_column(Iterator begin, Iterator end) {
            for (Iterator it = begin; it != end; ++it) { m_entries.push_back(*it); }
            m_bounds.push_back(m_entries.size());
        }

        void append_column() { bounds.push_back(entries.size()); }

        void push_back(ValueType e) {
            assert(0 < size());
            m_entries.push_back(e);
            ++bounds.back();
        }

        void pop_back() {
            assert(0 < size());
            m_entries.pop_back();
            --bounds.back();
        }

        template <typename Collection> void append_column(const Collection collection) {
            append_column(collection.cbegin(), collection.cend());
        }
    private:
        std::vector<size_t> m_bounds;
        std::vector<ValueType> m_entries;
    };
};

namespace NRipserComputeUtils {
    template <typename Entry> struct smaller_index {
        bool operator()(const Entry& a, const Entry& b) { return getIndex(a) < getIndex(b); }
    };

    class union_find {
        std::vector<index_t> parent;
        std::vector<uint8_t> rank;

    public:
        union_find(index_t n) : parent(n), rank(n, 0) {
            for (index_t i = 0; i < n; ++i) parent[i] = i;
        }

        index_t find(index_t x) {
            index_t y = x, z = parent[y];
            while (z != y) {
                y = z;
                z = parent[y];
            }
            y = parent[x];
            while (z != y) {
                parent[x] = z;
                x = y;
                y = parent[x];
            }
            return z;
        }
        void link(index_t x, index_t y) {
            x = find(x);
            y = find(y);
            if (x == y) return;
            if (rank[x] > rank[y])
                parent[y] = x;
            else {
                parent[x] = y;
                if (rank[x] == rank[y]) ++rank[y];
            }
        }
    };
    
    enum class INOUT_MODE {
        INNER_MODE = 0,
        FILE_MODE = 1,
    };

    struct RipserLib_D64_API RipsPrm {
        RipsPrm() : modulus(2), dim_max(0) {};
        void ParamCopy(boost::shared_ptr<RipsPrm> _rp) {
            modulus = _rp->modulus;
            dim_max = _rp->dim_max;
            threshold = _rp->threshold;
        }
        coefficient_t modulus = 2;
        index_t dim_max = 0;
        value_t threshold = std::numeric_limits<value_t>::max();
    };
    typedef boost::shared_ptr<RipsPrm> RipsPrmPtr;

    struct RipserLib_D64_API InputPrm {
        InputPrm(){};
        InputPrm(size_t _dim) : dim_point(_dim) {
            pcl = CPointCloudPtr<value_t>::Type(new CPointCloud<value_t>(dim_point));
        }
        void ParamCopy(boost::shared_ptr<InputPrm> _ip) {
            read_mode = _ip->read_mode;
            format = _ip->format;
            input_file = _ip->input_file;
        }
        void ResetPointCloud() {
            pcl.reset(new CPointCloud<value_t>());
        }
        INOUT_MODE read_mode = INOUT_MODE::FILE_MODE;
        FILEFORMAT format = POINT_CLOUD;
        std::wstring input_file = L"Pointcloud";
        size_t dim_point = 0;
        CPointCloudPtr<value_t>::Type pcl = CPointCloudPtr<value_t>::Type(new CPointCloud<value_t>());
        std::vector<value_t> distances;
    };
    typedef boost::shared_ptr<InputPrm> InputPrmPtr;

    struct RipserLib_D64_API OutputPrm {
        OutputPrm(){};
        INOUT_MODE write_mode = INOUT_MODE::FILE_MODE;
        std::wstring out_dir = L"Output";
        std::wstring basename = L"basename";

        std::vector<CPersistenceBarcodesPtr> Barcodes() {
            return out_barcodes;
        };
        void ParamCopy(boost::shared_ptr<OutputPrm> _op) {
            write_mode = _op->write_mode;
            out_dir = _op->out_dir;
            basename = _op->basename;
        }

        bool AddBarcodes(CPersistenceBarcodesPtr bar_ptr) {
            if (!bar_ptr) return false;
            out_barcodes.push_back(bar_ptr);
            return true;
        }

        void ClearBarcodes() {
            out_barcodes.clear();
        }
    private:
        std::vector<CPersistenceBarcodesPtr> out_barcodes = std::vector<CPersistenceBarcodesPtr>();
    };
    typedef boost::shared_ptr<OutputPrm> OutputPrmPtr;

    struct RipserLib_D64_API RipsComputePrm {
        RipsComputePrm() {};
        void ParamCopy(boost::shared_ptr<RipsComputePrm> _rcp) {
            rip_prm->ParamCopy(_rcp->rip_prm);
            input_prm->ParamCopy(_rcp->input_prm);
            output_prm->ParamCopy(_rcp->output_prm);
        }

        RipsPrmPtr rip_prm = RipsPrmPtr(new RipsPrm());
        InputPrmPtr input_prm = InputPrmPtr(new InputPrm());
        OutputPrmPtr output_prm = OutputPrmPtr(new OutputPrm());
    };
    typedef boost::shared_ptr<RipsComputePrm> RipsComputePrmPtr;
    
    RipserLib_D64_API bool ComputeRipPH(RipsComputePrmPtr prm);
    RipserLib_D64_API bool ComputeRipPHMultiFiles(int nthreads, std::vector<RipsComputePrmPtr> prm_vec);
}