//******************************************************************************************************************
//     The core C++ code is derived from Ripser, 
//     which is also available under an MIT license and copyright to Ulrich Baeur.
//        Bauer, Ulrich. "Ripser: a lean C++ code for the computation of Vietoris–Rips persistence barcodes." 
//        Software available at https ://github. com/Ripser/ripser (2017).
//
//******************************************************************************************************************

#pragma once
#
#include "RipserLib_D64.h"
#include <vector>
#include <unordered_map>

//#define USE_GOOGLE_HASHMAP
#ifdef USE_GOOGLE_HASHMAP
#include <sparsehash/sparse_hash_map>
template <class Key, class T> class hash_map : public google::sparse_hash_map<Key, T> {
public:
    inline void reserve(size_t hint) { this->resize(hint); }
};
#else
template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};
#endif


typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;
typedef index_t entry_t;

namespace NRipserUtils {
#ifdef USE_COEFFICIENTS
    struct entry_t {
        index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
        coefficient_t coefficient;
        entry_t(index_t _index, coefficient_t _coefficient) : index(_index), coefficient(_coefficient) {}
        entry_t(index_t _index) : index(_index), coefficient(1) {}
        entry_t() : index(0), coefficient(1) {}
    } __attribute__((packed));

    static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

    entry_t makeEntry(index_t _index, coefficient_t _coefficient) { return entry_t(_index, _coefficient); }
    index_t getIndex(entry_t e) { return e.index; }
    index_t getCoefficient(entry_t e) { return e.coefficient; }
    void setCoefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

    bool operator==(const entry_t& e1, const entry_t& e2) {
        return getIndex(e1) == getIndex(e2) && getCoefficient(e1) == getCoefficient(e2);
    }

    std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
        stream << get_index(e) << ":" << getCoefficient(e);
        return stream;
    }

#else
    const index_t getIndex(entry_t i);
    index_t getCoefficient(entry_t i);
    entry_t makeEntry(index_t _index, coefficient_t _value);
    void setCoefficient(index_t& e, const coefficient_t c);
#endif
}


namespace NRipserUtils {
    class CBinomialCoeffTable {
        public:
            CBinomialCoeffTable(index_t n, index_t k);
            index_t operator()(index_t n, index_t k) const;
        private:
            std::vector<std::vector<index_t>> m_binomial;
            index_t m_nmax, m_kmax;
    };

    class CDiameterIndexT : public std::pair<value_t, index_t> {
    public:
        CDiameterIndexT() : std::pair<value_t, index_t>() {}
        CDiameterIndexT(std::pair<value_t, index_t> p) : std::pair<value_t, index_t>(p) {}
    };

    RipserLib_D64_API bool isPrime(const coefficient_t n);
    
    std::vector<coefficient_t> MultiplicativeInverseVector(const coefficient_t m);
    
    index_t getNextVertex(index_t& v, const index_t idx, const index_t k,
        const CBinomialCoeffTable& binomial_coeff);

    template <typename OutputIterator>
    OutputIterator getSimplexVertices(index_t idx, const index_t dim, index_t v,
        const CBinomialCoeffTable& binomial_coeff, OutputIterator out) {
        --v;
        for (index_t k = dim + 1; k > 0; --k) {
            getNextVertex(v, idx, k, binomial_coeff);
            *out++ = v;
            idx -= binomial_coeff(v, k);
        }
        return out;
    }

    std::vector<index_t> VerticesOfSimplex(const index_t simplex_index, const index_t dim, const index_t n,
        const CBinomialCoeffTable& binomial_coeff);

    const entry_t& getEntry(const entry_t& e);
    value_t getDiameter(CDiameterIndexT i);
    index_t getIndex(CDiameterIndexT i);

    class diameter_entry_t : public std::pair<value_t, entry_t> {
    public:
        diameter_entry_t(std::pair<value_t, entry_t> p) : std::pair<value_t, entry_t>(p) {}
        diameter_entry_t(entry_t e) : std::pair<value_t, entry_t>(0, e) {}
        diameter_entry_t() : diameter_entry_t(0) {}
        diameter_entry_t(value_t _diameter, index_t _index, coefficient_t _coefficient)
            : std::pair<value_t, entry_t>(_diameter, makeEntry(_index, _coefficient)) {}
        diameter_entry_t(CDiameterIndexT _diameter_index, coefficient_t _coefficient)
            : std::pair<value_t, entry_t>(getDiameter(_diameter_index),
                makeEntry(getIndex(_diameter_index), _coefficient)) {}
        diameter_entry_t(CDiameterIndexT _diameter_index) : diameter_entry_t(_diameter_index, 1) {}
    };

    const entry_t& getEntry(const diameter_entry_t& p);
    entry_t& getEntry(diameter_entry_t& p);
    const index_t getIndex(const diameter_entry_t& p);
    const coefficient_t getCoefficient(const diameter_entry_t& p);
    const value_t& getDiameter(const diameter_entry_t& p);
    void setCoefficient(diameter_entry_t& p, const coefficient_t c);
}