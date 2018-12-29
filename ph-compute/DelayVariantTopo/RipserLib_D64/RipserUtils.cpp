#include "stdafx.h"
#include "RipserUtils.h"

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
    const index_t getIndex(entry_t i) { return i; }
    index_t getCoefficient(entry_t i) { return 1; }
    entry_t makeEntry(index_t _index, coefficient_t _value) { return entry_t(_index); }
    void setCoefficient(index_t& e, const coefficient_t c) {}
#endif
}
namespace NRipserUtils {
    CBinomialCoeffTable::CBinomialCoeffTable(index_t n, index_t k)
    {
        m_nmax = n;
        m_kmax = k;

        m_binomial.resize(n + 1);
        for (index_t i = 0; i <= n; i++) {
            m_binomial[i].resize(k + 1);
            for (index_t j = 0; j <= std::min(i, k); j++) {
                if (j == 0 || j == i)
                    m_binomial[i][j] = 1;
                else
                    m_binomial[i][j] = m_binomial[i - 1][j - 1] + m_binomial[i - 1][j];
            }
        }
    }

    index_t CBinomialCoeffTable::operator()(index_t n, index_t k) const
    {
        assert(n <= m_nmax);
        assert(k <= m_kmax);
        return m_binomial[n][k];
    }

    bool isPrime(const coefficient_t n) {
        if (!(n & 1) || n < 2) return n == 2;
        for (coefficient_t p = 3, q = n / p, r = n % p; p <= q; p += 2, q = n / p, r = n % p)
            if (!r) return false;
        return true;
    }

    std::vector<coefficient_t> MultiplicativeInverseVector(const coefficient_t m) {
        std::vector<coefficient_t> inverse(m);
        inverse[1] = 1;
        // m = a * (m / a) + m % a
        // Multipying with inverse(a) * inverse(m % a):
        // 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
        for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
        return inverse;
    }

    index_t getNextVertex(index_t& v, const index_t idx, const index_t k,
        const CBinomialCoeffTable& binomial_coeff) {
        if (binomial_coeff(v, k) > idx) {
            index_t count = v;
            while (count > 0) {
                index_t i = v;
                index_t step = count >> 1;
                i -= step;
                if (binomial_coeff(i, k) > idx) {
                    v = --i;
                    count -= step + 1;
                }
                else
                    count = step;
            }
        }
        assert(binomial_coeff(v, k) <= idx);
        assert(binomial_coeff(v + 1, k) > idx);
        return v;
    }

    std::vector<index_t> VerticesOfSimplex(const index_t simplex_index, const index_t dim, const index_t n,
        const CBinomialCoeffTable& binomial_coeff) {
        std::vector<index_t> vertices;
        getSimplexVertices(simplex_index, dim, n, binomial_coeff, std::back_inserter(vertices));
        return vertices;
    }

    const entry_t& getEntry(const entry_t& e) { return e; }
    value_t getDiameter(CDiameterIndexT i) { return i.first; };
    index_t getIndex(CDiameterIndexT i) { return i.second; };

    const entry_t& getEntry(const diameter_entry_t& p) { return p.second; }
    entry_t& getEntry(diameter_entry_t& p) { return p.second; }
    const index_t getIndex(const diameter_entry_t& p) { return getIndex(getEntry(p)); }
    const coefficient_t getCoefficient(const diameter_entry_t& p) { return getCoefficient(getEntry(p)); }
    const value_t& getDiameter(const diameter_entry_t& p) { return p.first; }
    void setCoefficient(diameter_entry_t& p, const coefficient_t c) { setCoefficient(getEntry(p), c); }

}