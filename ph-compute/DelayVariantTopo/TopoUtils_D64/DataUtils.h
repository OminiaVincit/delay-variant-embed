#pragma once
#include <boost\shared_ptr.hpp>

namespace NDataUtils {
    template <class PointT>
    class CPointCloud {
    public:
        CPointCloud() : m_dim(), m_pcl() {};
        CPointCloud(size_t _dim) : m_dim(_dim), m_pcl() {};

        size_t Dim() const { return m_dim; };
        size_t Nums() const { return m_pcl.size(); };
        std::vector<std::vector<PointT>> Pcl() const { return m_pcl; };

        bool AddElement(std::vector<PointT>& elem) {
            if (elem.size() == m_dim) {
                m_pcl.push_back(elem);
            }
            return false;
        }
		bool IsEmpty() {
			return (Nums() == 0);
		}
    private:
        size_t m_dim;
        std::vector<std::vector<PointT>> m_pcl;
    };

    template<class T>
    class CPointCloudPtr {
    public:
        typedef boost::shared_ptr<CPointCloud<T>> Type;
    };
    
}