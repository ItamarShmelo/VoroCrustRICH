#ifndef _KDTREE_FINDER_HPP
#define _KDTREE_FINDER_HPP

#include "ds/KDTree/KDTree.hpp"
#include "utils/IndexedVector.hpp"
#include "RangeFinder.hpp"

#define DIMENSIONS 3

class KDTreeFinder : public RangeFinder
{
public:
    template<typename RandomAccessIterator>
    KDTreeFinder(RandomAccessIterator first, RandomAccessIterator last, const Vector3D &ll ,const Vector3D &ur);
    inline KDTreeFinder(std::vector<Vector3D> &myPoints, const Vector3D &ll ,const Vector3D &ur): KDTreeFinder(myPoints.begin(), myPoints.end(), ll, ur){};
    ~KDTreeFinder();
    inline std::vector<size_t> range(const Vector3D &center, double radius) const override{
        std::vector<size_t> toReturn;
        for(const IndexedVector3D &vec : this->kdTree->range(_Sphere<IndexedVector3D>(IndexedVector3D(center.x, center.y, center.z, ILLEGAL_IDX), radius)))
        {
            toReturn.push_back(vec.index);
        }
        return toReturn;
    };
    inline size_t size() const override{return this->myPoints.size();};
    inline const Vector3D &getPoint(size_t index) const override{return this->myPoints[index];};

private:
    KDTree<IndexedVector3D, DIMENSIONS> *kdTree;
    std::vector<Vector3D> myPoints;
};

template<typename RandomAccessIterator>
KDTreeFinder::KDTreeFinder(RandomAccessIterator first, RandomAccessIterator last, const Vector3D &ll ,const Vector3D &ur)
{
    size_t index = 0;
    this->kdTree = new KDTree<IndexedVector3D, DIMENSIONS>(ll, ur);

    for(RandomAccessIterator it = first; it != last; it++)
    {
        const Vector3D &vec = *it;
        IndexedVector3D idx_vec = IndexedVector3D(vec.x, vec.y, vec.z, index);
        this->kdTree->insert(idx_vec);
        this->myPoints.push_back(vec);
        index++;
    }
}

KDTreeFinder::~KDTreeFinder()
{
    delete this->kdTree;
}

#endif // _KDTREE_FINDER_HPP