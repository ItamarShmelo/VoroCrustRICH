#ifndef _RANGE_TREE_FINDER_HPP
#define _RANGE_TREE_FINDER_HPP

#include "ds/BinaryTree.h"
#include "ds/RangeTree/RangeTree.h"
#include "ds/BinaryTree.cpp" // todo: not good
#include "ds/RangeTree/RangeTree.cpp" // todo: not good
#include "utils/IndexedVector.hpp"
#include "RangeFinder.hpp"

// todo, unnecessary?
template class BinaryTree<IndexedVector3D>;
template class RangeTree<IndexedVector3D>;

class RangeTreeFinder : public RangeFinder
{
public:
    template<typename RandomAccessIterator>
    RangeTreeFinder(RandomAccessIterator first, RandomAccessIterator last);
    inline RangeTreeFinder(std::vector<Vector3D> &myPoints): RangeTreeFinder(myPoints.begin(), myPoints.end()){};
    ~RangeTreeFinder();
    inline std::vector<size_t> range(const Vector3D &center, double radius) const override{
        std::vector<size_t> toReturn;
        for(const IndexedVector3D &vec : this->rangeTree->circularRange(center, radius))
        {
            toReturn.push_back(vec.index);
        }
        return toReturn;
    };
    inline size_t size() const override{return this->rangeTree->size();};
    inline const Vector3D &getPoint(size_t index) const override{return this->myPoints[index];};

private:
    RangeTree<IndexedVector3D> *rangeTree;
    std::vector<Vector3D> myPoints;
};

template<typename RandomAccessIterator>
RangeTreeFinder::RangeTreeFinder(RandomAccessIterator first, RandomAccessIterator last)
{
    std::vector<IndexedVector3D> data;
    size_t index = 0;
    for(RandomAccessIterator it = first; it != last; it++)
    {
        const Vector3D &vec = *it;
        data.push_back(IndexedVector3D(vec.x, vec.y, vec.z, index));
        this->myPoints.push_back(vec);
        index++;
    }
    this->rangeTree = new RangeTree<IndexedVector3D>(3);
    this->rangeTree->build(data.begin(), data.end());
}

RangeTreeFinder::~RangeTreeFinder()
{
    delete this->rangeTree;
}

#endif // _RANGE_TREE_FINDER_HPP