#ifndef _RANGE_TREE_FINDER_HPP
#define _RANGE_TREE_FINDER_HPP

#include <iostream> // todo remove

#include "ds/BinaryTree.h"
#include "ds/RangeTree/RangeTree.h"

// todo, unnecessary?
template class BinaryTree<Vector3D>;
template class RangeTree<Vector3D>;

#include "RangeFinder.hpp"

class RangeTreeFinder : public RangeFinder
{
public:
    template<typename RandomAccessIterator>
    RangeTreeFinder(RandomAccessIterator first, RandomAccessIterator last);
    inline RangeTreeFinder(std::vector<Vector3D> &myPoints): RangeTreeFinder(myPoints.begin(), myPoints.end()){};
    ~RangeTreeFinder();
    inline std::vector<size_t> range(const Vector3D &center, double radius) const override{
        std::cerr << "Warning: range queries in trees don't work!" << std::endl;
        // return this->rangeTree->circularRange(center, radius);
        return std::vector<size_t>();
    };
    inline const Vector3D &getPoint(size_t index) const override{
        std::cerr << "Warning: range queries in trees don't work!" << std::endl;
        // return this->rangeTree->circularRange(center, radius);
        return Vector3D(0, 0, 0);
    }
    inline size_t size() const override{return this->rangeTree->size();};

private:
    RangeTree<Vector3D> *rangeTree;
};

template<typename RandomAccessIterator>
RangeTreeFinder::RangeTreeFinder(RandomAccessIterator first, RandomAccessIterator last)
{
    this->rangeTree = new RangeTree<Vector3D>(3);
    this->rangeTree->build(first, last);
}

RangeTreeFinder::~RangeTreeFinder()
{
    delete this->rangeTree;
}

#endif // _RANGE_TREE_FINDER_HPP