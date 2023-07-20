#ifndef _RANGE_TREE_FINDER_HPP
#define _RANGE_TREE_FINDER_HPP

#include "ds/BinaryTree.h"
#include "ds/RangeTree/RangeTree.h"
#include "ds/BinaryTree.cpp" // todo: not good
#include "ds/RangeTree/RangeTree.cpp" // todo: not good
#include "RangeFinder.hpp"

typedef struct IndexedVector3D
{
    using coord_type = coord_t;

    coord_t x;
    coord_t y;
    coord_t z;
    size_t index;

    IndexedVector3D(coord_t x, coord_t y, coord_t z, size_t index): x(x), y(y), z(z), index(index){};
    IndexedVector3D(coord_t x, coord_t y, coord_t z): IndexedVector3D(x, y, z, 0){};
    explicit IndexedVector3D(): IndexedVector3D(0, 0, 0){};
    IndexedVector3D(const Vector3D &other): IndexedVector3D(other.x, other.y, other.z){};
    IndexedVector3D(const IndexedVector3D &other): IndexedVector3D(other.x, other.y, other.z, other.index){};
    inline IndexedVector3D &operator=(const IndexedVector3D &other){this->x = other.x; this->y = other.y; this->z = other.z; this->index = other.index; return *this;};
    inline bool operator==(const IndexedVector3D &other) const{return this->x == other.x and this->y == other.y and this->z == other.z and this->index == other.index;};
    inline bool operator<=(const IndexedVector3D &other) const{
        if(this->x < other.x) return true;
        if(this->x == other.x)
        {
            if(this->y < other.y) return true;
            if(this->y == other.y) return (this->z <= other.z);
        }
        return false;
    }
    inline coord_t &operator[](size_t idx){switch(idx){case 0: return this->x; case 1: return this->y;} return this->z;};
    inline const coord_t &operator[](size_t idx) const{switch(idx){case 0: return this->x; case 1: return this->y;} return this->z;};
} IndexedVector3D;


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