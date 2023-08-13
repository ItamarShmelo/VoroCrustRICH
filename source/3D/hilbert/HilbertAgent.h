/**
 * The hilbert agent is a class that is responsible for rearranging points and
 * getting information about the hilbert infastructure of the space.
*/

#ifndef _RICH_3DINTERSECT_H
#define _RICH_3DINTERSECT_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <boost/container/flat_set.hpp>
#include <mpi.h>
#include "utils/balance/balance.hpp"
#include "hilbertTypes.h"
#include "HilbertOrder3D.hpp"

#define POINT_SEND_TAG 115
#define FINISH_TAG 116
#define POINTS_EXCHANGE_RECEIVE_CYCLE 5

#define AVERAGE_INTERSECT 128

struct _Hilbert3DPoint
{
    coord_t x, y, z;
    hilbert_index_t idx;

    inline _Hilbert3DPoint(coord_t x, coord_t y, coord_t z, hilbert_index_t idx): x(x), y(y), z(z), idx(idx){};
    inline explicit _Hilbert3DPoint(): _Hilbert3DPoint(0, 0, 0, -1){};
    bool operator<(const _Hilbert3DPoint &other) const{return this->idx < other.idx;};
    bool operator==(const _Hilbert3DPoint &other) const{return this->idx == other.idx;};
};

class HilbertAgent
{
    template<typename T>
    using _set = boost::container::flat_set<T>;

public:
    HilbertAgent(const Vector3D &origin, const Vector3D &corner);
    HilbertAgent(const Vector3D &origin, const Vector3D &corner, int order): HilbertAgent(origin, corner){this->setOrder(order);};
    _set<size_t> getIntersectingCircle(const Vector3D &center, coord_t r) const;
    inline int getOrder() const{return this->order;};
    std::pair<Vector3D, Vector3D> getBoundingBox() const{return std::make_pair(this->myll, this->myur);};
    hilbert_index_t xyz2d(const Vector3D &point) const;
    Vector3D d2xyz(hilbert_index_t d) const;

    inline int getCellOwner(hilbert_index_t d) const
    {
        return std::min<size_t>(std::upper_bound(this->range.begin(), this->range.end(), d) - this->range.begin(), this->size - 1);
    };

    inline hilbert_index_t getMyHilbertMin() const{return this->myHilbertMin;};
    inline hilbert_index_t getMyHilbertMax() const{return this->myHilbertMax;};
    void determineBorders(const std::vector<Vector3D> &points);
    std::vector<Vector3D> determineBordersAndExchange(const std::vector<Vector3D> &points, std::vector<size_t> &self_index_, std::vector<int> &sentprocs_, std::vector<std::vector<size_t>> &sentpoints_);
    std::vector<Vector3D> pointsExchange(const std::vector<Vector3D> &points, std::vector<size_t> &self_index_, std::vector<int> &sentprocs_, std::vector<std::vector<size_t>> &sentpoints_, std::vector<double> &radiuses) const;
    void calculateBoundingBox();
    void buildHilbertOctTree();
    void setOrder(int order);

private:
    HilbertCurve3D curve;
    Vector3D ll, ur, myll, myur, dx, sidesLengths;
    int order;
    int rank, size;
    int hilbert_cells;
    int pointsPerRank;
    hilbert_index_t myHilbertMin, myHilbertMax;
    std::vector<hilbert_index_t> range;
    
    void pointsReceive(std::vector<Vector3D> &points, std::vector<double> &radiuses, bool blocking) const;
    inline int getOwner(const Vector3D &point) const{return this->getCellOwner(this->xyz2d(point));};
};

#endif // _RICH_3DINTERSECT_H
