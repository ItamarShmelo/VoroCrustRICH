/**
 * The hilbert agent is a class that is responsible for rearranging points and
 * getting information about the hilbert infastructure of the space.
*/

#ifndef _RICH_3DINTERSECT_H
#define _RICH_3DINTERSECT_H

// #include "utils/sort.hpp"
#include "hilbertTypes.h"
#include "HilbertOrder3D.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <mpi.h>

#define POINT_SEND_TAG 115
#define FINISH_TAG 116
#define POINTS_EXCHANGE_RECEIVE_CYCLE 5

class HilbertAgent
{
public:
    HilbertAgent(const Vector3D &origin, const Vector3D &corner, int order);
    std::set<hilbert_index_t> getIntersectingCircle(const Vector3D &center, coord_t r) const;
    std::vector<Vector3D> pointsExchange(const std::vector<Vector3D> &points, std::vector<size_t> &self_index_, std::vector<int> &sentprocs_, std::vector<std::vector<size_t>> &sentpoints_) const;
    inline int getOrder() const{return this->order;};
    std::pair<Vector3D, Vector3D> getBoundingBox() const{return std::make_pair(this->myll, this->myur);};
    hilbert_index_t xyz2d(const Vector3D &point) const;
    Vector3D d2xyz(hilbert_index_t d) const;
    inline int getCellOwner(hilbert_index_t d) const{return std::min<int>(this->size - 1, static_cast<int>(d / this->pointsPerRank));};
    inline hilbert_index_t getMyHilbertMin() const{return this->myHilbertMin;};
    inline hilbert_index_t getMyHilbertMax() const{return this->myHilbertMax;};
    void setBorders(const std::vector<Vector3D> &points);

private:
    HilbertCurve3D curve;
    Vector3D ll, ur, myll, myur, dx, sidesLengths;
    int order;
    int rank, size;
    int hilbert_cells;
    int pointsPerRank;
    hilbert_index_t myHilbertMin, myHilbertMax;

    void calculateBoundingBox();
    void pointsReceive(std::vector<Vector3D> &points, bool blocking) const;
    inline int getOwner(const Vector3D &point) const{return this->getCellOwner(this->xyz2d(point));};
};

#endif // _RICH_3DINTERSECT_H
