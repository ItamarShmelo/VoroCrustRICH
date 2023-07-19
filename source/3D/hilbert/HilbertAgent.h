/**
 * The hilbert agent is a class that is responsible for rearranging points and
 * getting information about the hilbert infastructure of the space.
*/

#ifndef _RICH_3DINTERSECT_H
#define _RICH_3DINTERSECT_H

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
    std::pair<Vector3D, Vector3D> getBoundingBox() const;
    hilbert_index_t xyz2d(const Vector3D &point) const;
    Vector3D d2xyz(hilbert_index_t d) const;

private:
    HilbertCurve3D curve;
    Vector3D ll, ur, dx;
    int order;
    int rank, size;

    void pointsReceive(std::vector<Vector3D> &points, bool blocking) const;
    int getOwner(const Vector3D &point) const;
};

#endif // _RICH_3DINTERSECT_H
