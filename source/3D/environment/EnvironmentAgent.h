#ifndef ENVIRONMENT_RICH_H
#define ENVIRONMENT_RICH_H

#ifdef RICH_MPI

#include <mpi.h>
#include "3D/hilbert/HilbertOrder3D.hpp"
#include "3D/elementary/Vector3D.hpp"
#include <boost/container/flat_set.hpp>

class EnvironmentAgent
{
public:
    template<typename T>
    using _set = boost::container::flat_set<T>;

    inline EnvironmentAgent(const Vector3D &ll, const Vector3D &ur, const MPI_Comm &comm = MPI_COMM_WORLD): ll(ll), ur(ur), comm(comm)
    {
        MPI_Comm_rank(this->comm, &this->rank);
        MPI_Comm_size(this->comm, &this->size);
    };
    virtual _set<int> getIntersectingRanks(const Vector3D &center, double radius) const = 0;
    virtual int getOwner(const Vector3D &point) const = 0;
    virtual void update(const std::vector<Vector3D> &newPoints) = 0;
    virtual void updateBorders(const std::vector<hilbert_index_t> &newRange, int newOrder) = 0;

    inline hilbert_index_t xyz2d(const Vector3D &point, int order) const
    {
        return EnvironmentAgent::xyz2d(point, this->ll, this->ur, order);
    }

    static inline hilbert_index_t xyz2d(const Vector3D &point, const Vector3D &ll, const Vector3D &ur, int order)
    {
        const Vector3D dx = ur - ll;
        return EnvironmentAgent::curve.Hilbert3D_xyz2d(Vector3D((point.x - ll.x) / dx.x, (point.y - ll.y) / dx.y, (point.z - ll.z) / dx.z), order);
    }

protected:
    static HilbertCurve3D curve;
    Vector3D ll, ur;
    MPI_Comm comm;
    int rank, size;
};

#endif // RICH_MPI

#endif // ENVIRONMENT_RICH_H