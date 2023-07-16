#include "HilbertAgent.h"

namespace
{
    /**
     * rounding up to the nearest hilbert corner in a specific axis.
    */
    inline coord_t getClosestCornerAbove(coord_t val, coord_t minCord, coord_t sideLength)
    {
        return ceil((val - minCord) / sideLength) * sideLength + minCord;
    }

    /**
     * rounding down to the nearest hilbert corner in a specific axis.
    */
    inline coord_t getClosestCornerBelow(coord_t val, coord_t minCord, coord_t sideLength)
    {
        return floor((val - minCord) / sideLength) * sideLength + minCord;
    }
}

HilbertAgent::HilbertAgent(const Vector3D &origin, const Vector3D &corner, int order): ll(origin), ur(corner), order(order), dx(corner - origin)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &this->size);
}

hilbert_index_t HilbertAgent::xyz2d(const Vector3D &point) const
{
    return curve.Hilbert3D_xyz2d(Vector3D((point.x - this->ll.x) / this->dx.x, 
                                            (point.y - this->ll.y) / this->dx.y,
                                            (point.z - this->ll.z) / this->dx.z),
                                this->order);
}

Vector3D HilbertAgent::d2xyz(hilbert_index_t d) const
{
    Vector3D scaledPoint = curve.Hilbert3D_d2xyz(d, this->order);
    return Vector3D(scaledPoint.x * (this->dx.x) + this->ll.x,
                    scaledPoint.y * (this->dx.y) + this->ll.y,
                    scaledPoint.z * (this->dx.z) + this->ll.z);
}

int HilbertAgent::getOwner(const Vector3D &point) const
{
    int hilbert_cells = pow(pow(2, order), 3);
    int pointsPerRank = hilbert_cells / this->size;
    hilbert_index_t d = this->xyz2d(point);
    return std::min<int>(this->size - 1, static_cast<int>(d / pointsPerRank));
}

std::pair<Vector3D, Vector3D> HilbertAgent::getBoundingBox() const
{
    int hilbert_cells = pow(pow(2, order), 3);
    int pointsPerRank = hilbert_cells / this->size;
    hilbert_index_t myMin = rank * pointsPerRank;
    hilbert_index_t myMax = (rank == size - 1)? hilbert_cells - 1 : (rank + 1) * pointsPerRank;
    
    Vector3D ll = this->d2xyz(myMin);
    Vector3D ur(ll);
    
    for(hilbert_index_t d = myMin; d <= myMax; d++)
    {
        Vector3D point = this->d2xyz(d);
        ll.x = std::min(ll.x, point.x);
        ll.y = std::min(ll.y, point.y);
        ll.z = std::min(ll.z, point.z);
        ur.x = std::max(ur.x, point.x);
        ur.y = std::max(ur.y, point.y);
        ur.z = std::max(ur.z, point.z);
    }
    return std::make_pair(ll, ur);
}

void HilbertAgent::pointsReceive(std::vector<Vector3D> &points, int &finished_ranks, bool blocking) const
{
    MPI_Status status;
    int received = 0;

    if(blocking and finished_ranks != this->size)
    {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }
    else
    {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &received, &status);
    }

    while((blocking and finished_ranks != this->size) or (!blocking and received))
    {
        switch(status.MPI_TAG)
        {
            case POINT_SEND_TAG:
                _3DPoint point;
                MPI_Recv(&point, sizeof(_3DPoint), MPI_BYTE, status.MPI_SOURCE, POINT_SEND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                points.push_back(Vector3D(point.x, point.y, point.z));
                break;

            case FINISH_TAG:
                int dummy;
                MPI_Recv(&dummy, 1, MPI_BYTE, status.MPI_SOURCE, FINISH_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                finished_ranks++;
                break;
            
            default:
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                break;
        }
        
        if(blocking and finished_ranks != this->size)
        {
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &received, &status);
        }
    }
}

std::vector<Vector3D> HilbertAgent::pointsExchange(const std::vector<Vector3D> &points, std::vector<size_t> &self_index_, std::vector<int> &sentprocs_, std::vector<std::vector<size_t>> &sentpoints_) const
{
    int finished_ranks = 0;
    std::vector<MPI_Request> requests;
    std::vector<Vector3D> new_points;

    for(size_t i = 0; i < points.size(); i++)
    {
        int node = this->getOwner(points[i]);
        _3DPoint point = {points[i].x, points[i].y, points[i].z};
        if(node != this->rank)
        {
            // point is not mine
            auto it = std::find(sentprocs_.begin(), sentprocs_.end(), node);
            int index = it - sentprocs_.begin();
            if(it == sentprocs_.end())
            {
                sentprocs_.push_back(node);
                sentpoints_.push_back(std::vector<size_t>());
            }
            sentpoints_[index].push_back(i);
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Isend(&point, sizeof(_3DPoint), MPI_BYTE, node, POINT_SEND_TAG, MPI_COMM_WORLD, &requests[requests.size() - 1]);
        }
        else
        {
            new_points.push_back(points[i]);
            self_index_.push_back(i);
        }

        if(i % POINTS_EXCHANGE_RECEIVE_CYCLE == 0)
        {
            this->pointsReceive(new_points, finished_ranks, false);
        }
    }
    // send a finish message
    for(int i = 0; i < this->size; i++)
    {
        int dummy = 0;
        MPI_Send(&dummy, 1, MPI_BYTE, i, FINISH_TAG, MPI_COMM_WORLD);
    }
    this->pointsReceive(new_points, finished_ranks, true);
    return new_points;
}

std::set<hilbert_index_t> HilbertAgent::getIntersectingCircle(const Vector3D &center, coord_t r) const
{
    Vector3D sidesLengths = this->dx / pow(2, this->order);
    std::set<hilbert_index_t> hilbertCells;

    coord_t _minX, _maxX;
    _minX = std::max(std::min(center.x - r, this->ur.x), this->ll.x);
    _maxX = std::max(std::min(center.x + r, this->ur.x), this->ll.x);
    _minX = getClosestCornerBelow(_minX, this->ll.x, sidesLengths.x);
    _maxX = getClosestCornerAbove(_maxX, this->ll.x, sidesLengths.x);

    for(coord_t _x = _minX; _x <= _maxX; _x += sidesLengths.x)
    {
        coord_t closestX =  (center.x < _x)? _x : ((center.x > _x + sidesLengths.x)? _x + sidesLengths.x : center.x);
        coord_t distanceXSquared = r * r - (closestX - center.x) * (closestX - center.x);
        distanceXSquared = std::max(static_cast<double>(0), distanceXSquared);

        coord_t _minY, _maxY;
        _minY = std::max(std::min(center.y - sqrt(distanceXSquared), this->ur.y), this->ll.y);
        _maxY = std::max(std::min(center.y + sqrt(distanceXSquared), this->ur.y), this->ll.y);
        _minY = getClosestCornerBelow(_minY, this->ll.y, sidesLengths.y);
        _maxY = getClosestCornerAbove(_maxY, this->ll.y, sidesLengths.y);

        for(coord_t _y = _minY; _y <= _maxY; _y += sidesLengths.y)
        {
            coord_t closestY =  (center.y < _y)? _y : ((center.y > _y + sidesLengths.y)? _y + sidesLengths.y : center.y);
            coord_t distanceYSquared = distanceXSquared - (closestY - center.y) * (closestY - center.y);
            distanceYSquared = std::max(static_cast<double>(0), distanceYSquared);
            
            coord_t _minZ, _maxZ;
            _minZ = std::max(std::min(center.z - sqrt(distanceYSquared), this->ur.z), this->ll.z);
            _maxZ = std::max(std::min(center.z + sqrt(distanceYSquared), this->ur.z), this->ll.z);
            _minZ = getClosestCornerBelow(_minZ, this->ll.z, sidesLengths.z);
            _maxZ = getClosestCornerAbove(_maxZ, this->ll.z, sidesLengths.z);
            
            for(coord_t _z = _minZ; _z <= _maxZ; _z += sidesLengths.z)
            {
                coord_t closestZ =  (center.z < _z)? _z : ((center.z > _z + sidesLengths.z)? _z + sidesLengths.z : center.z);

                if(((closestX - center.x) * (closestX - center.x) + (closestY - center.y) * (closestY - center.y) + (closestZ - center.z) * (closestZ - center.z)) <= r*r + __DBL_EPSILON__)
                {
                    // the testing point is inside the circle iff the whole cube intersects the circle
                    hilbert_index_t d = curve.Hilbert3D_xyz2d(Vector3D((_x - this->ll.x + (sidesLengths.x) / 2) / this->dx.x,
                                                                        (_y - this->ll.y + (sidesLengths.y) / 2) / this->dx.y,
                                                                        (_z - this->ll.z + (sidesLengths.z) / 2) / this->dx.z),
                                                            static_cast<int>(this->order));
                    hilbertCells.insert(d);
                }
            }
        }
    }
    
    return hilbertCells;
}