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
    this->range.resize(this->size);
    this->hilbert_cells = pow(pow(2, order), 3);
    this->pointsPerRank = hilbert_cells / this->size;
    this->sidesLengths = this->dx / pow(2, this->order);
    /*
    this->myHilbertMin = this->rank * this->pointsPerRank;
    this->myHilbertMax = (this->rank == this->size - 1)? this->hilbert_cells - 1 : (this->rank + 1) * this->pointsPerRank;
    */

    // this->calculateBoundingBox();
}

hilbert_index_t HilbertAgent::xyz2d(const Vector3D &point) const
{
    // todo: re-implement
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

void HilbertAgent::calculateBoundingBox()
{
    this->myll = this->ur;
    this->myur = this->ll;
    
    for(hilbert_index_t d = this->myHilbertMin; d <= this->myHilbertMax; d++)
    {
        Vector3D point = this->d2xyz(d);
        this->myll.x = std::min(this->myll.x, point.x);
        this->myll.y = std::min(this->myll.y, point.y);
        this->myll.z = std::min(this->myll.z, point.z);
        this->myur.x = std::max(this->myur.x, point.x);
        this->myur.y = std::max(this->myur.y, point.y);
        this->myur.z = std::max(this->myur.z, point.z);
    }
}

/*
void HilbertAgent::pointsReceive(std::vector<Vector3D> &points, bool blocking) const
{
    MPI_Status status;
    int received = 0;
    int finished_ranks = 0;
    bool sentFinished = false;

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
        if(status.MPI_TAG == POINT_SEND_TAG)
        {
            _3DPoint point;
            MPI_Recv(&point, sizeof(_3DPoint), MPI_BYTE, status.MPI_SOURCE, POINT_SEND_TAG, MPI_COMM_WORLD, &status);
            points.push_back(Vector3D(point.x, point.y, point.z));
        } else if(status.MPI_TAG == FINISH_TAG)
        {
            int dummy;
            MPI_Recv(&dummy, 1, MPI_BYTE, status.MPI_SOURCE, FINISH_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            finished_ranks++;
        } else
        {
            std::cerr << "Invalid tag arrived (tag " << status.MPI_TAG << ", source " << status.MPI_SOURCE << ", dest " << this->rank << ")" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &received, &status);

        if(blocking)
        {
            if(!received and !sentFinished)
            {
                // we are done! send a finish
                for(int i = 0; i < this->size; i++)
                {
                    int dummy = 0;
                    MPI_Send(&dummy, 1, MPI_BYTE, i, FINISH_TAG, MPI_COMM_WORLD);
                }
                sentFinished = true;
            }
            if(finished_ranks != this->size)
            {
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            }
        }
    }
}

std::vector<Vector3D> HilbertAgent::pointsExchange(const std::vector<Vector3D> &points, std::vector<size_t> &self_index_, std::vector<int> &sentprocs_, std::vector<std::vector<size_t>> &sentpoints_) const
{
    std::vector<Vector3D> new_points;

    for(size_t i = 0; i < points.size(); i++)
    {
        int node = this->getOwner(points[i]);
        _3DPoint point = {points[i].x, points[i].y, points[i].z};
        if(node != this->rank)
        {
            // point is not mine
            size_t index = std::find(sentprocs_.begin(), sentprocs_.end(), node) - sentprocs_.begin();
            if(index == sentprocs_.size())
            {
                sentprocs_.push_back(node);
                sentpoints_.push_back(std::vector<size_t>());
            }
            sentpoints_[index].push_back(i);
            MPI_Send(&point, sizeof(_3DPoint), MPI_BYTE, node, POINT_SEND_TAG, MPI_COMM_WORLD);
        }
        else
        {
            new_points.push_back(points[i]);
            self_index_.push_back(i);
        }

        if(i % POINTS_EXCHANGE_RECEIVE_CYCLE == 0)
        {
            this->pointsReceive(new_points, false);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    this->pointsReceive(new_points, true);
    return new_points;
}

*/

struct _Hilbert3DPoint
{
    coord_t x, y, z;
    hilbert_index_t idx;

    inline _Hilbert3DPoint(coord_t x, coord_t y, coord_t z, hilbert_index_t idx): x(x), y(y), z(z), idx(idx){};
    inline explicit _Hilbert3DPoint(): _Hilbert3DPoint(0, 0, 0, -1){};
    bool operator<(const _Hilbert3DPoint &other) const{return this->idx < other.idx;};
    bool operator==(const _Hilbert3DPoint &other) const{return this->idx == other.idx;};
};

std::vector<Vector3D> HilbertAgent::determineBordersAndExchange(const std::vector<Vector3D> &points)
{
    // todo: rename
    std::vector<_Hilbert3DPoint> indices;
    for(const Vector3D &point : points)
    {
        indices.push_back(_Hilbert3DPoint(point.x, point.y, point.z, this->xyz2d(point)));
    }
    BalanceJob<_Hilbert3DPoint> balance(indices);
    indices = balance.balance();

    if(!indices.empty())
    {
        this->myHilbertMin = indices[0].idx;
        this->myHilbertMax = indices[indices.size() - 1].idx;
    }
    else
    {
        this->myHilbertMin = this->myHilbertMax = 0;
    }
    if(this->myHilbertMax == static_cast<size_t>(this->hilbert_cells - 1))
    {
        this->myHilbertMax = this->hilbert_cells;
    }
    MPI_Allgather(&this->myHilbertMax, sizeof(hilbert_index_t), MPI_BYTE, &this->range[0], sizeof(hilbert_index_t), MPI_BYTE, MPI_COMM_WORLD);
    std::vector<Vector3D> result;
    for(const _Hilbert3DPoint &point : indices)
    {
        result.push_back(Vector3D(point.x, point.y, point.z));
    }
    return result;
}

boost::container::flat_set<size_t> HilbertAgent::getIntersectingCircle(const Vector3D &center, coord_t r) const
{
    boost::container::flat_set<size_t> hilbertCells;
    hilbertCells.reserve(AVERAGE_INTERSECT);

    coord_t _minX, _maxX;
    _minX = std::max(std::min(center.x - r, this->ur.x), this->ll.x);
    _maxX = std::max(std::min(center.x + r, this->ur.x), this->ll.x);
    _minX = getClosestCornerBelow(_minX, this->ll.x, this->sidesLengths.x);
    _maxX = getClosestCornerAbove(_maxX, this->ll.x, this->sidesLengths.x);

    for(coord_t _x = _minX; _x <= _maxX; _x += this->sidesLengths.x)
    {
        coord_t closestX =  (center.x < _x)? _x : ((center.x > _x + this->sidesLengths.x)? _x + this->sidesLengths.x : center.x);
        coord_t distanceXSquared = r * r - (closestX - center.x) * (closestX - center.x);
        distanceXSquared = std::max(static_cast<double>(0), distanceXSquared);

        coord_t _minY, _maxY;
        _minY = std::max(std::min(center.y - sqrt(distanceXSquared), this->ur.y), this->ll.y);
        _maxY = std::max(std::min(center.y + sqrt(distanceXSquared), this->ur.y), this->ll.y);
        _minY = getClosestCornerBelow(_minY, this->ll.y, this->sidesLengths.y);
        _maxY = getClosestCornerAbove(_maxY, this->ll.y, this->sidesLengths.y);

        for(coord_t _y = _minY; _y <= _maxY; _y += this->sidesLengths.y)
        {
            coord_t closestY =  (center.y < _y)? _y : ((center.y > _y + this->sidesLengths.y)? _y + this->sidesLengths.y : center.y);
            coord_t distanceYSquared = distanceXSquared - (closestY - center.y) * (closestY - center.y);
            distanceYSquared = std::max(static_cast<double>(0), distanceYSquared);
            
            coord_t _minZ, _maxZ;
            _minZ = std::max(std::min(center.z - sqrt(distanceYSquared), this->ur.z), this->ll.z);
            _maxZ = std::max(std::min(center.z + sqrt(distanceYSquared), this->ur.z), this->ll.z);
            _minZ = getClosestCornerBelow(_minZ, this->ll.z, this->sidesLengths.z);
            _maxZ = getClosestCornerAbove(_maxZ, this->ll.z, this->sidesLengths.z);
            
            for(coord_t _z = _minZ; _z <= _maxZ; _z += this->sidesLengths.z)
            {
                coord_t closestZ =  (center.z < _z)? _z : ((center.z > _z + this->sidesLengths.z)? _z + this->sidesLengths.z : center.z);

                if(((closestX - center.x) * (closestX - center.x) + (closestY - center.y) * (closestY - center.y) + (closestZ - center.z) * (closestZ - center.z)) <= r*r + __DBL_EPSILON__)
                {
                    // the testing point is inside the circle iff the whole cube intersects the circle
                    hilbertCells.insert(this->xyz2d(Vector3D(_x + (this->sidesLengths.x) / 2, _y + (this->sidesLengths.y) / 2, _z + (this->sidesLengths.z) / 2)));
                }
            }
        }
    }
    
    return hilbertCells;
}