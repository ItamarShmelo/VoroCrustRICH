#ifndef _BRUTE_FORCE_RANGE_HPP
#define _BRUTE_FORCE_RANGE_HPP

#include "RangeFinder.hpp"

class BruteForceFinder : public RangeFinder
{
public:
    template<typename RandomAccessIterator>
    BruteForceFinder(const RandomAccessIterator &first, const RandomAccessIterator &last):
                points(std::vector<Vector3D>()){this->points.insert(this->points.begin(), first, last);};
    BruteForceFinder(const std::vector<Vector3D> &points): BruteForceFinder(points.begin(), points.end()){};

    inline std::vector<Vector3D> range(const Vector3D &center, double radius) const override
    {
        std::vector<Vector3D> result;
        for(const Vector3D &point : this->points)
        {
            double distanceSquared = (point.x - center.x) * (point.x - center.x) + (point.y - center.y) * (point.y - center.y) + (point.z - center.z) * (point.y - center.y);
            if(distanceSquared <= radius * radius)
            {
                result.push_back(point);
            }
        }
        return result;
    }

    inline size_t size() const override{return this->points.size();};

private:
    std::vector<Vector3D> points;
};

#endif // _BRUTE_FORCE_RANGE_HPP