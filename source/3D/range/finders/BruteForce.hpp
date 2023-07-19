#ifndef _BRUTE_FORCE_RANGE_HPP
#define _BRUTE_FORCE_RANGE_HPP

#include <iostream> // todo remove

#include "RangeFinder.hpp"

class BruteForceFinder : public RangeFinder
{
public:
    template<typename RandomAccessIterator>
    BruteForceFinder(const RandomAccessIterator &first, const RandomAccessIterator &last):
                points(std::vector<Vector3D>()){this->points.insert(this->points.begin(), first, last); this->pointsSize = this->points.size();};
    BruteForceFinder(const std::vector<Vector3D> &points): BruteForceFinder(points.begin(), points.end()){};

    inline std::vector<size_t> range(const Vector3D &center, double radius) const override
    {
        std::vector<size_t> result;
        const Vector3D *_points = this->points.data();
	__builtin_prefetch(_points);
	for(size_t i = 0; i < this->pointsSize; i++)
        {
            const Vector3D &point = _points[i];
            double distanceSquared = (point.x - center.x) * (point.x - center.x) + (point.y - center.y) * (point.y - center.y) + (point.z - center.z) * (point.z - center.z);
            if(distanceSquared <= radius * radius)
            {
                result.push_back(i);
            }

        }
        return result;
    }

    inline const Vector3D &getPoint(size_t index) const override{return this->points[index];};

    inline size_t size() const override{return this->pointsSize;};

private:
    std::vector<Vector3D> points;
    size_t pointsSize;
};

#endif // _BRUTE_FORCE_RANGE_HPP
