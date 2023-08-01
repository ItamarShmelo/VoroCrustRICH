#ifndef _GEOMETRY_UTILS_RICH_HPP
#define _GEOMETRY_UTILS_RICH_HPP

typedef double coord_t;

template<typename T>
class _BoundingBox
{
public:
    T ll; // leftmost point of the box
    T ur; // rightmost point of the box

    inline _BoundingBox(const T &ll, const T &ur): ll(ll), ur(ur){};
    inline _BoundingBox(): _BoundingBox(T(), T()){};
    
    inline bool contains(const T &point, int dim) const
    {
        for(int i = 0; i < dim; i++)
        {
            if(point[i] < ll[i] | point[i] > ur[i])
            {
                return false;
            }
        }
        return true;
    }
};

template<typename T>
class _Sphere
{
public:
    T center;
    coord_t radius;
    int dim;

    _Sphere(const T &center, coord_t radius, int dim): center(center), radius(radius), dim(dim){};
    inline bool contains(const T &point) const
    {
        coord_t distance = 0;
        for(int i = 0; i < this->dim; i++)
        {
            distance += (point[i] - this->center[i]) * (point[i] - this->center[i]);
        }
        return distance <= (this->radius * this->radius);
    }
};

template<typename T>
bool SphereBoxIntersection(const _BoundingBox<T> &box, const _Sphere<T> &sphere)
{
    T closestPoint;
    coord_t distance = 0;
    for(int i = 0; i < sphere.dim; i++)
    {
        closestPoint[i] = (sphere.center[i] < box.ll[i])? box.ll[i] : ((sphere.center[i] > box.ur[i])? box.ur[i] : sphere.center[i]);
        distance += (closestPoint[i] - sphere.center[i]) * (closestPoint[i] - sphere.center[i]);
    }
    return distance <= sphere.radius * sphere.radius;
};

#endif // _GEOMETRY_UTILS_RICH_HPP