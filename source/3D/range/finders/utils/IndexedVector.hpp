#ifndef _INDEXED_VECTOR_HPP
#define _INDEXED_VECTOR_HPP

#include "3D/hilbert/hilbertTypes.h"

#define ILLEGAL_IDX -1

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
    inline bool operator<(const IndexedVector3D &other) const{return (*this) <= other;};
    inline coord_t &operator[](size_t idx){switch(idx){case 0: return this->x; case 1: return this->y;} return this->z;};
    inline const coord_t &operator[](size_t idx) const{switch(idx){case 0: return this->x; case 1: return this->y;} return this->z;};
    inline IndexedVector3D operator+(const IndexedVector3D &other) const{return IndexedVector3D(this->x + other.x, this->y + other.y, this->z + other.z, ILLEGAL_IDX);};
    inline IndexedVector3D operator*(coord_t scalar) const{return IndexedVector3D(this->x * scalar, this->y * scalar, this->z * scalar, ILLEGAL_IDX);};
    inline IndexedVector3D operator/(coord_t scalar) const{return this->operator*(1/scalar);};

} IndexedVector3D;


#endif // _INDEXED_VECTOR_HPP
