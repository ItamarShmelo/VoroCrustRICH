/*

#include <iostream>
#include <sstream>
#include <vector>
#include "OctTree.hpp"
#include "input.hpp"

struct _3DPoint
{
    double x, y, z;

    inline _3DPoint(double x, double y, double z): x(x), y(y), z(z){};
    inline _3DPoint(const _3DPoint &other): _3DPoint(other.x, other.y, other.z){};
    inline _3DPoint(): _3DPoint(0, 0, 0){};
    inline _3DPoint operator+(const _3DPoint &other) const{return _3DPoint(this->x + other.x, this->y + other.y, this->z + other.z);};
    inline _3DPoint operator-(const _3DPoint &other) const{return _3DPoint(this->x - other.x, this->y - other.y, this->z - other.z);};
    inline _3DPoint operator*(double scalar) const{return _3DPoint(this->x * scalar, this->y * scalar, this->z * scalar);};
    inline _3DPoint operator/(double scalar) const{return this->operator*(1/scalar);};
    inline _3DPoint &operator=(const _3DPoint &other){this->x = other.x; this->y = other.y; this->z = other.z; return (*this);};
    inline bool operator==(const _3DPoint &other) const{return this->x == other.x and this->y == other.y and this->z == other.z;};
    inline const double &operator[](size_t idx) const{if(idx == 0) return this->x; if(idx == 1) return this->y; return this->z;};
    inline double &operator[](size_t idx){if(idx == 0) return this->x; if(idx == 1) return this->y; return this->z;};
    friend std::ostream &operator<<(std::ostream &stream, const _3DPoint &point);
    friend std::istream &operator>>(std::istream &stream, _3DPoint &point);
};

std::ostream &operator<<(std::ostream &stream, const _3DPoint &point)
{
    stream << "(" << point.x << ", " << point.y << ", " << point.z << ")";
    return stream;;
}

std::istream &operator>>(std::istream &stream, _3DPoint &point)
{
    char nextChar;
    stream.get(nextChar);
    while(nextChar == ' ' or nextChar == '\n')
    {
        stream.get(nextChar);
    }
    if(nextChar != '(')
    {
        return stream;
    }
    const char delimeter[3] = {',', ',', ')'};

    for(int i = 0; i < 3; i++)
    {
        std::string str("");
        while(true)
        {
            stream.get(nextChar);
            if(nextChar == delimeter[i])
            {
                break;
            }
            str += nextChar;
        }
        point[i] = std::stod(str);
    }
    return stream;;
}

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        std::cerr << "Illegal number of arguments was given" << std::endl;
        return EXIT_FAILURE;
    }
    OctTree<_3DPoint> tree;
    std::vector<_3DPoint> values = readInput<_3DPoint>(argv[1]);

    const size_t goodNumber = 5000;
    std::vector<_3DPoint> goodValues(values.begin(), values.begin() + std::min<size_t>(goodNumber, values.size()));
    std::vector<_3DPoint> badvalues(values.begin() + std::min<size_t>(goodNumber, values.size()), values.end());
    tree.setBounds(_3DPoint(0, 0, 0), _3DPoint(1, 1, 1));

    for(int i = 0; i < goodValues.size(); i++)
    {
        assert(tree.insert(goodValues[i]) == true);
    }
    tree.print();
    for(int i = 0; i < goodValues.size(); i++)
    {
        assert(tree.find(goodValues[i]) == true);
    }
    for(int i = 0; i < badvalues.size(); i++)
    {
        assert(tree.find(badvalues[i]) == false);
    }
    std::cout << "All tests passed!" << std::endl;

    Sphere<_3DPoint> sphere(badvalues[0], (double)0.5, 3);
    for(const _3DPoint &point : tree.range(sphere))
    {
        std::cout << point << std::endl;
    }
    return EXIT_SUCCESS;
}
*/