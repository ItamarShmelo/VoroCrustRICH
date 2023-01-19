
#include "trees.hpp"
#include <boost/random.hpp>
#include <algorithm>

Trees::Trees(): kd_vertices(nullptr), 
                kd_edges(nullptr), 
                kd_faces(nullptr),
                vertices_points(),
                edges_points(),
                faces_points() {}
