
#include "trees.hpp"
#include <boost/random.hpp>
#include <algorithm>

Trees::Trees(): kd_vertices(nullptr), 
                kd_edges(nullptr), 
                kd_faces(nullptr),
                vertices_points(),
                edges_points(),
                faces_points() {}

ANNpointArray Trees::pointsFromVertices(std::vector<Vertex> const& vertices){
    std::size_t const Npoints = vertices.size();
    ANNpointArray points = annAllocPts(Npoints, 3);

    for(std::size_t i = 0; i<Npoints; ++i){
        Vector3D vertex = vertices[i]->vertex;
        
        points[i][0] = vertex.x;
        points[i][1] = vertex.y;
        points[i][2] = vertex.z;
    }
    
    return points;
}
