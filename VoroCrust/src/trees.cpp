
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

ANNpointArray Trees::superSampleEdges(std::vector<Edge> const& edges, std::size_t const Nsample){
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937, boost::uniform_01<>> rand_gen(rng, zeroone);

    std::vector<double> start_len(edges.size(), 0.0);

    double total_len = 0;
    for(unsigned int i=0; i<edges.size(); ++i){
        Edge const& edge = edges[i];
        start_len[i] = total_len;

        double const edge_len = abs(edge->vertex2->vertex - edge->vertex1->vertex);
        total_len += edge_len;
    }

    std::cout << "\nEdge Samples : \n---------------------\n";
    ANNpointArray points = annAllocPts(Nsample, 3);

    for(std::size_t i=0; i<Nsample; ++i){
        double const sample = rand_gen()*total_len;

        auto const iter_lower_bound = std::lower_bound(start_len.begin(), start_len.end(), sample);
        std::size_t edge_index = std::distance(start_len.begin(), iter_lower_bound) - 1;

        if(std::abs(start_len[edge_index] - sample) < 1e-14){
            //! TODO: make eps a user given parameter. 
            i--;
            continue;
        }

        Edge const& edge = edges[edge_index];
        Vector3D const& edge_vec = (edge->vertex2->vertex - edge->vertex1->vertex);

        double const factor = (sample - start_len[edge_index]) / abs(edge_vec);
        Vector3D const& point = edge->vertex1->vertex + factor*edge_vec;

        points[i][0] = point.x;
        points[i][1] = point.y;
        points[i][2] = point.z;        
    }

    return points;
}
