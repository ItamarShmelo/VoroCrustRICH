#ifndef TREES
#define TREES

#include "../../source/ANN/ANN.h"
#include "../../source/ANN/kd_tree.h"
#include "PLC/PL_Complex.hpp"
#include <memory>

class Trees {

    public:
    
        std::shared_ptr<ANNkd_tree> kd_vertices;
        std::shared_ptr<ANNkd_tree> kd_edges;
        std::shared_ptr<ANNkd_tree> kd_faces;

        ANNpointArray vertices_points;
        ANNpointArray edges_points;
        ANNpointArray faces_points;

        Trees();

        void loadPLC(PL_Complex const& plc, std::size_t const Nsample_edges, std::size_t const Nsample_faces);

        ANNpointArray pointsFromVertices(std::vector<Vertex> const& vertices);

        ANNpointArray superSampleEdges(std::vector<Edge> const& edges, std::size_t const Nsample);
        
        ANNpointArray superSampleFaces(std::vector<Face> const& faces, std::size_t const Nsample);

};

#endif // TREES
