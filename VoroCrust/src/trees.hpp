#ifndef TREES
#define TREES

#include "../../source/ANN/ANN.h"
#include "../../source/ANN/kd_tree.h"
#include "kdtree/kdtree.h"
#include "PLC/PL_Complex.hpp"
#include <memory>

struct kdtree_deleter{
    void operator ()(kdtree* p){
        kd_free(p);
    }
};

class Trees {

    public:
        //! \brief kd_tree holding information on vertices
        std::shared_ptr<ANNkd_tree> kd_vertices;
        //! \brief kd_tree holding information on edges
        std::shared_ptr<ANNkd_tree> kd_edges;
        //! \brief kd_tree holding information on vertices
        std::shared_ptr<ANNkd_tree> kd_faces;

        //! IMPORTANT: point arrays need to be saved since the kd_trees holds a pointer not a copy!
        //! \brief point array used to define `kd_vertices`
        ANNpointArray vertices_points;
        //! \brief point array used to define `kd_edges`
        ANNpointArray edges_points;
        //! \brief point array used to define `kd_faces`
        ANNpointArray faces_points;

        std::shared_ptr<kdtree> ball_kd_vertices;
        std::shared_ptr<kdtree> ball_kd_edges;
        std::shared_ptr<kdtree> ball_kd_faces;

        Trees();
        ~Trees() = default;

        //! \brief initialize the trees from a given PL_Complex
        //! \param plc
        //! \param Nsample_edges number of samples for the super sampling of plc.edges
        //! \param Nsample_faces number of samples for the super sampling of plc.faces
        void loadPLC(PL_Complex const& plc, std::size_t const Nsample_edges, std::size_t const Nsample_faces);

        //! \brief create a point arrey from a given vector of vertices
        //! \param vertices vector of Vertex to be turned into a point array
        ANNpointArray pointsFromVertices(std::vector<Vertex> const& vertices);

        //! \brief creates a point array by super sampling edges in a given edges vector
        //! \param edges vector of Eace to be super sampled for points
        //! \param Nsample number of points to sample
        ANNpointArray superSampleEdges(std::vector<Edge> const& edges, std::size_t const Nsample);
        
        //! \brief creates a point array by super sampling faces
        //! \param faces vector of Face to be super sampled for points
        //! \param Nsample number of points to sample
        ANNpointArray superSampleFaces(std::vector<Face> const& faces, std::size_t const Nsample);

};

#endif // TREES
