#ifndef TREES
#define TREES

#include "VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"
#include "PLC/PL_Complex.hpp"
#include <memory>

class Trees {

    public:
        //! \brief kd_tree holding information on sharp corners
        VoroCrust_KD_Tree_Boundary VC_kd_sharp_corners;
        //! \brief kd_tree holding information on sharp edges
        VoroCrust_KD_Tree_Boundary VC_kd_sharp_edges;
        //! \brief kd_tree holding information on surface patches
        VoroCrust_KD_Tree_Boundary VC_kd_faces;

        //! \brief ball kd_tree holding the information of the balls on the sharp corners
        VoroCrust_KD_Tree_Ball ball_kd_vertices;
        //! \brief ball kd_tree holding the information of the balls on the sharp edges
        VoroCrust_KD_Tree_Ball ball_kd_edges;
        //! \brief ball kd_tree holding the information of the balls on the faces.
        VoroCrust_KD_Tree_Ball ball_kd_faces;

        Trees();
        ~Trees() = default;

        //! \brief initialize the trees from a given PL_Complex
        //! \param plc
        //! \param Nsample_edges number of samples for the super sampling of plc.edges
        //! \param Nsample_faces number of samples for the super sampling of plc.faces
        void loadPLC(PL_Complex const& plc, std::size_t const Nsample_edges, std::size_t const Nsample_faces);

        //! \brief create a point arrey from a given vector of vertices
        //! \param vertices vector of Vertex to be turned into a point array
        std::tuple<std::vector<Vector3D>, std::vector<std::size_t>> pointsFromVertices(std::vector<Vertex> const& vertices);

        //! \brief creates a point array by super sampling edges in a given edges vector
        //! \param edges vector of Eace to be super sampled for points
        //! \param Nsample number of points to sample
        std::tuple<std::vector<Vector3D>, std::vector<Vector3D>, std::vector<std::size_t>, std::vector<std::size_t>> superSampleEdges(std::vector<Edge> const& edges, std::size_t const Nsample);
        
        //! \brief creates a point array by super sampling faces
        //! \param faces vector of Face to be super sampled for points
        //! \param Nsample number of points to sample
        std::tuple<std::vector<Vector3D>, std::vector<Vector3D>, std::vector<std::size_t>, std::vector<std::size_t>> superSampleFaces(std::vector<Face> const& faces, std::size_t const Nsample);
};

#endif // TREES
