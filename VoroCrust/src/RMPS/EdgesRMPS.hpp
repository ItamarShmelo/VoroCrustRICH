#ifndef EDGES_RMPS
#define EDGES_RMPS

#include "../PLC/PL_Complex.hpp"
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"
#include "../VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"
#include "../trees.hpp"
#include <vector>
#include <array>
#include <boost/random.hpp>

struct EligbleEdge {
    std::array<Vector3D, 2> edge;
    std::size_t crease_index;

    EligbleEdge() : edge(), crease_index(0) {}
    EligbleEdge(Vector3D const& v1, Vector3D const& v2, std::size_t const crease_index_) : edge({v1, v2}), crease_index(crease_index_) {}

    Vector3D const& operator [] (std::size_t const index) const {
        return edge[index];
    }

    Vector3D& operator [] (std::size_t const index)  {
        return edge[index];
    }
};

class EdgesRMPS {
    public:
        double const maxRadius;
        double const L_Lipschitz;
        double const alpha;
        double const sharpTheta;
        double const rejection_probability = 0.1;

        boost::variate_generator<boost::mt19937, boost::uniform_01<>> uni01_gen;
        
        std::vector<EligbleEdge> eligble_edges;

        EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_, double const sharpTheta_);
        ~EdgesRMPS() = default;

        //! \brief loads the sharp edges vector to the eligble edge vector
        void loadEdges(std::vector<Edge> const& sharp_edges);

        //! \brief do the RMPS sampling stage until there are no more eligble edges 
        void doSampling(VoroCrust_KD_Tree_Ball &edges_ball_tree, Trees const& trees);

    private:

        /*! \brief calculates the current total length of the eligble edges and the start len for each edge
            \return returns a pair <total_length, start_length> where start_length is a vector with the same size as the eligble edges
        */
        std::pair<double const, std::vector<double> const> calculateTotalLengthAndStartLengthOfEligbleEdges();

        //! \brief divides the eligble edges to two eligble edges at the midpoint
        void divideEligbleEdges();

        //! \brief checks if `p` is deeply covered by any ball in edges_ball_tree
        bool checkIfPointIsDeeplyCovered(Vector3D const& p, VoroCrust_KD_Tree_Ball const& edges_ball_tree, VoroCrust_KD_Tree_Ball const& corners_ball_tree) const;

        /*! \brief sample the eligble edges
            \return <success, sample> where success is a bool flag and is false if sampling failed
        */
        std::tuple<bool, std::size_t const, Vector3D const> sampleEligbleEdges(double const total_len, std::vector<double> const& start_len);

        //! \brief discard any eligble edges that meet the critrea for discardtion 
        void discardEligbleEdges(VoroCrust_KD_Tree_Ball &edges_ball_tree, Trees const& trees);

        //! \brief discard any eligble edge fully contained inside a corner ball
        void discardEligbleEdgesContainedInCornerBalls(VoroCrust_KD_Tree_Ball const& corners_ball_tree);

        //! \brief enfoce that the ball is co smooth with with every point on the crease that is in its radius.
        void enforceEdgeBallCoSmoothWithEligbleEdge(VoroCrust_KD_Tree_Boundary const& edges_boundary_tree, VoroCrust_KD_Tree_Ball &edges_ball_tree, std::size_t const ball_index);

        //
        double calculateSmoothnessLimitation(Vector3D const& center, Vector3D const& parallel, std::size_t const feature_index, VoroCrust_KD_Tree_Boundary const& edges_boundary_tree) const;

        bool isEligbleEdgeIsDeeplyCoveredInEdgeBall(EligbleEdge const& edge, VoroCrust_KD_Tree_Ball const& edges_ball_tree, std::size_t const ball_index);

        double calculateInitialRadius(Vector3D const& point, std::size_t const edge_index, VoroCrust_KD_Tree_Ball const& edges_ball_tree, VoroCrust_KD_Tree_Boundary const& edges_boundary_tree) const;

        bool isBallContainesExistingSample(Vector3D const& center, double const radius, VoroCrust_KD_Tree_Ball const& edges_ball_tree) const;
};

#endif // EDGES_RMPS