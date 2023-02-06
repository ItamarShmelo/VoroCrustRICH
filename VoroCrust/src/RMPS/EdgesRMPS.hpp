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

        boost::variate_generator<boost::mt19937, boost::uniform_01<>> uni01_gen;
        
        std::vector<EligbleEdge> eligble_edges;

        EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_);
        ~EdgesRMPS() = default;

        void loadEdges(std::vector<Edge> const& sharp_edges);

        void doSampling(VoroCrust_KD_Tree_Ball &edges_ball_tree, Trees const& trees);

    private:

        std::pair<double const, std::vector<double> const> calculateTotalLengthAndStartLengthOfEligbleEdges();

        void divideEligbleEdges();

        bool checkIfPointIsDeeplyCovered(Vector3D const& p, VoroCrust_KD_Tree_Ball const& edges_ball_tree) const;

        std::pair<bool, Vector3D> sampleEligbleEdges(double const total_len, std::vector<double> const& start_len);

        void discardEligbleEdges(Trees const& trees);

        void discardEligbleEdgesContainedInCornerBalls(VoroCrust_KD_Tree_Ball const& corners_ball_tree);

        bool isEdgeBallCoSmoothWithEligbleEdge(EligbleEdge const& edge, std::size_t const ball_index);

        
};

#endif // EDGES_RMPS