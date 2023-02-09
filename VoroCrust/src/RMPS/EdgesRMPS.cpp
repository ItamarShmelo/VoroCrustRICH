#include "EdgesRMPS.hpp"

EdgesRMPS::EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_, double const sharpTheta_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), alpha(alpha_), sharpTheta(sharpTheta_), uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()) {}

void EdgesRMPS::loadEdges(std::vector<Edge> const& sharp_edges){
    for(Edge const& edge : sharp_edges)
        eligble_edges.push_back(EligbleEdge(edge->vertex1->vertex, edge->vertex2->vertex, edge->crease_index));
}

std::pair<double const, std::vector<double> const> EdgesRMPS::calculateTotalLengthAndStartLengthOfEligbleEdges(){
    std::vector<double> start_len(eligble_edges.size(), 0.0);
    double total_len = 0.0;

    std::size_t num_of_eligble_edges = eligble_edges.size();
    for(std::size_t i=0; i<num_of_eligble_edges; ++i){
        EligbleEdge const& edge = eligble_edges[i];
        start_len[i] = total_len;
        double const edge_len = abs(edge[1] - edge[0]);
        total_len += edge_len;
    }    

    return std::pair<double const, std::vector<double> const> (total_len, start_len);
}

void EdgesRMPS::divideEligbleEdges(){
    std::vector<EligbleEdge> new_eligble_edges(2*eligble_edges.size(), EligbleEdge(Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 0.0), 0));

    std::size_t const num_of_eligble_edges = eligble_edges.size();
    for(std::size_t i=0; i<num_of_eligble_edges; ++i){
        EligbleEdge const& edge = eligble_edges[i];
        Vector3D const& midpoint = 0.5*(edge[1] + edge[0]);

        new_eligble_edges[2*i][0] = edge[0];
        new_eligble_edges[2*i][1] = midpoint;
        new_eligble_edges[2*i].crease_index = edge.crease_index;

        new_eligble_edges[2*i+1][0] = midpoint;
        new_eligble_edges[2*i+1][1] = edge[1];
        new_eligble_edges[2*i+1].crease_index = edge.crease_index;
    }

    eligble_edges = new_eligble_edges;
}

bool EdgesRMPS::checkIfPointIsDeeplyCovered(Vector3D const& p, VoroCrust_KD_Tree_Ball const& edges_ball_tree) const {
    std::size_t const nn_index = edges_ball_tree.nearestNeighbor(p);

    Vector3D const& q = edges_ball_tree.points[nn_index];
    double const r_q = edges_ball_tree.ball_radii[nn_index];

    double const r_max = (r_q + L_Lipschitz*distance(p, q)) / (1.0 - L_Lipschitz); 
    
    std::vector<int> suspects = edges_ball_tree.radiusSearch(p, r_max);
    
    for(int const i : suspects){
        Vector3D const& center = edges_ball_tree.points[i];
        double const r = edges_ball_tree.ball_radii[i];

        double const dist = distance(p, center);
        if(dist <= r*(1.0 - alpha)){
            return true;
        }
    }

    return false;
}

std::pair<bool, Vector3D> EdgesRMPS::sampleEligbleEdges(double const total_len, std::vector<double> const& start_len) {
    double const sample = uni01_gen()*total_len;

    auto const iter_lower_bound = std::lower_bound(start_len.begin(), start_len.end(), sample);
    std::size_t edge_index = std::distance(start_len.begin(), iter_lower_bound) - 1;
    
    if(std::abs(start_len[edge_index] - sample) < 1e-14){
        return std::pair<bool, Vector3D>(false, Vector3D(0.0, 0.0, 0.0));
    }
    
    EligbleEdge const& edge = eligble_edges[edge_index];
    Vector3D const& edge_vec = (edge[1] - edge[0]);
    double const factor = (sample - start_len[edge_index]) / abs(edge_vec);
    Vector3D const& point = edge[0] + factor*edge_vec;

    return std::pair<bool, Vector3D>(true, point);
}

void EdgesRMPS::discardEligbleEdgesContainedInCornerBalls(VoroCrust_KD_Tree_Ball const& corners_ball_tree){
    std::vector<std::size_t> to_discard;

    std::size_t const num_of_eligble_edges = eligble_edges.size();
    for(std::size_t i=0 ; i<num_of_eligble_edges; ++i){
        EligbleEdge const& edge = eligble_edges[i];
        std::size_t const nn_index = corners_ball_tree.nearestNeighborToSegment(edge.edge);

        Vector3D const& center = corners_ball_tree.points[nn_index];
        double const r = corners_ball_tree.ball_radii[nn_index];

        if(distance(edge[0], center) <= r && distance(edge[1], center) <= r){
            to_discard.push_back(i);
        }
    }

    for(std::size_t i=to_discard.size()-1; i>=0; --i){
        std::size_t const ind_to_discard = to_discard[i];
        eligble_edges.erase(eligble_edges.begin() + ind_to_discard);
    }
}
