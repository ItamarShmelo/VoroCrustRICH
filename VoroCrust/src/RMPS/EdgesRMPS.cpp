#include "EdgesRMPS.hpp"

EdgesRMPS::EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), alpha(alpha_), uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()) {}

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
