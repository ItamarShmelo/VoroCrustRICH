#include "EdgesRMPS.hpp"

EdgesRMPS::EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), alpha(alpha_), uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()) {}

void EdgesRMPS::loadEdges(std::vector<Edge> const& sharp_edges){
    for(Edge const& edge : sharp_edges)
        eligble_edges.push_back(EligbleEdge(edge->vertex1->vertex, edge->vertex2->vertex, edge->crease_index));
}

