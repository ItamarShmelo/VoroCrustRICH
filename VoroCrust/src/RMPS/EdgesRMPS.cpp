#include "EdgesRMPS.hpp"

EdgesRMPS::EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_, double const sharpTheta_, std::shared_ptr<PL_Complex> const& plc_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), alpha(alpha_), sharpTheta(sharpTheta_), uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()), plc(plc_), eligble_edges(), isDeleted() {}

void EdgesRMPS::loadEdges(std::vector<Edge> const& sharp_edges){
    if(not eligble_edges.empty()){
        std::cout << "eligble edges is not empty when loading edges" << std::endl;
        exit(1);
    }

    for(Edge const& edge : sharp_edges)
        eligble_edges.push_back(EligbleEdge(edge->vertex1->vertex, edge->vertex2->vertex, edge->crease_index, edge->index));

    isDeleted = std::vector<bool>(eligble_edges.size(), false);
}

std::pair<double const, std::vector<double> const> EdgesRMPS::calculateTotalLengthAndStartLengthOfEligbleEdges() const {
    std::size_t const num_of_eligble_edges = eligble_edges.size();

    std::vector<double> start_len(num_of_eligble_edges, 0.0);
    double total_len = 0.0;

    for(std::size_t i=0; i<num_of_eligble_edges; ++i){
        EligbleEdge const& edge = eligble_edges[i];
        start_len[i] = total_len;
        double const edge_len = abs(edge[1] - edge[0]);
        total_len += edge_len;
    }    

    return std::pair(total_len, start_len);
}

void EdgesRMPS::divideEligbleEdges(){
    std::size_t const num_of_eligble_edges = eligble_edges.size();

    std::vector<EligbleEdge> new_eligble_edges(2*num_of_eligble_edges, EligbleEdge(Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 0.0), 0, 0));

    for(std::size_t i=0; i<num_of_eligble_edges; ++i){
        EligbleEdge const& edge = eligble_edges[i];
        Vector3D const& midpoint = 0.5*(edge[1] + edge[0]);

        new_eligble_edges[2*i][0] = edge[0];
        new_eligble_edges[2*i][1] = midpoint;
        new_eligble_edges[2*i].crease_index = edge.crease_index;
        new_eligble_edges[2*i].plc_index = edge.plc_index;

        new_eligble_edges[2*i+1][0] = midpoint;
        new_eligble_edges[2*i+1][1] = edge[1];
        new_eligble_edges[2*i+1].crease_index = edge.crease_index;
        new_eligble_edges[2*i+1].plc_index = edge.plc_index;
    }

    eligble_edges = new_eligble_edges;
    isDeleted = std::vector<bool>(eligble_edges.size(), false);
}

bool EdgesRMPS::checkIfPointIsDeeplyCovered(Vector3D const& p, Trees const& trees) const {
    VoroCrust_KD_Tree_Ball const& corners_ball_tree = trees.ball_kd_vertices;
    VoroCrust_KD_Tree_Ball const& edges_ball_tree = trees.ball_kd_edges;

    if(not edges_ball_tree.empty()){
        auto const& [q, r_q] = edges_ball_tree.getBallNearestNeighbor(p);

        // maximal radius for centers which balls can deeply cover p 
        double const r_max = (r_q + L_Lipschitz*distance(p, q)) / (1.0 - L_Lipschitz); 
        
        auto const& suspects = edges_ball_tree.radiusSearch(p, r_max);
        
        for(auto const i : suspects){
            auto const& [center, r] = edges_ball_tree.getBall(i);

            double const dist = distance(p, center);
            if(dist <= r*(1.0 - alpha)){
                return true;
            }
        }
    }
    
    auto const& [center, radius] = corners_ball_tree.getBallNearestNeighbor(p);

    return distance(p, center) <= radius*(1.0-alpha);
}

std::tuple<bool, std::size_t const, Vector3D const> EdgesRMPS::sampleEligbleEdges(double const total_len, std::vector<double> const& start_len) {
    double const sample = uni01_gen()*total_len;

    // find on which edge the sample falls
    auto const& iter_lower_bound = std::lower_bound(start_len.begin(), start_len.end(), sample);
    std::size_t const edge_index = std::distance(start_len.begin(), iter_lower_bound) - 1;
    
    // if sample is too close to a vertex succes = false
    if(std::abs(start_len[edge_index] - sample) < 1e-14){
        return std::tuple<bool, std::size_t const, Vector3D const>(false, 0, Vector3D(0.0, 0.0, 0.0));
    }
    
    // find exact point
    EligbleEdge const& edge = eligble_edges[edge_index];
    Vector3D const& edge_vec = (edge[1] - edge[0]);
    double const factor = (sample - start_len[edge_index]) / abs(edge_vec);
    Vector3D const& point = edge[0] + factor*edge_vec;

    return std::tuple(true, edge_index, point);
}

void EdgesRMPS::discardEligbleEdgesContainedInCornerBalls(Trees const& trees){
    VoroCrust_KD_Tree_Ball const& corners_ball_tree = trees.ball_kd_vertices;

    std::size_t const num_of_eligble_edges = eligble_edges.size();
    for(std::size_t i=0 ; i<num_of_eligble_edges; ++i){
        EligbleEdge const& edge = eligble_edges[i];

        Crease const& edge_crease = plc->creases[edge.crease_index];
        
        if(edge_crease.front()->vertex1->isSharp){        
            auto const& [center, r] = corners_ball_tree.getBallNearestNeighbor(edge_crease.front()->vertex1->vertex);

            double const r_deeply = r*(1. - alpha);

            if(distance(edge[0], center) <= r_deeply && distance(edge[1], center) <= r_deeply){
                isDeleted[i] = true;
                continue;
            }
        }
        
        if(edge_crease.back()->vertex2->isSharp){        
            auto const& [center, r] = corners_ball_tree.getBallNearestNeighbor(edge_crease.back()->vertex2->vertex);

            double const r_deeply = r*(1. - alpha);

            if(distance(edge[0], center) <= r_deeply && distance(edge[1], center) <= r_deeply){
                isDeleted[i] = true;
            }
        }
    }
}

double EdgesRMPS::calculateSmoothnessLimitation(Vector3D const& p, EligbleEdge const& edge_sampled, Trees const& trees) const {
    VoroCrust_KD_Tree_Boundary const& corners_boundary_tree = trees.VC_kd_sharp_corners;
    VoroCrust_KD_Tree_Boundary const& edges_boundary_tree = trees.VC_kd_sharp_edges;
    VoroCrust_KD_Tree_Boundary const& faces_boundary_tree = trees.VC_kd_faces;

    // find nearest sharp corner
    double const dist_nearest_corner = corners_boundary_tree.distanceToNearestNeighbor(p, maxRadius);


    double initial_min_dist = std::min(dist_nearest_corner, maxRadius);
    std::vector<Vector3D> const parallel({edge_sampled[1] - edge_sampled[0]});

    double dist_non_cosmooth_edge = edges_boundary_tree.distanceToNearestNonCosmoothPoint(p, parallel, edge_sampled.crease_index, sharpTheta, initial_min_dist);

    initial_min_dist = std::min(initial_min_dist, dist_non_cosmooth_edge);

    double const dist_different_crease = edges_boundary_tree.distanceToNearestNeighborExcludingFeatures(p, {edge_sampled.crease_index}, initial_min_dist);

    dist_non_cosmooth_edge = std::min(dist_different_crease, dist_non_cosmooth_edge);

    initial_min_dist = std::min(dist_non_cosmooth_edge, initial_min_dist);

    // find nearest non cosmooth point on face
    std::vector<std::size_t> patches_to_exclude;
    
    double dist_non_cosmooth_face = std::numeric_limits<double>::max();
    for(auto const& patch_faces : plc->edges[edge_sampled.plc_index]->divided_faces){
        std::size_t const patch_index = patch_faces[0]->patch_index;
        patches_to_exclude.push_back(patch_index);
        
        std::vector<Vector3D> normals;
        normals.reserve(patch_faces.size());
        
        for(Face const& face : patch_faces){
            normals.push_back(face->calcNormal());
        }

        double const dist_neareset_non_cosmooth = faces_boundary_tree.distanceToNearestNonCosmoothPoint(p, normals, patch_index, sharpTheta, initial_min_dist);

        dist_non_cosmooth_face = std::min(dist_neareset_non_cosmooth, dist_non_cosmooth_face);
    }

    initial_min_dist = std::min(dist_non_cosmooth_edge, dist_non_cosmooth_face);

    double const dist_different_patch = faces_boundary_tree.distanceToNearestNeighborExcludingFeatures(p, patches_to_exclude, initial_min_dist);

    dist_non_cosmooth_face = std::min(dist_different_patch, dist_non_cosmooth_face);

    return std::min({dist_nearest_corner, dist_non_cosmooth_edge, dist_non_cosmooth_face});
}

bool EdgesRMPS::isEligbleEdgeDeeplyCoveredInEdgeBall(EligbleEdge const& edge, Trees const& trees, std::size_t const ball_index) const {
    auto const& [center, radius] = trees.ball_kd_edges.getBall(ball_index);

    double const r_deeply = radius * (1.0 - alpha);

    // if edge is deeply covered return true
    return (distance(center, edge[0]) <= r_deeply) && (distance(center, edge[1]) <= r_deeply);
}

bool EdgesRMPS::discardEligbleEdges(Trees const& trees){
    bool shrunkOtherStrataBalls = false;

    discardEligbleEdgesContainedInCornerBalls(trees);

    // discard in reverse order bacause each erase changes the indices
    for(long i = eligble_edges.size()-1; i >= 0; --i){
        if(isDeleted[i]) continue;
        EligbleEdge const& edge = eligble_edges[i];
        
        
        //! MIGHTBEUNECESSARY: this might be unecessary because I can take the nearest neighbor to one of the vertices
        // of the edge and use the r_max test (see FacesRMPS)
        auto const& [nn_edge_ball_center, nn_edge_ball_radius] = trees.ball_kd_edges.getBallNearestNeighbor(edge.edge[0]);

        double const r_max = (nn_edge_ball_radius + L_Lipschitz*distance(edge.edge[0], nn_edge_ball_center)) / (1.0 - L_Lipschitz); 

        auto const& balls_to_check_edges = trees.ball_kd_edges.radiusSearch(edge.edge[0], r_max);
        
        bool discard = false;
        for(auto const ball_index : balls_to_check_edges){
            //! MIGHTBEUNECESSARY: all relevent balls should already be in the same crease
            if(edge.crease_index == trees.ball_kd_edges.feature_index[ball_index]){
                discard = discard || isEligbleEdgeDeeplyCoveredInEdgeBall(edge, trees, ball_index);
            }                             
        }

        if(discard){
            isDeleted[i] = true;
        }
    }

    std::size_t not_deleted = 0;
    for(std::size_t i=0; i < isDeleted.size(); ++i){
        if(not isDeleted[i]) not_deleted++;
    }

    std::vector<EligbleEdge> new_eligble_edges(not_deleted, eligble_edges[0]);

    for(std::size_t i=0, j=0; i < isDeleted.size(); ++i){
        if(not isDeleted[i]){
            new_eligble_edges[j] = eligble_edges[i];
            ++j;
        }
    }
    
    eligble_edges = new_eligble_edges;

    return shrunkOtherStrataBalls;
}

double EdgesRMPS::calculateInitialRadius(Vector3D const& point, std::size_t const edge_index, Trees const& trees) const {

    VoroCrust_KD_Tree_Ball const& edges_ball_tree = trees.ball_kd_edges;
    
    EligbleEdge const& edge = eligble_edges[edge_index];
    
    // limitation from cosmoothness
    double const r_smooth = calculateSmoothnessLimitation(point, edge, trees);

    if(edges_ball_tree.empty()){
        return std::min(maxRadius, 0.49*r_smooth);
    }

    // limitation from Lipschitzness
    auto const& [q, r_q] = edges_ball_tree.getBallNearestNeighbor(point);
    double const dist = distance(point, q);

    return std::min({maxRadius, 0.49*r_smooth, r_q + L_Lipschitz*dist});
}

bool EdgesRMPS::doSampling(VoroCrust_KD_Tree_Ball &edges_ball_tree, Trees const& trees){

    bool resample = false;

    double total_len;
    std::vector<double> start_len;
    auto const& res = calculateTotalLengthAndStartLengthOfEligbleEdges();
    total_len = res.first;
    start_len = res.second; 
    
    std::cout << "total_len = " << total_len << ", num_of_eligble_edges = " << eligble_edges.size() << std::endl;

    std::size_t miss_counter = 0;

    while(not eligble_edges.empty()){
        auto const& [success, edge_index, p] = sampleEligbleEdges(total_len, start_len);

        if(!success) continue;

        if(miss_counter >= 100){
            // remake tree so seach is faster
            edges_ball_tree.remakeTree();

            divideEligbleEdges();
            //! IMPORTANT: resample needs to be on the right of the ||!!!!
            resample = discardEligbleEdges(trees) || resample;
    
            auto const& res = calculateTotalLengthAndStartLengthOfEligbleEdges();
            total_len = res.first;
            start_len = res.second;             
            miss_counter = 0;

            std::cout << "total_len = " << total_len << ", num_of_eligble_edges = " << eligble_edges.size() << std::endl;
            continue;
        }

        if(checkIfPointIsDeeplyCovered(p, trees)){
            miss_counter += 1;
            continue;
        }


        EligbleEdge const& edge = eligble_edges[edge_index];
        
        double radius = calculateInitialRadius(p, edge_index, trees);

        auto const& centers_in_new_ball_indices = edges_ball_tree.radiusSearch(p, radius);
        if(not centers_in_new_ball_indices.empty()){
            if(uni01_gen() < rejection_probability) continue;

            for(auto const center_index : centers_in_new_ball_indices){
                Vector3D const& center_in = edges_ball_tree.points[center_index];
                radius = std::min(radius, distance(p, center_in) / (1. - alpha));
            }
        }

        if(radius < 1e-8){
            std::cout << "trying to create a ball which is too small" << std::endl;
            std::cout << "at p = " << p.x << ", " << p.y << ", " << p.z << std::endl;
            std::cout << "radius = " << radius << std::endl;
            exit(1); 
        }

        std::cout << "edge sample " << edges_ball_tree.points.size() << ", r = " << radius << std::endl;

        edges_ball_tree.insert(p, edge[1]-edge[0], radius, edge.crease_index, edge.plc_index);
        
        miss_counter = 0;
    }

    return resample;
}
