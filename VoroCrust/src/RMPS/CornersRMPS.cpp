#include "CornersRMPS.hpp"
#include <boost/random.hpp>

CornersRMPS::CornersRMPS(double const maxRadius_, double const L_Lipschitz_, double const sharpTheta_, std::shared_ptr<PL_Complex> const& plc_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), sharpTheta(sharpTheta_), plc(plc_), eligble_corners() {}


void CornersRMPS::loadCorners(std::vector<Vertex> const& sharp_corners){
    eligble_corners.clear();
    for(Vertex const& corner_ptr : sharp_corners)
        eligble_corners.push_back(corner_ptr);
}

void CornersRMPS::doSampling(VoroCrust_KD_Tree_Ball &corner_ball_tree, Trees const& trees){

    Vector3D const empty_vector(0, 0, 0);

    while(not eligble_corners.empty()){

        // sample a vertex
        EligbleCorner const sample = newSample();

        double const radius = calculateInitialRadius(sample, trees);
        
        std::cout << "corner sample " << corner_ball_tree.points.size() << ", r = " << radius << std::endl;
        
        corner_ball_tree.insert(sample->vertex, empty_vector, radius, corner_ball_tree.points.size(), sample->index);
    }
}

EligbleCorner CornersRMPS::newSample(){
    // create a random number generator to sample from 0 - (eligble_corners.size()-1)
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_int_distribution<> int_distribution(0, eligble_corners.size() - 1);
    boost::variate_generator<boost::mt19937, boost::random::uniform_int_distribution<>> rand_gen(rng, int_distribution);

    // sample a random index;
    std::size_t const index = rand_gen();
    EligbleCorner sample = eligble_corners[index];

    eligble_corners.erase(eligble_corners.begin() + index);

    return sample;
}

double CornersRMPS::calculateInitialRadius(EligbleCorner const& corner, Trees const& trees) const {
    VoroCrust_KD_Tree_Ball const& corner_ball_tree = trees.ball_kd_vertices;

    double r_q = std::numeric_limits<double>::max();
    double dist_q = 0.0;
    
    if(not corner_ball_tree.points.empty()){
        // find nearest ball center
        int nearestBall_index = corner_ball_tree.nearestNeighbor(corner->vertex);

        Vector3D const& nearestBallCenter = corner_ball_tree.points[nearestBall_index];
        dist_q = distance(corner->vertex, nearestBallCenter); //||p-q||
        r_q = corner_ball_tree.ball_radii[nearestBall_index];
    }
    
    double const dist_q_prime = calculateSmoothnessLimitation(corner, trees);

    return std::min<double>({maxRadius, 0.49*dist_q_prime, r_q + L_Lipschitz*dist_q});
}

double CornersRMPS::calculateSmoothnessLimitation(EligbleCorner const& corner, Trees const& trees) const {
    VoroCrust_KD_Tree_Boundary const& corner_boundary_tree = trees.VC_kd_sharp_corners;
    VoroCrust_KD_Tree_Boundary const& edges_boundary_tree = trees.VC_kd_sharp_edges;
    VoroCrust_KD_Tree_Boundary const& faces_boundary_tree = trees.VC_kd_faces;

    // find nearest sharp corner
    int nearestSharpCorner_index = corner_boundary_tree.kNearestNeighbors(corner->vertex, 2)[1];
    Vector3D const& nearestSharpCorner = corner_boundary_tree.points[nearestSharpCorner_index];

    double const dist_nearest_corner = distance(corner->vertex, nearestSharpCorner); //||p-q^*||
    double dist_nearest_non_cosmooth_edge = std::numeric_limits<double>::max();
    double dist_nearest_non_cosmooth_face = std::numeric_limits<double>::max();

    // find neareset non cosmooth point on the edges
    std::vector<std::size_t> creases_exclude;
    for(Edge const& edge : corner->edges){
        if(not edge->isSharp) continue;
        
        std::size_t const crease_index = edge->crease_index;
        creases_exclude.push_back(crease_index);

        long nn_noncosmooth_on_edge_index = edges_boundary_tree.nearestNonCosmoothPointEdge(corner->vertex, edge->vertex2->vertex - edge->vertex1->vertex, crease_index, sharpTheta);

        // no noncosmooth point on the crease (very likely)
        if(nn_noncosmooth_on_edge_index < 0) continue;

        Vector3D const& nn_non_cosmooth_p = edges_boundary_tree.points[nn_noncosmooth_on_edge_index];

        dist_nearest_non_cosmooth_edge = std::min(dist_nearest_non_cosmooth_edge, distance(corner->vertex, nn_non_cosmooth_p));
    }

    long const nn_different_crease = edges_boundary_tree.nearestNeighborExcludingFeatures(corner->vertex, creases_exclude);
    // found any nearest neighbors
    if(nn_different_crease >= 0){
        Vector3D const& nn_different_crease_p = edges_boundary_tree.points[nn_different_crease];

        dist_nearest_non_cosmooth_edge = std::min(dist_nearest_non_cosmooth_edge, distance(corner->vertex, nn_different_crease_p));
    }
    

    //find neareset non cosmooth point on the face
    std::vector<std::size_t> patches_to_exclude;
    for(Face const& face : corner->faces){
        std::size_t const patch_index = face->patch_index;
        patches_to_exclude.push_back(patch_index); // might have multiple times the same patch


        long const nn_noncosmooth_on_face_index = faces_boundary_tree.nearestNonCosmoothPointFace(corner->vertex, face->calcNormal(), patch_index, sharpTheta, M_PI_2-sharpTheta);

        // no noncosmooth point on the patch (very likely)
        if(nn_noncosmooth_on_face_index < 0) continue;

        Vector3D const& nn_non_cosmooth_p = faces_boundary_tree.points[nn_noncosmooth_on_face_index];
        
        dist_nearest_non_cosmooth_face = std::min(dist_nearest_non_cosmooth_face, distance(corner->vertex, nn_non_cosmooth_p));
    }

    long const nn_different_patch = faces_boundary_tree.nearestNeighborExcludingFeatures(corner->vertex, patches_to_exclude);

    // found any nearest neighbors
    if(nn_different_patch >= 0){
        Vector3D const& nn_different_patch_p = faces_boundary_tree.points[nn_different_patch];

        dist_nearest_non_cosmooth_face = std::min(dist_nearest_non_cosmooth_face, distance(corner->vertex, nn_different_patch_p));
    }

    return std::min({dist_nearest_corner, dist_nearest_non_cosmooth_edge, dist_nearest_non_cosmooth_face});
}