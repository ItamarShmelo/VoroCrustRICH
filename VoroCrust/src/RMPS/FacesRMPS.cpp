#include "FacesRMPS.hpp"

FacesRMPS::FacesRMPS(double const maxRadius_, 
                     double const L_Lipschitz_, 
                     double const alpha_, 
                     double const sharpTheta_, 
                     std::shared_ptr<PL_Complex> const& plc_) : 
                                    maxRadius(maxRadius_), 
                                    L_Lipschitz(L_Lipschitz_), 
                                    alpha(alpha_), 
                                    sharpTheta(sharpTheta_), 
                                    uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()), 
                                    plc(plc_), 
                                    eligble_faces(),
                                    isDeleted() {}

void FacesRMPS::loadFaces(std::vector<Face> const& faces) {
    if(not eligble_faces.empty()){
        std::cout << "eligble faces is not empty when loading faces" << std::endl;
        exit(1);
    }

    for(Face const& face : faces){
        std::vector<Vector3D> vertices;
        for(Vertex const& vertex : face->vertices) {
            vertices.push_back(vertex->vertex);
        }

        eligble_faces.push_back(EligbleFace(vertices, face->patch_index, face->index, face->calcArea()));
    }

    isDeleted = std::vector<bool>(eligble_faces.size(), false);
}

std::pair<double const, std::vector<double> const> FacesRMPS::calculateTotalAreaAndStartAreaOfEligbleFaces() const {
    std::size_t const num_of_eligble_faces = eligble_faces.size();

    std::vector<double> start_area(num_of_eligble_faces, 0.0);
    double total_area = 0.0;

    for(std::size_t i=0; i<num_of_eligble_faces; ++i){
        EligbleFace const& face = eligble_faces[i];
        start_area[i] = total_area;
        total_area += face.area;
    }

    return std::pair<double const, std::vector<double> const>(total_area, start_area);
}

void FacesRMPS::divideEligbleFaces() {
    //! ASSUMPTION: This fundtion assumes that every face is a triangle!
    std::size_t const eligble_faces_size = eligble_faces.size();

    std::vector<EligbleFace> new_eligble_faces(4*eligble_faces_size, EligbleFace(std::vector<Vector3D>(3, Vector3D(0, 0, 0)), 0.0, 0.0, 0.0));

    for(std::size_t i=0; i < eligble_faces_size; ++i){
        EligbleFace const& face = eligble_faces[i];
        Vector3D const& midpoint_01 = 0.5*(face[0] + face[1]);
        Vector3D const& midpoint_12 = 0.5*(face[1] + face[2]);
        Vector3D const& midpoint_20 = 0.5*(face[2] + face[0]);

        double const sub_triangle_area = 0.25*face.area;

        new_eligble_faces[4*i][0] = face[0];
        new_eligble_faces[4*i][1] = midpoint_01;
        new_eligble_faces[4*i][2] = midpoint_20;

        new_eligble_faces[4*i+1][0] = face[1];
        new_eligble_faces[4*i+1][1] = midpoint_12;
        new_eligble_faces[4*i+1][2] = midpoint_01;

        new_eligble_faces[4*i+2][0] = face[2];
        new_eligble_faces[4*i+2][1] = midpoint_20;
        new_eligble_faces[4*i+2][2] = midpoint_12;

        new_eligble_faces[4*i+3][0] = midpoint_01;
        new_eligble_faces[4*i+3][1] = midpoint_12;
        new_eligble_faces[4*i+3][2] = midpoint_20;
        
        //! MAYBEINSERTTOLOOP:
        new_eligble_faces[4*i].patch_index = face.patch_index;
        new_eligble_faces[4*i+1].patch_index = face.patch_index;
        new_eligble_faces[4*i+2].patch_index = face.patch_index;
        new_eligble_faces[4*i+3].patch_index = face.patch_index;

        new_eligble_faces[4*i]  .plc_index = face.plc_index;
        new_eligble_faces[4*i+1].plc_index = face.plc_index;
        new_eligble_faces[4*i+2].plc_index = face.plc_index;
        new_eligble_faces[4*i+3].plc_index = face.plc_index;

        new_eligble_faces[4*i].area = sub_triangle_area;
        new_eligble_faces[4*i+1].area = sub_triangle_area;
        new_eligble_faces[4*i+2].area = sub_triangle_area;
        new_eligble_faces[4*i+3].area = sub_triangle_area;
    }

    eligble_faces = new_eligble_faces;
    isDeleted = std::vector<bool>(eligble_faces.size(), false);
}

bool FacesRMPS::checkIfPointIsDeeplyCovered(Vector3D const& p, Trees const& trees) const {
    VoroCrust_KD_Tree_Ball const& faces_ball_tree = trees.ball_kd_faces;

    if(not faces_ball_tree.empty()){
        std::size_t const nn_index = faces_ball_tree.nearestNeighbor(p);
        
        Vector3D const& q = faces_ball_tree.points[nn_index];
        double const r_q = faces_ball_tree.ball_radii[nn_index];

        // maximal radius for centers which balls can deeply cover p
        double const r_max = (r_q + L_Lipschitz*distance(p, q)) / (1.0 - L_Lipschitz); 

        auto const& suspects = faces_ball_tree.radiusSearch(p, r_max);

        for(auto const i : suspects) {
            Vector3D const& center = faces_ball_tree.points[i];
            double const r = faces_ball_tree.ball_radii[i];

            double const dist = distance(p, center);

            if(dist <= r*(1.0 - alpha)){
                return true;
            }
        }
    }
    
    VoroCrust_KD_Tree_Ball const& edges_ball_tree = trees.ball_kd_edges;
    
    std::size_t const nn_index = edges_ball_tree.nearestNeighbor(p);
    
    Vector3D const& q = edges_ball_tree.points[nn_index];
    double const r_q = edges_ball_tree.ball_radii[nn_index];

    // maximal radius for cetners which balls can deeply cover p
    double const r_max = (r_q + L_Lipschitz*distance(p, q)) / (1.0 - L_Lipschitz); 
    
    auto const& suspects = edges_ball_tree.radiusSearch(p, r_max);

    for(auto const i : suspects) {
        Vector3D const& center = edges_ball_tree.points[i];
        double const r = edges_ball_tree.ball_radii[i];

        double const dist = distance(p, center);
        if(dist <= r*(1.0 - alpha)){
            return true;
        }
    }

    VoroCrust_KD_Tree_Ball const& corners_ball_tree = trees.ball_kd_vertices;

    std::size_t const nn_corner_index = corners_ball_tree.nearestNeighbor(p);

    Vector3D const& center = corners_ball_tree.points[nn_corner_index];
    double const radius = corners_ball_tree.ball_radii[nn_corner_index];

    return distance(p, center) <= radius*(1.0-alpha);
}

std::tuple<bool, std::size_t const, Vector3D const> FacesRMPS::sampleEligbleFaces(double const total_area, std::vector<double> const& start_area) {
    double const sample = uni01_gen()*total_area;

    // find on which face the sample falls;
    auto const& iter_lower_bound = std::lower_bound(start_area.begin(), start_area.end(), sample);
    std::size_t const face_index = std::distance(start_area.begin(), iter_lower_bound) - 1;

    // sample point inside triangle
    //! CODEDUPLICATION: Done in the same way in Trees::superSampleFaces
    //! ASSUMESFACEISATRIANGLE: assumes face is a triangle

    // sample point in a triangle;
    double const sqrt_r1 = std::sqrt(uni01_gen());
    double const r2 = uni01_gen();

    // sample to close to boundary return success=false
    if(sqrt_r1 < 1e-14 || r2 < 1e-14) {
        return std::tuple<bool, std::size_t const, Vector3D const>(false, 0, Vector3D(0.0, 0.0, 0.0));
    }

    EligbleFace const& face = eligble_faces[face_index];
    Vector3D const& A = face[0];
    Vector3D const& B = face[1];
    Vector3D const& C = face[2];

    // sample point formula is taken from https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle-in-3d
    Vector3D const& point = (1.0-sqrt_r1)*A + (sqrt_r1*(1.0-r2))*B + (r2*sqrt_r1)*C;
    
    return std::tuple<bool, std::size_t const, Vector3D const>(true, face_index, point);
}

void FacesRMPS::discardEligbleFacesContainedInCornerBalls(Trees const& trees) {
    VoroCrust_KD_Tree_Ball const& corners_ball_tree = trees.ball_kd_vertices;


    std::size_t const eligble_faces_size = eligble_faces.size();

    for(std::size_t i=0; i<eligble_faces_size; ++i){
        EligbleFace const& face = eligble_faces[i];
        SurfacePatch const& patch = plc->patches[face.patch_index];

        for(std::size_t const corner_index : patch.patch_corners) {
            Vertex const& corner = plc->vertices[corner_index];

            std::size_t const nn_index = corners_ball_tree.nearestNeighbor(corner->vertex);

            Vector3D const& center = corners_ball_tree.points[nn_index];
            double const r = corners_ball_tree.ball_radii[nn_index]*(1. - alpha);

            if(face.isContainedInBall(center, r)){
                isDeleted[i] = true;
                break;
            }
        }
    }
}

bool EligbleFace::isContainedInBall(Vector3D const& center, double const r) const {

    for(Vector3D const& vertex : face)
        if(distance(center, vertex) > r)
            return false;
    
    return true;
}

void FacesRMPS::discardEligbleFacesContainedInEdgeBalls(Trees const& trees) {
    VoroCrust_KD_Tree_Ball const& edges_ball_tree = trees.ball_kd_edges;


    std::size_t const eligble_faces_size = eligble_faces.size();

    for(std::size_t i=0; i < eligble_faces_size; ++i){
        if(isDeleted[i]) continue;
        EligbleFace const& face = eligble_faces[i];
        // find the closes ball to the first vertex
        Vector3D const& v1 = face.face[0];        
        std::size_t const nn_index = edges_ball_tree.nearestNeighbor(v1);
        
        Vector3D const& q = edges_ball_tree.points[nn_index];
        double const r_q = edges_ball_tree.ball_radii[nn_index];
        
        // maximal radius for centers which balls can cover v1 hence the face
        double const r_max = (r_q + L_Lipschitz*distance(v1, q)) / (1.0 - L_Lipschitz);
        auto const& suspects = edges_ball_tree.radiusSearch(v1, r_max);

        for(auto const j : suspects) {
            Vector3D const& center = edges_ball_tree.points[j];
            double const r = edges_ball_tree.ball_radii[j]*(1. - alpha);

            if(face.isContainedInBall(center, r)){
                isDeleted[i] = true;
                break;
            }
        }
    }
}

double FacesRMPS::calculateSmoothnessLimitation(Vector3D const& p, EligbleFace const& face_sampled, Trees const& trees) const {
    // Corners
    VoroCrust_KD_Tree_Boundary const& corners_boundary_tree = trees.VC_kd_sharp_corners;

    auto const nearestSharpCorner_index = corners_boundary_tree.nearestNeighbor(p);
    Vector3D const& nearestSharpCorner = corners_boundary_tree.points[nearestSharpCorner_index];

    double const dist_nearest_corner = distance(p, nearestSharpCorner);
    
    // Edges
    VoroCrust_KD_Tree_Boundary const& edges_boundary_tree = trees.VC_kd_sharp_edges;

    auto const nearestPointOnSharpEdge_index = edges_boundary_tree.nearestNeighbor(p);
    Vector3D const& nearestPointOnSharpEdge = edges_boundary_tree.points[nearestPointOnSharpEdge_index];

    double const dist_nearest_sharp_edge = distance(p, nearestPointOnSharpEdge);

    // Faces
    double dist_non_cosmooth_face = std::numeric_limits<double>::max();

    VoroCrust_KD_Tree_Boundary const& faces_boundary_tree = trees.VC_kd_faces;
    Face const& plc_face = plc->faces[face_sampled.plc_index];

    long const nn_non_cosmooth_index = faces_boundary_tree.nearestNonCosmoothPoint(p, plc_face->calcNormal(), face_sampled.patch_index, sharpTheta);

    if(nn_non_cosmooth_index >= 0){
        Vector3D const& nn_non_cosmooth_p = faces_boundary_tree.points[nn_non_cosmooth_index];

        dist_non_cosmooth_face = std::min(dist_non_cosmooth_face, distance(p, nn_non_cosmooth_p));
    }

    long const nn_different_patch = faces_boundary_tree.nearestNeighborExcludingFeatures(p, {face_sampled.patch_index});

    if(nn_different_patch >= 0){
        Vector3D const& nn_different_patch_p = faces_boundary_tree.points[nn_different_patch];

        dist_non_cosmooth_face = std::min(dist_non_cosmooth_face, distance(p, nn_different_patch_p));
    }

    return std::min({dist_nearest_corner, dist_nearest_sharp_edge, dist_non_cosmooth_face});
}

bool FacesRMPS::isEligbleFaceDeeplyCoveredInFaceBall(EligbleFace const& face, Trees const& trees, std::size_t const ball_index) const {
    VoroCrust_KD_Tree_Ball const& faces_ball_tree = trees.ball_kd_faces;

    Vector3D const& center = faces_ball_tree.points[ball_index];
    double const r_deeply = faces_ball_tree.ball_radii[ball_index] * (1.0 - alpha);

    return face.isContainedInBall(center, r_deeply);
}

bool FacesRMPS::discardEligbleFaces(Trees const& trees) {
    bool shrunkOtherStataBalls = false;

    discardEligbleFacesContainedInCornerBalls(trees);
    discardEligbleFacesContainedInEdgeBalls(trees);

    VoroCrust_KD_Tree_Ball const& faces_ball_tree = trees.ball_kd_faces;
    // erase in reverse so indices wont change
    // might be better to put the indices in a vector then discard all at once
    for(long i=eligble_faces.size()-1; i>=0 ; --i){
        if(isDeleted[i]) continue;
        EligbleFace const& face = eligble_faces[i];

        std::size_t const i_nn_face_ball = faces_ball_tree.nearestNeighbor(face.face[0]);
        Vector3D const& nn_face_ball_center = faces_ball_tree.points[i_nn_face_ball];
        double const nn_face_ball_radius = faces_ball_tree.ball_radii[i_nn_face_ball];

        double const r_max = (nn_face_ball_radius + L_Lipschitz*distance(face.face[0], nn_face_ball_center)) / (1.0 - L_Lipschitz); 

        auto const& balls_to_check_faces = faces_ball_tree.radiusSearch(face.face[0], r_max);

        bool discard = false;
        
        for(auto const ball_index : balls_to_check_faces) {
            //! MIGHTBEUNECESSARY: all relevent balls should already be in the same patch
            discard = discard || isEligbleFaceDeeplyCoveredInFaceBall(face, trees, ball_index);
        }

        if(discard){
            isDeleted[i] = true;
        }
    }

    std::size_t not_deleted = 0;
    for(std::size_t i=0; i < isDeleted.size(); ++i){
        if(not isDeleted[i]) not_deleted++;
    }

    std::vector<EligbleFace> new_eligble_faces(not_deleted, EligbleFace(std::vector<Vector3D>(3, Vector3D(0, 0, 0)), 0.0, 0.0, 0.0));

    for(std::size_t i=0, j=0; i < isDeleted.size(); ++i){
        if(not isDeleted[i]){
            new_eligble_faces[j] = eligble_faces[i];
            ++j;
        }
    }

    eligble_faces = new_eligble_faces;

    return shrunkOtherStataBalls;
}

double FacesRMPS::calculateInitialRadius(Vector3D const& p, std::size_t const face_index, Trees const& trees) const {

    VoroCrust_KD_Tree_Ball const& faces_ball_tree = trees.ball_kd_faces;

    EligbleFace const& face = eligble_faces[face_index];

    //limitation from cosmoothness
    double const r_smooth = calculateSmoothnessLimitation(p, face, trees);

    if(faces_ball_tree.empty()){
        return std::min(maxRadius, 0.49*r_smooth);
    }

    // limitation from Lipschitzness
    std::size_t const nn_face_ball = faces_ball_tree.nearestNeighbor(p);
    Vector3D const& q = faces_ball_tree.points[nn_face_ball];
    double const r_q = faces_ball_tree.ball_radii[nn_face_ball];

    double const dist  = distance(p, q);

    return std::min({maxRadius, 0.49*r_smooth, r_q + L_Lipschitz*dist});
}

bool FacesRMPS::doSampling(VoroCrust_KD_Tree_Ball &faces_ball_tree, Trees const& trees){
    bool resample = false;

    double total_area;
    std::vector<double> start_area;

    auto const& res = calculateTotalAreaAndStartAreaOfEligbleFaces();
    total_area = res.first;
    start_area = res.second;

    std::cout << "total_area = " << total_area << ", num_of_eligble_faces = " << eligble_faces.size() << std::endl;

    std::size_t miss_counter = 0;

    while(not eligble_faces.empty()) {
        
        auto const& [success, face_index, p] = sampleEligbleFaces(total_area, start_area);

        if(!success) continue;

        if(miss_counter >= 100){
            // remake tree so search is faster
            faces_ball_tree.remakeTree();
            
            divideEligbleFaces();
            //! IMPORTANT: resample needs to be on the right of the ||!!!!
            resample = discardEligbleFaces(trees) || resample; 

            auto const& res = calculateTotalAreaAndStartAreaOfEligbleFaces();
            total_area = res.first;
            start_area = res.second;
            
            miss_counter = 0;
            std::cout << "total_area = " << total_area << ", num_of_eligble_faces = " << eligble_faces.size() << std::endl;
            continue;
        }

        if(checkIfPointIsDeeplyCovered(p, trees)){
            miss_counter += 1;
            continue;
        }
        

        EligbleFace const& face = eligble_faces[face_index];

        double radius = calculateInitialRadius(p, face_index, trees);

        auto const& centers_in_new_ball_indices = faces_ball_tree.radiusSearch(p, radius);

        if(not centers_in_new_ball_indices.empty()){
            if(uni01_gen() < rejection_probability) continue;

            for(auto const center_index : centers_in_new_ball_indices){
                Vector3D const& center_in = faces_ball_tree.points[center_index];
                radius = std::min(radius, distance(p, center_in) / (1. - alpha));
            }
        }

        if(radius < 1e-8){
            std::cout << "trying to create a ball which is too small" << std::endl;
            std::cout << "at p = " << p.x << ", " << p.y << ", " << p.z << std::endl;
            std::cout << "radius = " << radius << std::endl;
            exit(1);
        }

        std::cout << "face sample " << faces_ball_tree.points.size() << ", r = " << radius << std::endl;

        Face const& plc_face = plc->faces[face.plc_index];
        faces_ball_tree.insert(p, plc_face->calcNormal(), radius, face.patch_index, face.plc_index);

        miss_counter = 0;
    }

    return resample;
}