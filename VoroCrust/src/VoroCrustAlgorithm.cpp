#include "VoroCrustAlgorithm.hpp"
#include <cmath>
#include <iostream>
#include <boost/random.hpp>

VoroCrustAlgorithm::VoroCrustAlgorithm( PL_Complex const& plc_, 
                                        double const sharpTheta_, 
                                        double const flatTheta_, 
                                        double const maxRadius_,
                                        double const L_Lipschitz_,
                                        double const alpha_): plc(std::make_shared<PL_Complex>(plc_)), 
                                                              trees(),
                                                              sharpTheta(sharpTheta_),
                                                              flatTheta(flatTheta_),
                                                              maxRadius(maxRadius_),
                                                              L_Lipschitz(L_Lipschitz_),
                                                              alpha(alpha_),
                                                              maximal_num_iter(15),
                                                              cornersDriver(maxRadius_, L_Lipschitz_, sharpTheta_, plc),
                                                              edgesDriver(maxRadius_, L_Lipschitz_, alpha_, sharpTheta_, plc),
                                                              facesDriver(maxRadius_, L_Lipschitz_, alpha_, sharpTheta_, plc),
                                                              sliverDriver(L_Lipschitz_) {

    if(sharpTheta > M_PI_2){
        std::cout << "ERROR: sharpTheta > pi/2" << std::endl;
        exit(1);
    }           

    if(L_Lipschitz >= 1){
        std::cout << "ERROR: L_Lipschitz >= 1" << std::endl;
        exit(1);
    }

    if(L_Lipschitz <= 0){
        std::cout << "ERROR: L_Lipschitz <= 0 " << std::endl;
        exit(1);
    }

}

void VoroCrustAlgorithm::run() {
    
    //! TODO: Maybe all this need to be in the plc under detect features?
    if(not plc->checkAllVerticesAreUnique()) exit(1);

    if(not plc->checkIfALLFacesAreFlat()) exit(1);

    if(not plc->checkAllVerticesAreOnFace()) exit(1);

    plc->detectFeatures(sharpTheta, flatTheta);

    //! TODO: make sampling size a user input!
    //! IMPORTANT: sampling size can effect the convergence of the algorithm because the radius is determined using proximity to the sampled points on different features. Make sure that the sampling size is compatible to the size of the smallest polygon in the data. One wants the sampling to be "dense" in the edges and faces.
    trees.loadPLC(*plc, 1e5, 1e6);    

    //! TODO: init eligable edges vertices and faces
    cornersDriver.loadCorners(plc->sharp_corners);
    cornersDriver.doSampling(trees.ball_kd_vertices, trees);
    trees.ball_kd_vertices.remakeTree();
    // sliver elimination loop
    for(std::size_t iteration = 0; iteration < maximal_num_iter; ++iteration){
        // enfore lipchitzness on vertices
        std::cout << "\nCorners Lipchitzness\n--------------\n" << std::endl;
        enforceLipschitzness(trees.ball_kd_vertices);    
        
        do {
            std::cout << "\nEdgesRMPS\n--------------\n" << std::endl;
            edgesDriver.loadEdges(plc->sharp_edges);
            edgesDriver.doSampling(trees.ball_kd_edges, trees);
            trees.ball_kd_edges.remakeTree();
        } while(enforceLipschitzness(trees.ball_kd_edges));
            
        do {
            std::cout << "\nFacesRMPS\n--------------\n" << std::endl;
            facesDriver.loadFaces(plc->faces);
            facesDriver.doSampling(trees.ball_kd_faces, trees);
            trees.ball_kd_faces.remakeTree();
        } while(enforceLipschitzness(trees.ball_kd_faces));
        
        if(not sliverDriver.eliminateSlivers(trees)) break;

        enforceLipschitzness(trees.ball_kd_edges);
        enforceLipschitzness(trees.ball_kd_faces);
    }
}

bool VoroCrustAlgorithm::enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree){
    long const num_of_points = ball_tree.points.size();
    
    bool isBallsShrunk = false;
    std::size_t number_of_balls_shrunk = 0;

    std::cout << "Running enforceLipschitzness " << std::endl;
    // go through all pairs i != j and enforce r_i <= r_j + L * ||p_i - p_j||
    //! MAYBE: only consider the overlapping balls or balls up to radius r_i?
    for(long i = 0; i<num_of_points; ++i){
        Vector3D const& p_i = ball_tree.points[i];
        double const r_i = ball_tree.ball_radii[i];
        std::vector<int> const& suspects = ball_tree.radiusSearch(p_i, 1.5*r_i);

        for(int const j : suspects){
            if(i == j) continue;
            double const dist = distance(p_i, ball_tree.points[j]);
            
            if(r_i > ball_tree.ball_radii[j] + L_Lipschitz*dist){
                isBallsShrunk = true;
                number_of_balls_shrunk++;
                ball_tree.ball_radii[i] = ball_tree.ball_radii[j] + L_Lipschitz*dist;
            }
        }
    }

    std::cout << "enforceLipschitzness : " << number_of_balls_shrunk << " balls shrunk" << std::endl;
    return isBallsShrunk;
}

std::vector<Seed> VoroCrustAlgorithm::getSeeds() const {
    return sliverDriver.getSeeds(trees);
}

std::pair<std::vector<Seed>, std::vector<Seed>> VoroCrustAlgorithm::determineIfSeedsAreInsideOrOutside(std::vector<Seed> const& seeds) const {

    // if(seeds.size() % 2 != 0){
    //     std::cout << "ERROR: seeds is not even" << std::endl;
    //     exit(1);
    // }

    if(seeds.empty()){
        std::cout << "ERROR: seeds is empty" << std::endl;
        exit(1);
    }

    std::vector<Seed> in_seeds, out_seeds;
    in_seeds.reserve(seeds.size());
    out_seeds.reserve(seeds.size());

    int i = 0;
    for(auto const& seed : seeds) {
        auto const location = plc->determineLocation(seed.p);

        if(location == PL_Complex::Location::OUT){
            out_seeds.push_back(seed);
        } else {
            in_seeds.push_back(seed);
        }
    }

    return std::make_pair(in_seeds, out_seeds);
}

VoroCrust_KD_Tree_Ball makeSeedBallTree(std::vector<Seed> const& seeds){
    
    std::size_t const seeds_size = seeds.size();

    std::vector<Vector3D> points;
    std::vector<double> radii;

    points.reserve(seeds_size);
    radii.reserve(seeds_size);

    for(auto const& seed : seeds){
        points.push_back(seed.p);
        radii.push_back(seed.radius);
    }

    return VoroCrust_KD_Tree_Ball(points, std::vector<Vector3D>(seeds_size, Vector3D(0.0, 0.0, 0.0)), radii, std::vector<std::size_t>(seeds_size, 0), std::vector<std::size_t>(seeds_size, 0));
}

std::pair<std::vector<Vector3D>, std::vector<Vector3D>> VoroCrustAlgorithm::calcVolumeSeedsUniform(std::vector<Seed> const& seeds, std::size_t const num_points_x, std::size_t const num_points_y, std::size_t const num_points_z) const {
    auto const [ll_x, ll_y, ll_z, ur_x, ur_y, ur_z] = plc->getBoundingBox(); // upper right

    auto const len_x = ur_x - ll_x;
    auto const len_y = ur_y - ll_y;
    auto const len_z = ur_z - ll_z;

    auto const step_x = len_x / num_points_x;
    auto const step_y = len_y / num_points_y;
    auto const step_z = len_z / num_points_z;

    auto const total_num_points = num_points_x * num_points_y * num_points_z;
    auto [in_seeds_boundary, out_seeds_boundary] = determineIfSeedsAreInsideOrOutside(seeds);

    auto const in_seeds_boundary_size = in_seeds_boundary.size();
    std::vector<Vector3D> in_seeds, out_seeds;
    
    in_seeds.reserve(in_seeds_boundary_size + total_num_points + 100);
    out_seeds.reserve(out_seeds_boundary.size()+100);

    auto const& in_seeds_tree = makeSeedBallTree(in_seeds_boundary);

    for(std::size_t i=0; i < num_points_x; ++i){
        for(std::size_t j=0; j < num_points_y; ++j){
            for(std::size_t k=0; k < num_points_z; ++k){
                Vector3D const seed(ll_x + step_x*i, ll_y + step_y*j, ll_z + step_z*k);

                auto const location = plc->determineLocation(seed);

                if(location == PL_Complex::Location::IN){
                    // if seed is inside check that it is not contained in a boundary seed
                    bool add_seed = true;
                    auto index = in_seeds_tree.nearestNeighbor(seed);

                    auto const& p_in_seed = in_seeds_tree.points[index];
                    auto const r_in_seed = in_seeds_tree.ball_radii[index];

                    if(distance(p_in_seed, seed) > r_in_seed){
                        if(add_seed) in_seeds.push_back(seed);
                    }
                } 
            }
        }
    }
    
    for(auto const out_seed : out_seeds_boundary) {
        out_seeds.push_back(out_seed.p);
    }

    for(auto const in_seed : in_seeds_boundary){
        in_seeds.push_back(in_seed.p);
    }

    return std::make_pair(in_seeds, out_seeds);
}

std::string VoroCrustAlgorithm::repr() const {
    std::ostringstream s;
    
    s << "VoroCrustAlgorithm : \n--------------------------------\n\n";
    s << "PLC : \n------------\n" << plc->repr() << std::endl;

    return s.str();
}

std::pair<std::vector<Vector3D>, std::vector<Vector3D>> VoroCrustAlgorithm::calcVolumeSeedsNonUniform(std::vector<Seed> const& seeds) const {
    Vector3D const empty_vec(0.0, 0.0, 0.0);

    auto const [ll_x, ll_y, ll_z, ur_x, ur_y, ur_z] = plc->getBoundingBox();
    auto const len_x = ur_x - ll_x;
    auto const len_y = ur_y - ll_y;
    auto const len_z = ur_z - ll_z;

    auto const& [in_seeds_boundary, out_seeds_boundary] = determineIfSeedsAreInsideOrOutside(seeds);
    auto const& in_seeds_tree = makeSeedBallTree(in_seeds_boundary);
    
    VoroCrust_KD_Tree_Ball volume_seeds_tree;

    // lightweight dart-throwing
    boost::random::variate_generator uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>());

    std::size_t miss_counter = 0;
    std::size_t num_of_samples = 0;
    while(miss_counter < 100){
        Vector3D const p(ll_x+len_x*uni01_gen(), ll_y+len_y*uni01_gen(), ll_z+len_z*uni01_gen());
        
        if(plc->determineLocation(p) == PL_Complex::Location::OUT){
            miss_counter++;
            continue;
        }

        //! CODEDUPLICATION:
    
        auto index_s = in_seeds_tree.nearestNeighbor(p);
        auto const& p_s = in_seeds_tree.points[index_s];
        auto const r_s = in_seeds_tree.ball_radii[index_s];

        if(distance(p_s, p) < r_s){
            miss_counter++;
            continue;
        }
        
        
        double r_volume = std::numeric_limits<double>::max();
        if(not volume_seeds_tree.points.empty()){
            auto index_nearest = volume_seeds_tree.nearestNeighbor(p);
            auto const& p_nearest = volume_seeds_tree.points[index_nearest];
            auto const r_nearest = volume_seeds_tree.ball_radii[index_nearest];

            if(distance(p_nearest, p) < r_nearest){
                miss_counter++;
                continue;
            }

            r_volume = r_nearest + L_Lipschitz*distance(p, p_nearest);
        }
        
        double const r = std::min({maxRadius, r_s + L_Lipschitz*distance(p, p_s), r_volume});
        volume_seeds_tree.insert(p, empty_vec, r, 0, 0);
        miss_counter=0;
        num_of_samples++;
        std::cout << "sample : "  << num_of_samples << "\n";
    }

    volume_seeds_tree.remakeTree();
    
    auto in_seeds = std::move(volume_seeds_tree.points);
    auto in_seeds_boundry_vec = std::move(in_seeds_tree.points);

    in_seeds.insert(in_seeds.end(), in_seeds_boundry_vec.begin(), in_seeds_boundry_vec.end());

    std::vector<Vector3D> out_seeds;
    out_seeds.reserve(out_seeds_boundary.size()+100);

    for(auto const out_seed : out_seeds_boundary) {
        out_seeds.push_back(out_seed.p);
    }

    return std::pair(in_seeds, out_seeds);
}