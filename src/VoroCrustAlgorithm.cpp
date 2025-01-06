#include "VoroCrustAlgorithm.hpp"
#include <cmath>
#include <iostream>
#include <boost/random.hpp>

bool enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree, double const L_Lipschitz){
    std::size_t const num_of_points = ball_tree.size();
    
    std::size_t number_of_balls_shrunk = 0;

    std::cout << "Running enforceLipschitzness " << std::endl;
    // For each ball i enforce lipschitzness (going through all balls up to distance r_i/L_Lipschitz)
    for(std::size_t i = 0; i<num_of_points; ++i){
        auto const [p_i, r_i] = ball_tree.getBall(i);
        
        auto const& suspects = ball_tree.radiusSearch(p_i, r_i/L_Lipschitz);
        for(auto const j : suspects){

            auto const& [p_j, r_j] = ball_tree.getBall(j);
            double const dist = distance(p_i, p_j);
            
            if(ball_tree.ball_radii[i] > r_j + L_Lipschitz*dist){
                number_of_balls_shrunk++;
                ball_tree.ball_radii[i] = std::min(ball_tree.ball_radii[i], r_j + L_Lipschitz*dist);
            }
        }
    }

    std::cout << "enforceLipschitzness : " << number_of_balls_shrunk << " balls shrunk" << std::endl;
    return number_of_balls_shrunk > 0;
}

VoroCrust_KD_Tree_Ball makeSeedBallTree(std::vector<Seed> const& seeds){
    
    std::size_t const seeds_size = seeds.size();

    std::vector<Vector3D> points;
    std::vector<double> radii;

    points.reserve(seeds_size+1);
    radii.reserve(seeds_size+1);

    for(auto const& seed : seeds){
        points.push_back(seed.p);
        radii.push_back(seed.radius);
    }

    return VoroCrust_KD_Tree_Ball(points, 
                                  std::vector<Vector3D>(seeds_size, Vector3D(0.0, 0.0, 0.0)), 
                                  std::vector<std::size_t>(seeds_size, 0), 
                                  std::vector<std::size_t>(seeds_size, 0),
                                  radii);
}

std::pair<std::vector<Seed>, std::vector<Seed>> determineSeedsInOut(std::vector<Seed> const& seeds, PL_ComplexPtr const& plc) {
    if(seeds.empty()){
        throw std::runtime_error("determineSeedsInOut: seeds is empty");
    }

    std::vector<Seed> in_seeds{};
    std::vector<Seed> out_seeds{};
    
    in_seeds.reserve(seeds.size());
    out_seeds.reserve(seeds.size());

    for(auto const& seed : seeds){
        if(plc->determineLocation(seed.p) == PL_Complex::Location::IN) {
            in_seeds.push_back(seed);
        } else {
            out_seeds.push_back(seed);
        }
    }

    std::pair<std::vector<Seed>, std::vector<Seed>> inout_seeds = std::make_pair(in_seeds, out_seeds);

    return inout_seeds;
}

std::pair<std::vector<Seed>, std::vector<Seed>> randomSampleVolumeSeeds(PL_ComplexPtr const& plc, std::pair<std::vector<Seed>, std::vector<Seed>> const& inout_seeds, double const maxSize, Trees const& trees, double const L_Lipschitz) {

    Vector3D const empty_vec(0.0, 0.0, 0.0);
    
    if(inout_seeds.first.empty() or inout_seeds.second.empty()){
        throw std::runtime_error("randomSampleVolumeSeeds: inout_seeds is empty");
        }

    auto const [ll_x, ll_y, ll_z, ur_x, ur_y, ur_z] = plc->getBoundingBox();
        auto const len_x = ur_x - ll_x;
        auto const len_y = ur_y - ll_y;
        auto const len_z = ur_z - ll_z;

    auto const& in_seeds_tree = makeSeedBallTree(inout_seeds.first);
    auto const& out_seeds_tree = makeSeedBallTree(inout_seeds.second);

    VoroCrust_KD_Tree_Ball in_volume_seeds_tree;
    VoroCrust_KD_Tree_Ball out_volume_seeds_tree;

        // lightweight dart-throwing
        boost::random::variate_generator uni01_gen(boost::mt19937_64(std::time(nullptr)), boost::random::uniform_01<>());

        std::size_t num_of_samples = 0;
    PL_Complex::Location location;
    
    bool lipschitzness_in_volume  = false;
    bool lipschitzness_out_volume = false;
    
        do {
            std::size_t miss_counter = 0;
            while(miss_counter < 1000){
                Vector3D const p(ll_x + len_x*uni01_gen(),
                                 ll_y + len_y*uni01_gen(),
                                 ll_z + len_z*uni01_gen());
                
            location = plc->determineLocation(p);

            bool in_boundary_ball = in_seeds_tree.isContainedInBall(p) || trees.ball_kd_faces.isContainedInBall(p);

                if(not trees.ball_kd_vertices.empty()) in_boundary_ball = in_boundary_ball || trees.ball_kd_vertices.isContainedInBall(p);
                if(not trees.ball_kd_edges.empty()) in_boundary_ball = in_boundary_ball || trees.ball_kd_edges.isContainedInBall(p);

                if(in_boundary_ball){
                    ++miss_counter;
                    continue;
                }

            if(location == PL_Complex::Location::IN){
                auto const& [p_s, r_s] = in_seeds_tree.getBallNearestNeighbor(p);

                double r_volume = std::numeric_limits<double>::max();
                if(not in_volume_seeds_tree.empty()){
                    if(in_volume_seeds_tree.isContainedInBall(p)){
                        ++miss_counter;
                        continue;
                    }

                    auto const& [p_nearest, r_nearest] = in_volume_seeds_tree.getBallNearestNeighbor(p);

                    r_volume = r_nearest + L_Lipschitz*distance(p, p_nearest);
                }

                double const r = std::min({maxSize, r_s + L_Lipschitz*distance(p, p_s), r_volume});

                in_volume_seeds_tree.insert(p, empty_vec, r, 0, 0);
            } else {
                auto const& [p_s, r_s] = out_seeds_tree.getBallNearestNeighbor(p);

                double r_volume = std::numeric_limits<double>::max();
                if(not out_volume_seeds_tree.empty()){
                    if(out_volume_seeds_tree.isContainedInBall(p)){
                        ++miss_counter;
                        continue;
                    }

                    auto const& [p_nearest, r_nearest] = out_volume_seeds_tree.getBallNearestNeighbor(p);

                    r_volume = r_nearest + L_Lipschitz*distance(p, p_nearest);
                }

                double const r = std::min({maxSize, r_s + L_Lipschitz*distance(p, p_s), r_volume});

                out_volume_seeds_tree.insert(p, empty_vec, r, 0, 0);
            }

                num_of_samples++;
                
                if(num_of_samples % 10000 == 0){
                    std::cout << "sample: " << num_of_samples << ", miss_counter: " << miss_counter << std::endl;
                in_volume_seeds_tree.remakeTree();
                out_volume_seeds_tree.remakeTree();
                }

                miss_counter = 0;
            }

        in_volume_seeds_tree.remakeTree();
        out_volume_seeds_tree.remakeTree();

        lipschitzness_in_volume  = enforceLipschitzness(in_volume_seeds_tree, L_Lipschitz);
        lipschitzness_out_volume = enforceLipschitzness(out_volume_seeds_tree, L_Lipschitz);
    } while(lipschitzness_in_volume || lipschitzness_out_volume);

    std::vector<Seed> in_volume_temp_seeds = getSeedsFromBallTree(in_volume_seeds_tree);
    std::vector<Seed> out_volume_temp_seeds = getSeedsFromBallTree(out_volume_seeds_tree);
        
    std::pair<std::vector<Seed>, std::vector<Seed>> inout_volume_seeds = std::pair<std::vector<Seed>, std::vector<Seed>>(in_volume_temp_seeds, out_volume_temp_seeds);

    return inout_volume_seeds;
}

std::vector<Seed> getSeedsFromBallTree(VoroCrust_KD_Tree_Ball const& ball_tree){
    auto const tree_size = ball_tree.size();

    std::vector<Seed> seeds;
    seeds.reserve(tree_size + 1);

    for(std::size_t i=0; i<tree_size; ++i){
        seeds.emplace_back(ball_tree.getBall(i));
    }

    return seeds;
}

void dumpSeeds(std::string const& dirname, std::vector<Seed> const& seeds){
    std::filesystem::path dirpath(dirname);
    std::filesystem::create_directory(dirpath);

    auto const num_of_seeds = seeds.size();

    std::vector<double> x, y, z, r;

    x.reserve(num_of_seeds+1);
    y.reserve(num_of_seeds+1);
    z.reserve(num_of_seeds+1);
    r.reserve(num_of_seeds+1);

    for(auto const& seed : seeds){
        auto const& p = seed.p;
        x.push_back(p.x);
        y.push_back(p.y);
        z.push_back(p.z);
        r.push_back(seed.radius);
    }

    dump_vector(dirpath / "x.txt", x);
    dump_vector(dirpath / "y.txt", y);
    dump_vector(dirpath / "z.txt", z);
    dump_vector(dirpath / "r.txt", r);
}

std::vector<Seed> load_dumpSeeds(std::filesystem::path const& dirname){
    if(not std::filesystem::is_directory(dirname)){
        throw std::runtime_error("Seeds dump directory does not exist");
    }

    auto const x = load_dump_vector<double>(dirname / "x.txt");
    auto const y = load_dump_vector<double>(dirname / "y.txt");
    auto const z = load_dump_vector<double>(dirname / "z.txt");
    auto const r = load_dump_vector<double>(dirname / "r.txt");

    auto const num_of_seeds = x.size();

    if(y.size() != num_of_seeds || z.size() != num_of_seeds || r.size() != num_of_seeds){
        std::cout << "ERROR: x,y,z,r sizes do not match!" << std::endl;
        throw std::runtime_error("ERROR: x,y,z,r sizes do not match!");
    }

    std::vector<Seed> seeds;
    seeds.reserve(num_of_seeds+1);

    for(std::size_t i=0; i<num_of_seeds; ++i){
        seeds.emplace_back(Vector3D(x[i], y[i], z[i]), r[i]);
    }

    return seeds;
}

