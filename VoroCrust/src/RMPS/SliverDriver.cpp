#include "SliverDriver.hpp"
#include <iostream>
#include "../../../source/misc/utils.hpp"

SliverDriver::SliverDriver(double const L_Lipschitz_) : L_Lipschitz(L_Lipschitz_), r_new_corner_balls(), r_new_edge_balls(), r_new_face_balls(), number_of_slivers_eliminated(0) {}

std::vector<BallInfo> SliverDriver::groupOverlappingBalls(BallInfo const& ball_info, Trees const& trees) const {
    //! MAYBE: I need to change r_max when looking for lower dimensional balls;
    auto const& [p, radius] = getBall(ball_info, trees);

    std::vector<BallInfo> overlapping_balls;

    double const r_max = (2.0 / (1.0 - L_Lipschitz)) * radius;
    
    VoroCrust_KD_Tree_Ball const& corners_ball_tree = trees.ball_kd_vertices;
    std::vector<std::size_t> const& overlapping_corner_balls_indices = corners_ball_tree.getOverlappingBalls(p, radius, r_max);

    VoroCrust_KD_Tree_Ball const& edges_ball_tree = trees.ball_kd_edges;
    std::vector<std::size_t> const& overlapping_edge_balls_indices = edges_ball_tree.getOverlappingBalls(p, radius, r_max);

    VoroCrust_KD_Tree_Ball const& faces_ball_tree = trees.ball_kd_faces;
    std::vector<std::size_t> const& overlapping_face_balls_indices = faces_ball_tree.getOverlappingBalls(p, radius, r_max);


    for(std::size_t const ball_index : overlapping_corner_balls_indices){
        overlapping_balls.push_back(BallInfo(ball_index, Dim::CORNER));
    }

    for(std::size_t const ball_index : overlapping_edge_balls_indices){
        overlapping_balls.push_back(BallInfo(ball_index, Dim::EDGE));
    }

    for(std::size_t const ball_index : overlapping_face_balls_indices){
        overlapping_balls.push_back(BallInfo(ball_index, Dim::FACE));
    }

    auto const& it = std::find(overlapping_balls.begin(), overlapping_balls.end(), ball_info);

    if(it == overlapping_balls.end()){
        std::cout << "current ball is not in overlapping balls!" << std::endl;
        exit(1);
    }

    // erase the current ball from the overlapping balls
    int const index = it - overlapping_balls.begin();
    overlapping_balls.erase(overlapping_balls.begin() + index);

    return overlapping_balls;
}

std::vector<Triplet> SliverDriver::formTripletsOfOverlappingBalls(std::vector<BallInfo> const& overlapping_balls, Trees const& trees) const {
    std::vector<std::pair<std::size_t, std::size_t>> triplets;

    std::size_t overlapping_balls_size = overlapping_balls.size();

    for(std::size_t i=0 ; i < overlapping_balls_size; ++i){
        auto const& [p_i, r_i] = getBall(overlapping_balls[i], trees);
        for(std::size_t j=i+1; j < overlapping_balls_size; ++j){
            auto const& [p_j, r_j] = getBall(overlapping_balls[j], trees);
            
            if(distance(p_i, p_j) < r_i + r_j) triplets.push_back(Triplet(i, j));
        }
    }

    return triplets;
}

Ball SliverDriver::getBall(BallInfo const& ball_info, Trees const& trees) const {
    if(ball_info.dim == Dim::FACE){
        std::size_t const index = ball_info.index;
        return Ball(trees.ball_kd_faces.points[index], trees.ball_kd_faces.ball_radii[index]);
    }

    
    if(ball_info.dim == Dim::EDGE){
        std::size_t const index = ball_info.index;
        return Ball(trees.ball_kd_edges.points[index], trees.ball_kd_edges.ball_radii[index]);
    }

    
    if(ball_info.dim == Dim::CORNER){
        std::size_t const index = ball_info.index;
        return Ball(trees.ball_kd_vertices.points[index], trees.ball_kd_vertices.ball_radii[index]);
    }

    std::cout << "ERROR: getBall" << std::endl;
    exit(1);
}

void SliverDriver::eliminateSliversForBallsInBallTree(Dim const dim, Trees const& trees) {
    VoroCrust_KD_Tree const * ball_tree_ptr;

    switch(dim){
        case Dim::CORNER:
            ball_tree_ptr = &trees.ball_kd_vertices;
            break;
        case Dim::EDGE:
            ball_tree_ptr = &trees.ball_kd_edges;
            break;
        case Dim::FACE:
            ball_tree_ptr = &trees.ball_kd_faces;
            break;
        default:
            std::cout << "error in dim in eliminateSliversForBallsInBallTree" << std::endl;
            break;
    }
    
    std::size_t const num_of_balls = ball_tree_ptr->points.size();
    for(std::size_t i=0; i < num_of_balls; ++i){
        dealWithBall(BallInfo(i, dim), trees);
    }
}

void SliverDriver::dealWithBall(BallInfo const& ball_info, Trees const& trees) {
    std::vector<BallInfo> const& overlapping_balls = groupOverlappingBalls(ball_info, trees);
    std::vector<Triplet>  const& triplets = formTripletsOfOverlappingBalls(overlapping_balls, trees);

    for(Triplet const& triplet : triplets){
        dealWithTriplets(ball_info, triplet, overlapping_balls, trees);
    }
}

void SliverDriver::dealWithTriplets(BallInfo const& ball_info_1, Triplet const& triplet, std::vector<BallInfo> const& overlapping_balls, Trees const& trees) {
    BallInfo const& ball_info_2 = overlapping_balls[triplet.first];
    BallInfo const& ball_info_3 = overlapping_balls[triplet.second];

    Ball const& ball_1 = getBall(ball_info_1, trees);
    Ball const& ball_2 = getBall(ball_info_2, trees);
    Ball const& ball_3 = getBall(ball_info_3, trees);

    auto const& [seed_p, seed_m] = calculateIntersectionSeeds(ball_1, ball_2, ball_3);

    for(BallInfo const& ball_info_4 : overlapping_balls){
        if(ball_info_4 == ball_info_2 || ball_info_4 == ball_info_3) continue;
        auto const& [p4, r4] = getBall(ball_info_4, trees);

        bool const is_seed_p_covered = distance(seed_p, p4) < r4; // check if seed_p is covered by ball_4
        bool const is_seed_m_covered = distance(seed_m, p4) < r4; // check if seed_m is covered by ball_4

        // check if there is a half covered seed pair
        double r_new = std::numeric_limits<double>::max();

        //! WARNING: EPSILONTICA
        if(is_seed_p_covered && !is_seed_m_covered){
            r_new = distance(p4, seed_p) * (1.0 - 1e-14);
        } else if(is_seed_m_covered && !is_seed_p_covered){
            r_new = distance(p4, seed_m) * (1.0 - 1e-14);
        }

        setRadiusOfBall(r_new, ball_info_4, trees);
    }
}
