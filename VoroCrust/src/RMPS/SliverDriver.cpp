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
