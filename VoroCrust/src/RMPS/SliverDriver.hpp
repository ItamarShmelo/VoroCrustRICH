#ifndef SLIVER_DRIVER
#define SLIVER_DRIVER

#include "../trees.hpp"
#include <array>

enum Dim {
    CORNER = 0,
    EDGE = 1,
    FACE = 2
};

struct BallInfo {
    std::size_t index;
    Dim dim;

    BallInfo(std::size_t const index_, Dim const dim_) : index(index_), dim(dim_) {}

    friend bool operator==(BallInfo const& lhs, BallInfo const& rhs);
};

using Ball = std::pair<Vector3D, double const>;
using Triplet = std::pair<std::size_t, std::size_t>; // Triplet is p, Triplet.first, Triplet.second
using InfoQuartet = std::array<BallInfo const, 4>;
using BallQuartet = std::array<Ball const, 4>;

inline std::tuple<double const, double const, double const, double const> getLineCoeff(double const x1, double const y1, double const z1, double const r1, double const x2, double const y2, double const z2, double const r2);

inline std::pair<double const, double const> getZDependency(double const a1, double const b1, double const c1, double const k1, double const a3, double const b3, double const c3, double const k3);

class SliverDriver {
    public:
        double const L_Lipschitz;
        std::vector<double> r_new_corner_balls;
        std::vector<double> r_new_edge_balls;
        std::vector<double> r_new_face_balls;

        std::size_t number_of_slivers_eliminated;
        
        SliverDriver(double const L_Lipschitz_);

        bool eliminateSlivers(Trees &trees);

        std::vector<Vector3D> getSeeds(Trees const& trees) const;

    private:
        void eliminateSliversForBallsInBallTree(Dim const dim, Trees const& trees);

        void dealWithBall(BallInfo const& ball_info, Trees const& trees);

        void dealWithTriplets(BallInfo const& ball_info_1, Triplet const& triplet, std::vector<BallInfo> const& overlapping_balls, Trees const& trees);

        std::pair<Vector3D, Vector3D> calculateIntersectionSeeds(Ball const& ball_1, Ball const& ball_2, Ball const& ball_3) const;

        std::vector<BallInfo> groupOverlappingBalls(BallInfo const& ball, Trees const& trees) const;

        std::vector<Triplet> formTripletsOfOverlappingBalls(std::vector<BallInfo> const& overlapping_balls, Trees const& trees) const;

        void dealWithHalfCoveredSeeds(InfoQuartet const& balls_info, BallQuartet const& balls, Trees const& trees);

        Ball getBall(BallInfo const& ball_info, Trees const& trees) const;

        void setRadiusOfBall(double const r_new, BallInfo const& ball_info, Trees const& trees);

};

#endif // SLIVER_DRIVER