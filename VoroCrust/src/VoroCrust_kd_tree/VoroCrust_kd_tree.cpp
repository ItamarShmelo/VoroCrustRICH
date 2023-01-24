#include "VoroCrust_kd_tree.hpp"
#include <functional>
#include <numeric>



std::shared_ptr<Node> newNode(int index, int axis){
    std::shared_ptr<Node> node = std::make_shared<Node>();

    node->index = index;
    node->axis = axis;

    return node;
}

VoroCrust_KD_Tree::VoroCrust_KD_Tree() : root(nullptr), points() {}

VoroCrust_KD_Tree::VoroCrust_KD_Tree(std::vector<Vector3D> const& points) : points(points), root(nullptr) {
    makeTree(points);
}

void VoroCrust_KD_Tree::clear(){
    root.reset();
    points.clear();
}

void VoroCrust_KD_Tree::makeTree(std::vector<Vector3D> const& points_){
    clear();

    points = points_;
    std::vector<int> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    root = buildRecursive(indices.data(), static_cast<int>(points.size()), 0);
}
