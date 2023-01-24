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

std::shared_ptr<Node> VoroCrust_KD_Tree::buildRecursive(int * indices, int npoints, int depth){
    if(npoints <= 0){
        return nullptr;
    }

    int const axis = depth % 3;
    int const mid  = (npoints - 1) / 2;

    std::nth_element(indices, indices + mid, indices + npoints, [&](int lhs, int rhs){
        return points[lhs][axis] < points[rhs][axis];
    });

    std::shared_ptr<Node> node = newNode(indices[mid], axis);

    node->left = buildRecursive(indices, mid, depth + 1);
    node->right = buildRecursive(indices + mid + 1, npoints - mid - 1, depth+1);

    return node;
}

int VoroCrust_KD_Tree::nearestNeighbor(Vector3D const& query) const {
    int guess = 0;
    double minDist = std::numeric_limits<double>::max();

    nearestNeighborRecursive(query, root, &guess, &minDist);

    return guess;
}

void VoroCrust_KD_Tree::nearestNeighborRecursive(Vector3D const& query, std::shared_ptr<Node> const& node, int *guess, double *minDist) const {
    if(node.get() == nullptr) return;

    Vector3D const& train = points[node->index];

    double const dist = distance(train, query);

    if(dist < *minDist){
        *minDist = dist;
        *guess = node->index;
    }

    int const axis = node->axis;
    std::shared_ptr<Node> node_first = query[axis] < train[axis] ? node->left : node->right;

    nearestNeighborRecursive(query, node_first, guess, minDist);

    double const diff = fabs(query[axis] - train[axis]);
    if(diff < *minDist){
        std::shared_ptr<Node> node_second = query[axis] < train[axis] ? node->right : node->left;
        nearestNeighborRecursive(query, node_second, guess, minDist);
    }
}
