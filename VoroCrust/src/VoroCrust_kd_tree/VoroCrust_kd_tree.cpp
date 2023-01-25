#include "VoroCrust_kd_tree.hpp"
#include <functional>
#include <numeric>
#include <iostream>
 
//! TODO: index should be std::size_t and everything else (except axis) which is unsigned should be std::size_t

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
    root.reset(); // root = nullptr
    points.clear();
}

void VoroCrust_KD_Tree::makeTree(std::vector<Vector3D> const& points_){
    clear();

    points = points_;
    std::vector<int> indices(points.size()); 
    std::iota(indices.begin(), indices.end(), 0); // after iota : indices[0] = 0, indices[1] = 1 etc..
    
    root = buildRecursive(indices.data(), static_cast<int>(points.size()), 0);
}

std::shared_ptr<Node> VoroCrust_KD_Tree::buildRecursive(int * indices, int npoints, int depth){
    if(npoints <= 0){
        return nullptr;
    }

    int const axis = depth % DIM; // which axis now divides the data
    int const mid  = (npoints - 1) / 2;
    
    // find the median point, nth_element also weak sorts the array around `mid` i.e.
    //  for i < mid : points[indices[i]][axis] <= points[indices[mid]][axis] and
    //  for i > mid : points[indices[mid]][axis] <= points[indices[i]][axis]
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

    // if distance to node is less the current minimal distance updae `minDist` and `guess`
    if(dist < *minDist){
        *minDist = dist;
        *guess = node->index;
    }

    int const axis = node->axis;
    // which node subtree to search tree
    std::shared_ptr<Node> node_first = query[axis] < train[axis] ? node->left : node->right;
    
    nearestNeighborRecursive(query, node_first, guess, minDist);

    double const diff = fabs(query[axis] - train[axis]);
    // if distance to current dividing axis to other node is `more` then `minDist` 
    // then we can discard this subtree 
    if(diff < *minDist){
        std::shared_ptr<Node> node_second = query[axis] < train[axis] ? node->right : node->left;
        nearestNeighborRecursive(query, node_second, guess, minDist);
    }
}

void VoroCrust_KD_Tree::insert(Vector3D const& point){
    if (root == nullptr)
    {
        root = newNode(0, 0);
        points.push_back(point);
        return;
    }
    
    insertRecursive(point, root);
}

void VoroCrust_KD_Tree::insertRecursive(Vector3D const& point, std::shared_ptr<Node> const& node){
    int const axis = node->axis;
    Vector3D const& train = points[node->index];

    // if point[axis] is on the left of node go left
    if (point[axis] < train[axis]){
        if(node->left == nullptr){ // check if left node is null ptr if so add new Node here 
            node->left = newNode(points.size(), (axis+1) % DIM);
            points.push_back(point);
            return;
        }

        insertRecursive(point, node->left);
    } else { // go right 
        if(node->right == nullptr){  // check if right node is null ptr if so add new Node here 
            node->right = newNode(points.size(), (axis+1) % DIM);
            points.push_back(point);
            return;
        }

        insertRecursive(point, node->right);
    }
}

void VoroCrust_KD_Tree::remakeTree(){
    root.reset(); // root = nullptr

    std::vector<int> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    root = buildRecursive(indices.data(), static_cast<int>(points.size()), 0);
}

bool VoroCrust_KD_Tree::operator==(VoroCrust_KD_Tree const& t) const {
    if(points.size() != t.points.size()) return false;
    return equalRecursive(root, t.root, t);
}

bool VoroCrust_KD_Tree::equalRecursive(std::shared_ptr<Node> const& node1, std::shared_ptr<Node> const& node2, VoroCrust_KD_Tree const& t) const {
    if(node1 == nullptr || node2 == nullptr){
        if(node1 == nullptr && node2 == nullptr) // assert that if one node is nullptr then both are
            return true;
        return false;
    }

    if(points[node1->index] == t.points[node2->index]){ // check by value
        return (equalRecursive(node1->left, node2->left, t) && (equalRecursive(node1->right, node2->right, t)));
    }

    return false;
}

/* BOUNDARY KD TREE */
VoroCrust_KD_Tree_Boundary::VoroCrust_KD_Tree_Boundary() : VoroCrust_KD_Tree(), vectors() {}

VoroCrust_KD_Tree_Boundary::VoroCrust_KD_Tree_Boundary(std::vector<Vector3D> const& points, std::vector<Vector3D> const& vecs) : VoroCrust_KD_Tree(points), vectors(vecs){
    if (points.size() != vecs.size()){
        std::cout << "ERROR : points.size() != vecs.size() in initialization of VoroCrust_KD_Tree_Boundary" << std::endl;
        exit(1);
    }
}

void VoroCrust_KD_Tree_Boundary::insert(Vector3D const& point, Vector3D const& vec){
    this->VoroCrust_KD_Tree::insert(point);
    vectors.push_back(vec);
}

/* BALL KD TREE */
VoroCrust_KD_Tree_Ball::VoroCrust_KD_Tree_Ball() : VoroCrust_KD_Tree(), ball_radii() {}

VoroCrust_KD_Tree_Ball::VoroCrust_KD_Tree_Ball(std::vector<Vector3D> const& points, std::vector<double> const& radii) : VoroCrust_KD_Tree(points), ball_radii(radii) {
    if (points.size() != ball_radii.size()){
        std::cout << "ERROR : points.size() !=  ball_redii.size() in initialization of VoroCrust_KD_Tree_Ball" << std::endl;
        exit(1);
    }
}

void VoroCrust_KD_Tree_Ball::insert(Vector3D const& point, double radius){
    this->VoroCrust_KD_Tree::insert(point);
    ball_radii.push_back(radius);
}
