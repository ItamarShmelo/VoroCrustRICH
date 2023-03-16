#include "VoroCrust_kd_tree.hpp"
#include <functional>
#include <numeric>
#include <iostream>
 
//! TODO: index should be std::size_t and everything else (except axis) which is unsigned should be std::size_t

NodePtr newNode(int index, int axis){
    NodePtr node = std::make_shared<Node>();

    node->index = index;
    node->axis = axis;

    return node;
}

VoroCrust_KD_Tree::VoroCrust_KD_Tree() : root(nullptr), points() {}

VoroCrust_KD_Tree::VoroCrust_KD_Tree(std::vector<Vector3D> const& points) : root(nullptr), points(points) {
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

NodePtr VoroCrust_KD_Tree::buildRecursive(int * indices, int npoints, int depth){
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

    NodePtr node = newNode(indices[mid], axis);

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

void VoroCrust_KD_Tree::nearestNeighborRecursive(Vector3D const& query, NodePtr const& node, int *guess, double *minDist) const {
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
    NodePtr const& node_first = query[axis] < train[axis] ? node->left : node->right;
    
    nearestNeighborRecursive(query, node_first, guess, minDist);

    double const diff = fabs(query[axis] - train[axis]);
    // if distance to current dividing axis to other node is `more` then `minDist` 
    // then we can discard this subtree 
    if(diff < *minDist){
        NodePtr const& node_second = query[axis] < train[axis] ? node->right : node->left;
        nearestNeighborRecursive(query, node_second, guess, minDist);
    }
}

std::vector<int> VoroCrust_KD_Tree::kNearestNeighbors(Vector3D const& query, unsigned int const k) const {
    if(k > points.size()){
        std::cout << "Error: kNearestNeighbors, `k` can't be larger then number of points" << std::endl;
        exit(1);
    }

    std::vector<double> minDist(k, std::numeric_limits<double>::max());

    std::vector<int> indices(k, -1);

    kNearestNeighborsRecursive(query, k, root, indices, minDist);
    return indices;
}

void VoroCrust_KD_Tree::kNearestNeighborsRecursive(Vector3D const& query, int const k, NodePtr const& node, std::vector<int>& indices, std::vector<double> &minDist) const{
    if(node.get() == nullptr) return;

    Vector3D const& train = points[node->index];

    double const dist = distance(train, query);
    int i;
    for(i = k; i > 0; --i){
        if(dist > minDist[i-1])
            break;
    }

    if(i != k){
        indices.insert(indices.begin() + i, node->index);
        minDist.insert(minDist.begin() + i, dist);

        indices.pop_back();
        minDist.pop_back();
    }

    int const axis = node->axis;
    NodePtr node_first = query[axis] < train[axis] ? node->left : node->right;

    kNearestNeighborsRecursive(query, k, node_first, indices, minDist);

    double const diff = fabs(query[axis] - train[axis]);
    //! TODO: make the epsilon a constant throughout the program or use machine limits
    if(diff < minDist.back()){
        NodePtr const& node_second = query[axis] < train[axis] ? node->right : node->left;
        kNearestNeighborsRecursive(query, k, node_second, indices, minDist);
    }
}

std::vector<int> VoroCrust_KD_Tree::radiusSearch(Vector3D const& query, double const radius) const {
    std::vector<int> indices;
    radiusSearchRecursive(query, radius, root, indices);
    return indices;
}

void VoroCrust_KD_Tree::radiusSearchRecursive(Vector3D const& query, double const radius, NodePtr const& node, std::vector<int> &indices) const {
    if(node.get() == nullptr) return;

    Vector3D const& train = points[node->index];
    
    double const dist = distance(query, train);

    if(dist < radius) indices.push_back(node->index);

    int const axis = node->axis;
    NodePtr node_first = query[axis] < train[axis] ? node->left : node->right;

    radiusSearchRecursive(query, radius, node_first, indices);

    double const diff = fabs(query[axis] - train[axis]);
    
    if (diff < radius){
        NodePtr const& node_second = query[axis] < train[axis] ? node->right : node->left;
        radiusSearchRecursive(query, radius, node_second, indices);
    }
}

int VoroCrust_KD_Tree::nearestNeighborToSegment(std::array<Vector3D, 2> const& segment) const {
    int guess = 0;
    double minDist = std::numeric_limits<double>::max();

    nearestNeighborToSegmentRecursive(segment, root, guess, minDist);
    return guess;
}

void VoroCrust_KD_Tree::nearestNeighborToSegmentRecursive(std::array<Vector3D, 2> const& segment, NodePtr const& node, int &guess, double &minDist) const{
    if(node.get() == nullptr) return;

    Vector3D const& train = points[node->index];

    double const dist = distancePointToSegment(segment, train);

    if(dist < minDist){
        minDist = dist;
        guess = node->index;
    }

    int const axis = node->axis;

    Vector3D const& point0 = segment[0];
    Vector3D const& point1 = segment[1];

    double const diff = std::min({fabs(point0[axis] - train[axis]), fabs(point1[axis] - train[axis])});

    if(point0[axis] <= train[axis] && point1[axis] <= train[axis]){
        nearestNeighborToSegmentRecursive(segment, node->left, guess, minDist);

        if(diff < minDist){
            nearestNeighborToSegmentRecursive(segment, node->right, guess, minDist);
        }
    } else if(point0[axis] >= train[axis] && point1[axis] >= train[axis]){
        nearestNeighborToSegmentRecursive(segment, node->right, guess, minDist);

        if(diff < minDist){
            nearestNeighborToSegmentRecursive(segment, node->left, guess, minDist);
        }
    } else {
        nearestNeighborToSegmentRecursive(segment, node->left, guess, minDist);
        nearestNeighborToSegmentRecursive(segment, node->right, guess, minDist);
    }
}

double VoroCrust_KD_Tree::distancePointToSegment(std::array<Vector3D, 2> const& segment, Vector3D const& point) const {
    // taken from https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

    Vector3D const& v1 = segment[0] - point; // x1-p
    Vector3D const& v2 = segment[1] - segment[0]; // x2-x1

    // proj_{segment}(p) = x1 + (x2-x1)*t
    double const t = - ScalarProd(v1, v2) / ScalarProd(v2, v2);

    // if t >= 1 than the projection is after x2 so return the distance to x2
    if(t >= 1){
        return distance(point, segment[1]);
    }
    
    // if t <= 0 than the projection is before x1 so return the distance to x1
    if(t <= 0){
        return distance(point, segment[0]);
    }
    
    // if 0<t<1 then the projection is on the segment, return the abs(p-proj(p)) (i.e the magnitude of the perpinicular part)
    return abs(CrossProduct(v1, v2))/ abs(v2);
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

void VoroCrust_KD_Tree::insertRecursive(Vector3D const& point, NodePtr const& node){
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

bool VoroCrust_KD_Tree::equalRecursive(NodePtr const& node1, NodePtr const& node2, VoroCrust_KD_Tree const& t) const {
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

VoroCrust_KD_Tree_Boundary::VoroCrust_KD_Tree_Boundary(std::vector<Vector3D> const& points, std::vector<Vector3D> const& vecs, std::vector<std::size_t> const& feature_index_, std::vector<std::size_t> const& plc_index_) : VoroCrust_KD_Tree(points), vectors(vecs), feature_index(feature_index_), plc_index(plc_index_) {
    if (points.size() != vecs.size()){
        std::cout << "ERROR : points.size() != vecs.size() in initialization of VoroCrust_KD_Tree_Boundary" << std::endl;
        exit(1);
    }

    if (points.size() != feature_index.size()){
        std::cout << "ERROR : points.size() !=  feature_index.size() in initialization of VoroCrust_KD_Tree_Ball" << std::endl;
        exit(1);
    }
    
    if (points.size() != plc_index.size()){
        std::cout << "ERROR : points.size() !=  plc_index.size() in initialization of VoroCrust_KD_Tree_Ball" << std::endl;
        exit(1);
    }
}

void VoroCrust_KD_Tree_Boundary::insert(Vector3D const& point, Vector3D const& vec, std::size_t const f_index, std::size_t const plc_index_){
    this->VoroCrust_KD_Tree::insert(point);
    vectors.push_back(vec);
    feature_index.push_back(f_index);
    plc_index.push_back(plc_index_);
}

int VoroCrust_KD_Tree_Boundary::nearestNonCosmoothPointEdge(Vector3D const& query, Vector3D const& vec, std::size_t const f_index, double const angle) const {
    int guess = -1;
    double minDist = std::numeric_limits<double>::max();

    nearestNonCosmoothPointEdgeRecursive(query, vec, f_index, angle, root, guess, minDist);

    return guess;
}


//! \brief calculates the angles assuming v1, v2 have the same orientation
double CalcAngleAssumedSameOrientation(Vector3D const& v1, Vector3D const& v2){
    return acos(std::abs(ScalarProd(v1, v2)) / abs(v1) / abs(v2));
}

void VoroCrust_KD_Tree_Boundary::nearestNonCosmoothPointEdgeRecursive(Vector3D const& query, Vector3D const& vec, std::size_t const f_index, double const angle, NodePtr const& node, int &guess, double &minDist) const {
    if(node.get() == nullptr) return;

    Vector3D const& train = points[node->index];

    if(f_index == feature_index[node->index]){    
        double const dist = distance(train, query);
    
        if(dist < minDist){
            Vector3D const& v_pq = query - train;
            Vector3D const& v_sigma_q = vectors[node->index];
            Vector3D const& v_tau_p = vec;

            //! FORDEBUGGING: should be in the if 
            // double const angle_v_pq_v_sigma_q = CalcAngleAssumedSameOrientation(v_sigma_q, v_pq);
            double const angle_v_sigma_q_v_tau_p = CalcAngle(v_sigma_q, v_tau_p);
            
            // cosmoothness test
            if(angle_v_sigma_q_v_tau_p > angle){
                minDist = dist;
                guess = node->index;
            }
        }
    }

    int const axis = node->axis;

    NodePtr const& node_first = query[axis] < train[axis] ? node->left : node->right;

    nearestNonCosmoothPointEdgeRecursive(query, vec, f_index, angle, node_first, guess, minDist);
    
    double const diff = fabs(query[axis] - train[axis]);
    // if distance to current dividing axis to other node is `more` then `minDist` 
    // then we can discard this subtree 
    if(diff < minDist){
        NodePtr const& node_second = query[axis] < train[axis] ? node->right : node->left;
        nearestNonCosmoothPointEdgeRecursive(query, vec, f_index, angle, node_second, guess, minDist);
    }
}

int VoroCrust_KD_Tree_Boundary::nearestNonCosmoothPointFace(Vector3D const& query, Vector3D const& vec, std::size_t const f_index, double const angle1) const {
    int guess = -1;
    double minDist = std::numeric_limits<double>::max();

    nearestNonCosmoothPointFaceRecursive(query, vec, f_index, angle1, root, guess, minDist);
    
    return guess;
}


//! CODEDUPLICATION: This is same exact as nearestNonCosmoothPointEdge right now!!!!
void VoroCrust_KD_Tree_Boundary::nearestNonCosmoothPointFaceRecursive(Vector3D const& query, Vector3D const& vec, std::size_t const f_index, double const angle1, NodePtr const& node, int &guess, double &minDist) const {

    if(node.get() == nullptr) return;

    int const index = node->index;
    Vector3D const& train = points[index];

    if(f_index == feature_index[index]){
        double const dist = distance(train, query);

        if(dist < minDist){
            Vector3D const& v_pq = query - train;
            Vector3D const& v_sigma_q = vectors[node->index];
            Vector3D const& v_tau_p = vec;
            
            //! FORDEBUGGING: should be in the if 
            // double const angle_v_pq_v_sigma_q = CalcAngleAssumedSameOrientation(v_sigma_q, v_pq);
            double const angle_v_sigma_q_v_tau_p = CalcAngle(v_sigma_q, v_tau_p);

            // cosmoothness test
            if(angle_v_sigma_q_v_tau_p > angle1){
                minDist = dist;
                guess = node->index;
            }
        }
    }


}


int VoroCrust_KD_Tree_Boundary::nearestNeighborExcludingFeatures(Vector3D const& query, std::vector<std::size_t> const& to_exclude) const {
    int guess = -1;
    double minDist = std::numeric_limits<double>::max();

    nearestNeighborExcludingFeaturesRecursive(query, to_exclude, root, guess, minDist);

    return guess;
}

void VoroCrust_KD_Tree_Boundary::nearestNeighborExcludingFeaturesRecursive(Vector3D const& query, std::vector<std::size_t> const& to_exclude, NodePtr const& node, int &guess, double &minDist) const {
    
    if(node.get() == nullptr) return;

    int const index = node->index;
    Vector3D const& train = points[index];

    if(std::find(to_exclude.begin(), to_exclude.end(), feature_index[index]) == to_exclude.end()){
        double const dist = distance(train, query);

        if(dist < minDist){
            minDist = dist;
            guess = index;
        }
    }

    int const axis = node->axis;

    NodePtr const& node_first = query[axis] < train[axis] ? node->left : node->right;
    
    nearestNeighborExcludingFeaturesRecursive(query, to_exclude, node_first, guess, minDist);
    
    double const diff = fabs(query[axis] - train[axis]);
    // if distance to current dividing axis to other node is `more` then `minDist` 
    // then we can discard this subtree 
    if(diff < minDist){
        NodePtr const& node_second = query[axis] < train[axis] ? node->right : node->left;
        nearestNeighborExcludingFeaturesRecursive(query, to_exclude, node_second, guess, minDist);
    }
}

/* BALL KD TREE */
VoroCrust_KD_Tree_Ball::VoroCrust_KD_Tree_Ball() : VoroCrust_KD_Tree_Boundary(), ball_radii() {}

VoroCrust_KD_Tree_Ball::VoroCrust_KD_Tree_Ball(std::vector<Vector3D> const& points, std::vector<Vector3D> const& vecs, std::vector<double> const& radii, std::vector<std::size_t> const& feature_index_,
std::vector<std::size_t> const& plc_index_) : VoroCrust_KD_Tree_Boundary(points, vecs, feature_index_, plc_index_), ball_radii(radii) {
    if (points.size() != ball_radii.size()){
        std::cout << "ERROR : points.size() !=  ball_redii.size() in initialization of VoroCrust_KD_Tree_Ball" << std::endl;
        exit(1);
    }
}

void VoroCrust_KD_Tree_Ball::insert(Vector3D const& point, Vector3D const& vec, double radius, std::size_t const f_index, std::size_t const& plc_index_){
    this->VoroCrust_KD_Tree_Boundary::insert(point, vec, f_index, plc_index_);
    ball_radii.push_back(radius);
}

std::vector<std::size_t> VoroCrust_KD_Tree_Ball::getOverlappingBalls(Vector3D const& center, double const radius, double const r_max) const {
    
    std::vector<int> suspects = radiusSearch(center, r_max);

    std::vector<std::size_t> overlapping_balls_indices;

    for(std::size_t const i : suspects){
        if(distance(center, points[i]) < radius + ball_radii[i]){
            overlapping_balls_indices.push_back(i);
        }
    }

    return overlapping_balls_indices;
}

