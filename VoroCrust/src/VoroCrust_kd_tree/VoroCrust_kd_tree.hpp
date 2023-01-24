#ifndef VOROCRUST_KD_TREE
#define VOROCRUST_KD_TREE

#include <memory>
#include <vector>
#include <algorithm>
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"

#define DIM 3

struct Node
{
    int index;
    int axis;
    std::shared_ptr<Node> left, right;

    Node() : index(0), axis(0), left(nullptr), right(nullptr) {}
};

std::shared_ptr<Node> newNode(int index, int axis);

class VoroCrust_KD_Tree {
    public:  
        std::shared_ptr<Node> root;
        
        std::vector<Vector3D> points;

        VoroCrust_KD_Tree();

        VoroCrust_KD_Tree(std::vector<Vector3D> const& points);

        void makeTree(std::vector<Vector3D> const& points_);

        void clear();

        std::shared_ptr<Node> buildRecursive(int *indices, int npoints, int depth);

        void insert(Vector3D const& point);

        void insertRecursive(Vector3D const& point, std::shared_ptr<Node> const& node);

        void remakeTree();

        int nearestNeighbor(Vector3D const& query) const;

        void nearestNeighborRecursive(Vector3D const& query, std::shared_ptr<Node> const& node, int *guess, double *minDist) const;

        bool operator==(VoroCrust_KD_Tree const& tree) const;

        bool equalRecursive(std::shared_ptr<Node> const& node1, std::shared_ptr<Node> const& node2, VoroCrust_KD_Tree const& t) const;
};



#endif /* VOROCRUST_KD_TREE */