#ifndef VOROCRUST_KD_TREE
#define VOROCRUST_KD_TREE

#include <memory>
#include <vector>
#include <algorithm>
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"

#define DIM 3

/*! \brief a Node in a VoroCrust_KD_Tree
*/
struct Node;
using NodePtr = std::shared_ptr<Node>;

struct Node
{
    //! \brief index of point in the VoroCrust_KD_Tree::points vector
    int index;
    //! \brief depth % DIM. depth of node is depth in the tree
    int axis;
    //! \brief pointers to left and right nodes in tree
    NodePtr left, right;

    Node() : index(0), axis(0), left(nullptr), right(nullptr) {}
};

NodePtr newNode(int index, int axis);

class VoroCrust_KD_Tree {
    public:
        //! \brief pointer to root node of the tree  
        NodePtr root;
        
        //! \brief holds the actual points of the tree
        std::vector<Vector3D> points;

        //! \brief construct an empty VoroCrust_KD_Tree
        VoroCrust_KD_Tree();

        //! \brief construct a VoroCrust_KD_Tree from a given vector of points
        VoroCrust_KD_Tree(std::vector<Vector3D> const& points);
        
        //! \brief make the VoroCrust_KD_Tree from a given vector of points
        void makeTree(std::vector<Vector3D> const& points_);

        //! \brief clear and erase the tree
        void clear();

        //! \brief builds the tree recursively
        NodePtr buildRecursive(int *indices, int npoints, int depth);

        //! \brief insert a new point to the tree
        void insert(Vector3D const& point);

        //! \brief recursively checks where to insert the new point
        void insertRecursive(Vector3D const& point, NodePtr const& node);

        //! \brief remakes the tree (used after a lot of new insertions)
        void remakeTree();

        //! \brief finds nearest neighbor in tree to `query`
        int nearestNeighbor(Vector3D const& query) const;

        //! \brief finds nearest neighbor recursively in tree to `query`
        void nearestNeighborRecursive(Vector3D const& query, NodePtr const& node, int *guess, double *minDist) const;

        /*! \brief finds the `k` nearest neighbors to `query` in the tree 
            \param query
            \param k number of nearest neighbors */
        std::vector<int> kNearestNeighbors(Vector3D const& query, int const k) const;
        
        //! \brief finds the `k` nearest neighbors in the tree recursively
        void kNearestNeighborsRecursive(Vector3D const& query, int const k, NodePtr const& node, std::vector<int>& indices, std::vector<double> &minDist) const;

        std::vector<int> radiusSearch(Vector3D const& query, double const radius) const;

        void radiusSearchRecursive(Vector3D const& query, double const radius, NodePtr const& node, std::vector<int> &indices) const;

        //! \brief checks if two trees are equal
        bool operator==(VoroCrust_KD_Tree const& tree) const;
        
        //! \brief checks if two trees are equal recursively
        bool equalRecursive(NodePtr const& node1, NodePtr const& node2, VoroCrust_KD_Tree const& t) const;
};

//! \brief a boundary tree for F_C, F_E and T_S
class VoroCrust_KD_Tree_Boundary : public VoroCrust_KD_Tree {
    public:
        std::vector<Vector3D> vectors;

        VoroCrust_KD_Tree_Boundary();

        VoroCrust_KD_Tree_Boundary(std::vector<Vector3D> const& points) : VoroCrust_KD_Tree(points) {}

        VoroCrust_KD_Tree_Boundary(std::vector<Vector3D> const& points, std::vector<Vector3D> const& vecs);

        ~VoroCrust_KD_Tree_Boundary() = default;

        void insert(Vector3D const& point, Vector3D const& vec);
};

class VoroCrust_KD_Tree_Ball : public VoroCrust_KD_Tree {
    public:
        std::vector<double> ball_radii;

        VoroCrust_KD_Tree_Ball();

        VoroCrust_KD_Tree_Ball(std::vector<Vector3D> const& points, std::vector<double> const& radii);

        ~VoroCrust_KD_Tree_Ball() = default;

        void insert(Vector3D const& point, double radius);
};



#endif /* VOROCRUST_KD_TREE */