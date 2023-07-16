/**
 * Implements a d-dimensional data structure for fast range queries, called range trees.
 * See more here: https://www.cs.umd.edu/class/fall2020/cmsc420-0201/Lects/lect17-range-tree.pdf
 * 
 * @author maor miz
*/
#ifndef _RICH_RANGETREE_H_
#define _RICH_RANGETREE_H_

#include <algorithm>
#include "../BinaryTree.h"

template<typename T>
class RangeTree : public BinaryTree<T>
{
protected:
    class RangeNode : public BinaryTree<T>::Node
    {
    public:
        RangeNode(): BinaryTree<T>::Node(), subtree(nullptr){};
        RangeNode(const T &value): BinaryTree<T>::Node(value), subtree(nullptr){};
        ~RangeNode() override{delete subtree;};

        RangeTree *subtree;
    };

    inline typename RangeTree<T>::RangeNode *getRoot() override{return this->actualRoot;};
    inline const typename RangeTree<T>::RangeNode *getRoot() const override{return this->actualRoot;};
    inline void setRoot(typename BinaryTree<T>::Node *other) override{if(other == nullptr) return; this->actualRoot = dynamic_cast<RangeNode*>(other); assert(this->actualRoot != nullptr);};

public:
    RangeTree(RangeNode *root, const typename BinaryTree<T>::Compare &compare, int currentDim, int dimensions): BinaryTree<T>(root, compare), actualRoot(root), dim(dimensions), currentDim(currentDim){};
    RangeTree(const typename BinaryTree<T>::Compare &compare, int currentDim, int dimensions): RangeTree<T>(nullptr, compare, currentDim, dimensions){};
    RangeTree(int currentDim, int dimensions): RangeTree<T>([currentDim](const T &a, const T &b){return (a[currentDim] <= b[currentDim]);}, currentDim, dimensions){};
    RangeTree(int dimensions): RangeTree<T>(0, dimensions){};
    ~RangeTree() override{this->deleteSubtree(this->getRoot());};

    template<typename InputIterator>
    inline RangeTree(const InputIterator &first, const InputIterator &last, int dimensions): RangeTree<T>(){for(InputIterator it = first; it < last; it++) this->insert(*it);}
;
    template<typename RandomAccessIterator>
    inline void build(RandomAccessIterator &first, RandomAccessIterator &last){assert(this->treeSize == 0); std::sort(first, last, this->compare); this->setRoot(this->buildHelper(first, last));};

    inline std::vector<T> range(const std::vector<std::pair<typename T::coord_type, typename T::coord_type>> &range) const{std::vector<T> result; this->rangeHelper(range, this->getRoot(), 0, result); return result;};
    inline std::vector<T> circularRange(const T &center, typename T::coord_type radius) const/*{std::vector<T> result; this->circularRangeHelper(center, radius, this->getRoot(), 0, result); return result;}*/;
    std::vector<T> circularRangeRectangular(const T &center, typename T::coord_type radius) const;

private:
    RangeNode *actualRoot = nullptr;
    int dim, currentDim;

    template<typename RandomAccessIterator>
    RangeNode *buildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last);

    void rangeHelper(const std::vector<std::pair<typename T::coord_type, typename T::coord_type>> &range, const RangeNode *node, int coord, std::vector<T> &result) const;
    void circularRangeHelper(const T &center, typename T::coord_type radius, const RangeNode *node, int coord, std::vector<T> &result) const;
};

#endif // _RICH_RANGETREE_H_