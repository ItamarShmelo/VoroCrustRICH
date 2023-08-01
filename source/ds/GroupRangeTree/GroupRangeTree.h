#ifndef _RICH_GROUP_RANGETREE_H
#define _RICH_GROUP_RANGETREE_H

#include <vector>
#include <algorithm>
#include "../GroupTree/GroupTree.h"

template<typename T, int N>
class GroupRangeTree : public GroupTree<T, N>
{
protected:
    class GroupRangeNode : public GroupTree<T, N>::Node
    {
    public:
        GroupRangeNode(): GroupTree<T, N>::Node(), subtree(nullptr){};
        GroupRangeNode(const T &value): GroupTree<T, N>::Node(value), subtree(nullptr){};
        template<typename Container>
        GroupRangeNode(const Container &container): GroupTree<T, N>::Node(container), subtree(nullptr){};
        ~GroupRangeNode() override{delete this->subtree;};

        GroupRangeTree<T, N> *subtree;
    };

    void splitNode(typename GroupTree<T, N>::Node *node) override;
    GroupRangeNode *tryInsert(const T &value);

private:
    inline typename GroupRangeTree<T, N>::GroupRangeNode *getRoot() override{return this->actualRoot;};
    inline const typename GroupRangeTree<T, N>::GroupRangeNode *getRoot() const override{return this->actualRoot;};
    inline void setRoot(typename GroupTree<T, N>::Node *other) override{if(other == nullptr) return; this->actualRoot = dynamic_cast<GroupRangeNode*>(other); assert(this->actualRoot != nullptr);};

    template<typename RandomAccessIterator>
    GroupRangeNode *fastBuildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last);
    void rangeHelper(const std::vector<std::pair<typename T::coord_type, typename T::coord_type>> &range, const GroupRangeNode *node, int coord, std::vector<T> &result) const;

    void recreateSubtree(GroupRangeNode *node);

    inline void deleteHelper(typename GroupTree<T, N>::Node *node) override{if(node == nullptr) return; deleteHelper(node->left); deleteHelper(node->right); delete dynamic_cast<GroupRangeNode*>(node);};
    inline void deleteTree() override{this->deleteHelper(this->getRoot()); this->setRoot(nullptr); this->treeSize = 0;};

public:
    GroupRangeTree(GroupRangeNode *root, const typename GroupTree<T, N>::Compare &compare, int currentDim, int dimensions): GroupTree<T, N>(root, compare), actualRoot(root), dim(dimensions), currentDim(currentDim){};
    GroupRangeTree(const typename GroupTree<T, N>::Compare &compare, int currentDim, int dimensions): GroupRangeTree(nullptr, compare, currentDim, dimensions){};
    GroupRangeTree(int currentDim, int dimensions): GroupRangeTree([currentDim](const T &a, const T &b){return (a[currentDim] <= b[currentDim]);}, currentDim, dimensions){};
    GroupRangeTree(int dimensions): GroupRangeTree<T, N>(0, dimensions){};
    template<typename InputIterator>
    inline GroupRangeTree(const InputIterator &first, const InputIterator &last, int dimensions): GroupRangeNode(){for(InputIterator it = first; it != last; it++) this->insert(*it);};
    ~GroupRangeTree() override{this->deleteTree();};

    inline std::vector<T> range(const std::vector<std::pair<typename T::coord_type, typename T::coord_type>> &range) const{std::vector<T> result; this->rangeHelper(range, this->getRoot(), 0, result); return result;};
    std::vector<T> circularRange(const T &center, typename T::coord_type radius) const; 

    template<typename RandomAccessIterator>
    inline void build(RandomAccessIterator first, RandomAccessIterator last){assert(this->treeSize == 0); std::sort(first, last, this->compare); this->setRoot(this->fastBuildHelper(first, last));};
    template<typename Container>
    inline void build(Container &&container){this->build(container.begin(), container.end());};
    template<typename Container>
    inline void build(Container &container){this->build(container.begin(), container.end());};

private:
    GroupRangeNode *actualRoot = nullptr;
    int dim, currentDim;
};

#endif // _RICH_GROUP_RANGETREE_H