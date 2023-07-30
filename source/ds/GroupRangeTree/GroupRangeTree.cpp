#include "GroupRangeTree.h"


struct _3DPoint
{
    using coord_type = double;

    coord_type x, y, z;

    inline _3DPoint(double x, double y, double z): x(x), y(y), z(z){};
    explicit inline _3DPoint(): _3DPoint(0, 0, 0){};
    inline _3DPoint(const _3DPoint &other){this->x = other.x; this->y = other.y; this->z = other.z;};

    inline coord_type &operator[](int idx){if(idx == 0) return this->x; if(idx == 1) return this->y; return this->z;};
    inline const coord_type &operator[](int idx) const{if(idx == 0) return this->x; if(idx == 1) return this->y; return this->z;};
    inline bool operator==(const _3DPoint &other) const{return this->x == other.x and this->y == other.y and this->z == other.z;};
    inline bool operator<=(const _3DPoint &other) const
    {
        if(this->x < other.x) return true;
        if(this->x == other.x and this->y < other.y) return true;
        if(this->x == other.x and this->y == other.y and this->z <= other.z) return true;
        return false;
    }
    friend std::ostream &operator<<(std::ostream &stream, const _3DPoint &point);
};

template<typename T, int N>
void GroupRangeTree<T, N>::recreateSubtree(GroupRangeNode *node)
{
    if(node == nullptr)
    {
        return;
    }
    if(this->currentDim == this->dim - 1)
    {
        return;
    }
    delete node->subtree;
    node->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
    node->subtree->build(this->getAllDecendants(node));

    // rebuild also my children:
    this->recreateSubtree(dynamic_cast<GroupRangeNode*>(node->left));
    this->recreateSubtree(dynamic_cast<GroupRangeNode*>(node->right));
}

template<typename T, int N>
template<typename RandomAccessIterator>
typename GroupRangeTree<T, N>::GroupRangeNode *GroupRangeTree<T, N>::fastBuildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last)
{
    if(first == last)
    {
        return nullptr;
    }

    GroupRangeNode *node;
    if(last - first <= N)
    {
        node = new GroupRangeNode();
        node->numValues = last - first;
        int j = 0;
        for(RandomAccessIterator it = first; it != last; it++)
        {
            node->values[j++] = *it;
        }
    }
    else
    {
        RandomAccessIterator mid = first + (last - first) / 2;
        node = new GroupRangeNode(*mid);
        GroupRangeNode *left = this->fastBuildHelper(first, mid);
        node->left = left;
        GroupRangeNode *right = this->fastBuildHelper(mid + 1, last); 
        node->right = right;
    }
    if(this->currentDim == this->dim - 1)
    {
        node->subtree = nullptr;
    }
    else
    {
        node->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
        node->subtree->build(this->getAllDecendants(node));
    }
    this->treeSize += node->numValues;
    
    if(node->left != nullptr) node->left->parent = node;
    if(node->right != nullptr) node->right->parent = node;

    if(node->left == nullptr and node->right == nullptr)
    {
        // a leaf, run recursive update from it
        this->recursiveUpdateNodeInfo(node);
    }
    return node;
}

template<typename T, int N>
void GroupRangeTree<T, N>::splitNode(typename GroupTree<T, N>::Node *node)
{
    // todo - call parent splitNode?
    if(node == nullptr)
    {
        return;
    }
    // assert(node->numValues > 1);
    
    std::sort(node->values.begin(), node->values.begin() + node->numValues, this->compare);
    int numLeft = node->numValues / 2;
    int numRight = node->numValues - numLeft - 1;
    std::array<T, N> left;
    std::array<T, N> right;
    for(int i = 0; i < numLeft; i++)
    {
        left[i] = node->values[i];
    }
    for(int i = numLeft + 1; i < node->numValues; i++)
    {
        right[i - (numLeft + 1)] = node->values[i];
    }

    if(numLeft > 0 and node->left == nullptr)
    {
        node->left = new GroupRangeNode(left);
        node->left->numValues = numLeft;
        node->left->parent = node;
    }
    if(numRight > 0 and node->right == nullptr)
    {
        node->right = new GroupRangeNode(right);
        node->right->numValues = numRight;
        node->right->parent = node;
    }

    node->numValues = 1;
    node->values[0] = node->values[numLeft];
    node->minInSubtree = node->maxInSubtree = node->values[0];

    this->updateNodeInfo(node->left);
    this->updateNodeInfo(node->right);
    this->recreateSubtree(dynamic_cast<GroupRangeNode*>(node)); // rebuild the sub-trees
}

template<typename T, int N>
typename GroupRangeTree<T, N>::GroupRangeNode *GroupRangeTree<T, N>::tryInsert(const T &value)
{
    if(this->getRoot() == nullptr)
    {
        this->setRoot(new GroupRangeNode(value));
        if(this->currentDim != this->dim - 1)
        {
            this->getRoot()->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
            this->getRoot()->subtree->insert(value);
        }
        return this->getRoot();
    }
    GroupRangeNode *current = this->getRoot();
    while(true)
    {
        if(this->compare(value, current->values[0]))
        {
            if(current->left == nullptr)
            {
                break;
            }
            current = dynamic_cast<GroupRangeNode*>(current->left);
        }
        else
        {
            if(current->right == nullptr)
            {
                break;
            }
            current = dynamic_cast<GroupRangeNode*>(current->right);
        }
    }

    // check if value already exists:
    for(int i = 0; i < current->numValues; i++)
    {
        if(value == current->values[i])
        {
            return current;
        }
    }

    GroupRangeNode *addTo;
    if(current->numValues < N)
    {
        addTo = current;
    }
    else
    {
        this->splitNode(current);
        if(this->compare(value, current->values[0]))
        {
            if(current->left == nullptr)
            {
                current->left = new GroupRangeNode(value);
                current->left->parent = current;
            }
            addTo = dynamic_cast<GroupRangeNode*>(current->left);
        }
        else
        {
            if(current->right == nullptr)
            {
                current->right = new GroupRangeNode(value);
                current->right->parent = current;
            }
            addTo = dynamic_cast<GroupRangeNode*>(current->right);
        }
        addTo->numValues -= 1;
    }
    addTo->values[addTo->numValues] = value;
    addTo->numValues++;
    
    // add the new value to its ancestors subtrees:
    if(this->currentDim != this->dim - 1)
    {
        GroupRangeNode *node = addTo;
        while(node != nullptr)
        {
            if(node->subtree == nullptr)
            {
                assert(node == addTo);
                // should not happen, unless `node` is `addTo` (because its ancesors' subtrees should be exist)
                node->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
            }
            node->subtree->insert(value);
            node = dynamic_cast<GroupRangeNode*>(node->parent);
        }
    }
    return addTo;
}

template<typename T, int N>
void GroupRangeTree<T, N>::rangeHelper(const std::vector<std::pair<typename T::coord_type, typename T::coord_type>> &range, const GroupRangeNode *node, int coord, std::vector<T> &result) const
{
    if(node == nullptr or coord >= this->dim)
    {
        return;
    }
    if(node->minInSubtree[coord] > range[coord].second or node->maxInSubtree[coord] < range[coord].first)
    {
        return;
    }
    if(coord == this->dim - 1)
    {
        // search over the tree to find the maching points
        this->rangeHelper(range, dynamic_cast<const GroupRangeNode*>(node->left), coord, result);
        for(int i = 0; i < node->numValues; i++)
        {
            const T &value = node->values[i];
            if(value[coord] >= range[coord].first and value[coord] <= range[coord].second)
            {
                result.push_back(value);
            }
        }
        this->rangeHelper(range, dynamic_cast<const GroupRangeNode*>(node->right), coord, result);
    }
    else
    {
        if(node->minInSubtree[coord] >= range[coord].first and node->maxInSubtree[coord] <= range[coord].second)
        {
            // the whole subtree matches this coordinate! move on to the next one
            assert(node->subtree != nullptr); // if the subtree was nullptr, then coord was dim-1.
            node->subtree->rangeHelper(range, node->subtree->getRoot(), coord + 1, result);
        }
        else
        {
            for(int i = 0; i < node->numValues; i++)
            {
                const T &value = node->values[i];
                int j;
                for(j = coord; j < this->dim; j++) if(value[j] < range[j].first or value[j] > range[j].second) break;
                if(j == this->dim) result.push_back(value);
            }
            this->rangeHelper(range, dynamic_cast<const GroupRangeNode*>(node->right), coord, result);
            this->rangeHelper(range, dynamic_cast<const GroupRangeNode*>(node->left), coord, result);
        }
    }
}
