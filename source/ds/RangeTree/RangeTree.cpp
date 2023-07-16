#include "RangeTree.h"

template<typename T>
template<typename RandomAccessIterator>
typename RangeTree<T>::RangeNode *RangeTree<T>::buildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last)
{
    if(first == last)
    {
        return nullptr;
    }
    RandomAccessIterator mid = first + (last - first) / 2;
    RangeNode *node = new RangeNode(*mid);
    ++this->treeSize;
    RangeNode *left = this->buildHelper(first, mid);
    node->left = left;
    RangeNode *right = this->buildHelper(mid + 1, last);
    node->right = right;
    
    if(left != nullptr) left->parent = node;
    if(right != nullptr) right->parent = node;

    this->updateNodeInfo(node);
    if(this->currentDim == this->dim - 1)
    {
        node->subtree = nullptr;
    }
    else
    {
        node->subtree = new RangeTree<T>(this->currentDim + 1, this->dim);
        node->subtree->build(this->getAllDecendants(node));
    }
    return node;
}

template<typename T>
void RangeTree<T>::rangeHelper(const std::vector<std::pair<typename T::coord_type, typename T::coord_type>> &range, const RangeNode *node, int coord, std::vector<T> &result) const
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
        this->rangeHelper(range, dynamic_cast<const RangeNode*>(node->left), coord, result);
        if(node->value[coord] >= range[coord].first and node->value[coord] <= range[coord].second)
        {
            result.push_back(node->value);
        }
        this->rangeHelper(range, dynamic_cast<const RangeNode*>(node->right), coord, result);
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
            int i;
            for(i = coord; i < this->dim; i++) if(node->value[i] < range[i].first or node->value[i] > range[i].second) break;
            if(i == this->dim) result.push_back(node->value);
            this->rangeHelper(range, dynamic_cast<const RangeNode*>(node->right), coord, result);
            this->rangeHelper(range, dynamic_cast<const RangeNode*>(node->left), coord, result);
        }
    }
}

template<typename T>
void RangeTree<T>::circularRangeHelper(const T &center, typename T::coord_type radius, const RangeNode *node, int coord, std::vector<T> &result) const
{
    if(node == nullptr or coord >= this->dim)
    {
        return;
    }
    
    typename T::coord_type cornerMax = center[coord] + radius;
    typename T::coord_type cornerMin = center[coord] - radius;

    if(node->minInSubtree[coord] > cornerMax or node->maxInSubtree[coord] < cornerMin)
    {
        return;
    }

    if(coord == this->dim - 1)
    {
        typename T::coord_type distanceSquared = typename T::coord_type();
        for(int i = 0; i < this->dim; i++)
        {
            distanceSquared += ((node->value[i] - center[i]) * (node->value[i] - center[i]));
        }
        if(distanceSquared <= radius * radius)
        {
            result.push_back(node->value);
        }
        if(node->minInSubtree[coord] <= cornerMax)
        {
            this->circularRangeHelper(center, radius, dynamic_cast<const RangeNode*>(node->left), coord, result);
        }
        if(node->maxInSubtree[coord] >= cornerMin)
        {
            this->circularRangeHelper(center, radius, dynamic_cast<const RangeNode*>(node->right), coord, result);
        }
    }
    else
    {
        if(node->minInSubtree[coord] >= cornerMin and node->maxInSubtree[coord] <= cornerMax)
        {
            // the whole subtree matches this coordinate! move on to the next one
            assert(node->subtree != nullptr); // if the subtree was nullptr, then coord was dim-1.
            node->subtree->circularRangeHelper(center, radius, node->subtree->getRoot(), coord + 1, result);
        }
        else
        {
            typename T::coord_type distanceSquared = typename T::coord_type();
            for(int i = 0; i < this->dim; i++)
            {
                distanceSquared += ((node->value[i] - center[i]) * (node->value[i] - center[i]));
            }
            if(distanceSquared <= radius * radius)
            {
                result.push_back(node->value);
            }
            this->circularRangeHelper(center, radius, dynamic_cast<const RangeNode*>(node->right), coord, result);
            this->circularRangeHelper(center, radius, dynamic_cast<const RangeNode*>(node->left), coord, result);
        }
    }
}

template<typename T>
std::vector<T> RangeTree<T>::circularRangeRectangular(const T &point, typename T::coord_type radius) const
{
    // can have better implementations in the future...
    std::vector<std::pair<typename T::coord_type, typename T::coord_type>> range;
    for(int i = 0; i < this->dim; i++)
    {
        range.push_back(std::pair<typename T::coord_type, typename T::coord_type>(point[i] - radius, point[i] + radius));
    }

    std::vector<T> result;
    for(const T &possiblePoint : this->range(range))
    {
        typename T::coord_type distanceSquared = typename T::coord_type();
        for(int i = 0; i < this->dim; i++)
        {
            distanceSquared += ((possiblePoint[i] - point[i]) * (possiblePoint[i] - point[i]));
        }
        if(distanceSquared <= radius * radius)
        {
            result.push_back(possiblePoint);
        } 
    }
    return result;
}

template<typename T>
inline std::vector<T> RangeTree<T>::circularRange(const T &center, typename T::coord_type radius) const
{
    // a simple implementation - 
    std::vector<T> result;
    std::vector<std::pair<typename T::coord_type, typename T::coord_type>> rectangle;
    for(int i = 0; i < this->dim; i++)
    {
        rectangle.push_back(std::make_pair<typename T::coord_type, typename T::coord_type>(center[i] - radius, center[i] + radius));
    }

    for(const auto &point : this->range(rectangle))
    {
        typename T::coord_type squaredDistance = typename T::coord_type();
        for(int i = 0; i < this->dim; i++)
        {
            squaredDistance += (center[i] - point[i]) * (center[i] - point[i]);
        }
        if(squaredDistance <= radius * radius)
        {
            result.push_back(point);
        }
    }
    return result;
}