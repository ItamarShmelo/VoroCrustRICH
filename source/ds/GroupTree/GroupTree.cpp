#include "GroupTree.h"

template<typename T, int N>
void GroupTree<T, N>::Node::getAllDecendantsHelper(std::vector<T> &vec)
{
    if(this->left != nullptr) this->left->getAllDecendantsHelper(vec);
    for(int i = 0; i < this->numValues; i++)
    {
        vec.push_back(this->values[i]);
    }
    if(this->right != nullptr) this->right->getAllDecendantsHelper(vec);
}

template<typename T, int N>
template<typename RandomAccessIterator>
typename GroupTree<T, N>::Node *GroupTree<T, N>::fastBuildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last)
{
    if(first == last)
    {
        return nullptr;
    }

    Node *node, *left = nullptr, *right = nullptr;
    if(last - first <= N)
    {
        node = new Node();
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
        node = new Node(*mid);
        Node *left = this->fastBuildHelper(first, mid);
        node->left = left;
        this->updateNodeInfo(left);
        Node *right = this->fastBuildHelper(mid + 1, last);
        node->right = right;
        this->updateNodeInfo(right);
    }
    this->treeSize += node->numValues;
    
    if(left != nullptr) left->parent = node;
    if(right != nullptr) right->parent = node;

    return node;
}

template<typename T, int N>
template<typename Container>
GroupTree<T, N>::Node::Node(const Container &container): numValues(std::min<size_t>(container.size(), N))
{
    for(int i = 0; i < this->numValues; i++)
    {
        this->values[i] = container[i];
    }
    this->left = this->right = this->parent = nullptr;
}

template<typename T, int N>
GroupTree<T, N>::Node::Node(const std::array<T, N> &values, int numValues): numValues(numValues)
{
    for(int i = 0; i < this->numValues; i++)
    {
        this->values[i] = reinterpret_cast<T>(values[i]);
    }
    this->left = this->right = this->parent = nullptr;
}

template<typename T, int N>
const typename GroupTree<T, N>::Node *GroupTree<T, N>::tryFind(const T &value) const
{
    const Node *node = this->getRoot();
    while(node != nullptr)
    {
        if(node->numValues == 1 and node->values[0] == value)
        {
            return node;
        }
        
        if(node->numValues > 1)
        {
            for(int i = 0; i < node->numValues; i++)
            {
                if(node->values[i] == value)
                {
                    return node;
                }
            }
            return nullptr; // we reached a leaf and `value` is not here
        }
        else
        {
            if(this->compare(value, node->values[0]))
            {
                node = node->left;
            }
            else
            {
                node = node->right;
            }
        }
    }
    return nullptr;
}

template<typename T, int N>
void GroupTree<T, N>::splitNode(Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    assert(node->numValues > 1);
    
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
        node->left = new Node(left);
        node->left->numValues = numLeft;
        node->left->parent = node;
    }
    if(numRight > 0 and node->right == nullptr)
    {
        node->right = new Node(right);
        node->right->numValues = numRight;
        node->right->parent = node;
    }

    node->numValues = 1;
    node->values[0] = node->values[numLeft];
    node->minInSubtree = node->maxInSubtree = node->values[0];

    this->updateNodeInfo(node->left);
    this->updateNodeInfo(node->right);
}

template<typename T, int N>
void GroupTree<T, N>::updateNodeInfo(Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    if(node->left == nullptr and node->right == nullptr)
    {
        node->minInSubtree = node->maxInSubtree = node->values[0];
        for(int i = 0; i < node->numValues; i++)
        {
            node->minInSubtree = std::min(node->minInSubtree, node->values[i], this->compare);
            node->maxInSubtree = std::max(node->maxInSubtree, node->values[i], this->compare);
        }
    }
    if(node->parent == nullptr)
    {
        node->depth = 0;
        return;
    }
    node->depth = node->parent->depth + 1;
    node->parent->height = std::max(node->parent->height, node->height + 1);
    node->parent->minInSubtree = std::min(node->minInSubtree, node->parent->minInSubtree, this->compare);
    node->parent->maxInSubtree = std::max(node->maxInSubtree, node->parent->maxInSubtree, this->compare);
}

template<typename T, int N>
bool GroupTree<T, N>::insert(const T &value)
{
    Node *node = this->tryInsert(value);
    if(node == nullptr)
    {
        return false;
    }
    ++this->treeSize;
    this->recursiveUpdateNodeInfo(node);
    return true;
}

template<typename T, int N>
typename GroupTree<T, N>::Node *GroupTree<T, N>::tryInsert(const T &value)
{
    if(this->getRoot() == nullptr)
    {
        this->setRoot(new Node(value));
        return this->getRoot();
    }
    Node *current = this->getRoot();
    while(true)
    {
        if(this->compare(value, current->values[0]))
        {
            if(current->left == nullptr)
            {
                break;
            }
            current = current->left;
        }
        else
        {
            if(current->right == nullptr)
            {
                break;
            }
            current = current->right;
        }
    }

    for(int i = 0; i < current->numValues; i++)
    {
        if(value == current->values[i])
        {
            return current;
        }
    }

    Node *addTo;
    if(current->numValues < N)
    {
        addTo = current;
    }
    else
    {
        this->splitNode(current);
        addTo = (this->compare(value, current->values[0]))? current->left : current->right;
    }
    addTo->values[addTo->numValues] = value;
    addTo->numValues++;
    return addTo;
}

#ifdef DEBUG_MODE
template<typename T, int N>
void GroupTree<T, N>::printHelper(const Node *node, int indent) const
{
    if(node == nullptr)
    {
        std::cout << "nullptr" << std::endl;
    }
    else
    {
        if(node->numValues == 1)
        {
            std::cout << node->values[0] << " (minInSubtree: " << node->minInSubtree << ", maxInSubtree: " << node->maxInSubtree << ")" << std::endl;
        }
        else
        {
            std::cout << "{";
            for(int i = 0; i < node->numValues - 1; i++)
            {
                std::cout << node->values[i] << ", ";
            }
            std::cout << node->values[node->numValues - 1] << "} " << " (minInSubtree: " << node->minInSubtree << ", maxInSubtree: " << node->maxInSubtree << ")" << std::endl;
        }
    }
    if(node != nullptr)
    {
        for(int i = 0; i < indent; i++){std::cout << "\t";}
        std::cout << "[L] ";
        printHelper(node->left, indent + 1);
        for(int i = 0; i < indent; i++){std::cout << "\t";}
        std::cout << "[R] ";
        printHelper(node->right, indent + 1);
    }
}
#endif // DEBUG_MODE