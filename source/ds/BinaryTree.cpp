#include "BinaryTree.h"

template<typename T>
void BinaryTree<T>::deleteSubtree(Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    this->deleteSubtree(node->left);
    node->left = nullptr;
    this->deleteSubtree(node->right);
    node->right = nullptr;
    delete node;
}

template<typename T>
const typename BinaryTree<T>::Node *BinaryTree<T>::findParent(const T &value) const
{
    const Node *current = this->getRoot();
    if(current->value == value)
    {
        return nullptr;
    }
    while(current != nullptr)
    {
        if(this->compare(value, current->value))
        {
            if(current->left == nullptr or current->left->value == value)
            {
                return current;
            }
            current = current->left;
        }
        else
        {
            if(current->right == nullptr or current->right->value == value)
            {
                return current;
            }
            current = current->right;
        }
    }
    return nullptr;
}

template<typename T>
const typename BinaryTree<T>::Node *BinaryTree<T>::findNode(const T &value) const
{
    if(this->getRoot() == nullptr or (this->getRoot() != nullptr and value == this->getRoot()->value))
    {
        return this->getRoot();
    }
    const Node *parent = this->findParent(value);
    if(parent == nullptr)
    {
        return nullptr;
    }
    if(parent->left != nullptr and parent->left->value == value)
    {
        return parent->left;
    }
    if(parent->right != nullptr and parent->right->value == value)
    {
        return parent->right;
    }
    return nullptr;
}


template<typename T>
const typename BinaryTree<T>::Node *BinaryTree<T>::findClosestAbove(const T &value) const
{
    if(this->getRoot() != nullptr and this->getRoot()->value == value)
    {
        return this->getRoot();
    }
    const Node *current = this->findParent(value);
    if(current->right != nullptr and current->right->value == value)
    {
        return current->right;
    }
    if(current->left != nullptr and current->left->value == value)
    {
        return current->left;
    }
    if(this->compare(value, current->value))
    {
        // value < current->value
        return current;
    }
    while(current != nullptr)
    {
        if(current->parent != nullptr and current == current->parent->left)
        {
            // find where for the first time, the node is a left child of its parent
            return current->parent;
        }
        current = current->parent;
    }
    return nullptr;
}

template<typename T>
const typename BinaryTree<T>::Node *BinaryTree<T>::findSuccessor(const T &value) const
{
    const Node *closest = this->findClosestBelow(value);
    if(closest == nullptr)
    {
        return nullptr;
    }
    if(!(closest->value == value))
    {
        return closest; // value is not inside the tree anyway
    }
    // closest contains the value of value
    if(closest->right != nullptr)
    {
        const Node *node = closest->right;
        while(node->left != nullptr)
        {
            node = node->left;
        }
        return node;
    }
    // closest->right is null, we need to go up
    const Node *node = closest;
    while(node->parent != nullptr and node == node->parent->right) node = node->parent;
    return node->parent;
}

template<typename T>
const typename BinaryTree<T>::Node *BinaryTree<T>::findClosestBelow(const T &value) const
{
    if(this->getRoot() != nullptr and this->getRoot()->value == value)
    {
        return this->getRoot();
    }
    const Node *current = this->findParent(value);
    if(current->right != nullptr and current->right->value == value)
    {
        return current->right;
    }
    if(current->left != nullptr and current->left->value == value)
    {
        return current->left;
    }
    if(!this->compare(value, current->value))
    {
        // value > current->value
        return current;
    }
    while(current != nullptr)
    {
        if(current->parent != nullptr and current == current->parent->right)
        {
            // find where for the first time, the node is a right child of its parent
            return current->parent;
        }
        current = current->parent;
    }
    return nullptr;
}

template<typename T>
const typename BinaryTree<T>::Node *BinaryTree<T>::findPredecessor(const T &value) const
{
    const Node *closest = this->findClosestAbove(value);
    if(closest == nullptr)
    {
        return nullptr;
    }
    if(!(closest->value == value))
    {
        return closest; // value is not inside the tree anyway
    }
    // closest contains the value of value
    if(closest->left != nullptr)
    {
        const Node *node = closest->left;
        while(node->right != nullptr)
        {
            node = node->right;
        }
        return node;
    }
    // closest->right is null, we need to go up
    const Node *node = closest;
    while(node->parent != nullptr and node == node->parent->left) node = node->parent;
    return node->parent;
}


template<typename T>
typename BinaryTree<T>::Node *BinaryTree<T>::_insert(const T &value)
{
    if(this->getRoot() == nullptr)
    {
        this->setRoot(new Node(value));
        return this->getRoot();
    }
    Node *parent = this->findParent(value);    
    if(parent == nullptr and this->getRoot()->value == value)
    {
        ++this->getRoot()->duplications;
        // this->updateNodeInfo(this->getRoot());
        return this->getRoot();
    }
    assert(parent != nullptr);

    if(parent->left != nullptr and parent->left->value == value)
    {
        ++parent->left->duplications;
        // this->updateNodeInfo(parent->left);
        return parent->left;
    }
    if(parent->right != nullptr and parent->right->value == value)
    {
        ++parent->right->duplications;
        // this->updateNodeInfo(parent->right);
        return parent->right;
    }
    assert(parent->left == nullptr or parent->right == nullptr);

    Node *node = new Node(value);
    node->parent = parent;

    if(this->compare(value, parent->value))
    {
        parent->left = node; // should be hang on left
    }
    else
    {
        parent->right = node; // should be hang on right
    }
    return node;
}

template<typename T>
bool BinaryTree<T>::insert(const T &value)
{
    Node *node = this->_insert(value);
    if(node == nullptr)
    {
        return false;
    }
    ++this->treeSize;
    this->recursiveUpdateNodeInfo(node);
    return true;
}

template<typename T>
typename BinaryTree<T>::Node *BinaryTree<T>::_remove(Node *node)
{
    if(node->left == nullptr or node->right == nullptr)
    {
        Node *probablyNotNullChild = (node->left == nullptr)? node->right : node->left;
        if(node->right != nullptr) node->right->parent = node->parent;
        if(node->left != nullptr) node->left->parent = node->parent;
        if(node->parent == nullptr)
        {
            // tree will now be empty
            this->setRoot(probablyNotNullChild);
            return node;
        }
        if(node == node->parent->left) node->parent->left = probablyNotNullChild;
        else if(node == node->parent->right) node->parent->right = probablyNotNullChild;
        return node;
    }
    Node *successor = this->findSuccessor(node->value);
    assert(successor != nullptr); // successor can not be null since the right child exists
    // replace value, and duplications
    node->value = successor->value;
    node->duplications = successor->duplications;
    Node *newPlace = _remove(successor);
    return newPlace;
}

template<typename T>
bool BinaryTree<T>::remove(const T &value)
{
    Node *node = this->findNode(value);
    if(node == nullptr)
    {
        return false; // value is not in the tree
    }
    --this->treeSize;
    if(node->duplications > 1)
    {
        --node->duplications;
        return true;
    }
    Node *fixFrom = this->_remove(node);
    assert(node != nullptr);
    Node *fixFromParent = fixFrom->parent;
    this->recursiveUpdateNodeInfo(fixFromParent);
    delete fixFrom;
    if(this->treeSize == 0)
    {
        this->setRoot(nullptr);
    }
    return true;
}

template<typename T>
void BinaryTree<T>::updateNodeInfo(Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    int myLeftHeight = (node->left == nullptr)? -1 : node->left->height;
    int myRightHeight = (node->right == nullptr)? -1 : node->right->height;
    node->height = std::max(myLeftHeight, myRightHeight) + 1;
    node->leftSize = (node->left == nullptr)? 0 : node->left->leftSize + node->left->rightSize + node->left->duplications;
    node->rightSize = (node->right == nullptr)? 0 : node->right->leftSize + node->right->rightSize + node->right->duplications;
    node->minInSubtree = (node->left == nullptr)? node->value : node->left->minInSubtree;
    node->maxInSubtree = (node->right == nullptr)? node->value : node->right->maxInSubtree;
}

template<typename T>
void BinaryTree<T>::fastRotateRight(Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    Node *newParent = node->left;
    assert(newParent != nullptr); // node->left can not be null if we want a left rotation!
    newParent->parent = node->parent;
    if(node->parent != nullptr)
    {
        if(node == node->parent->left) node->parent->left = newParent;
        else node->parent->right = newParent;
    }
    else
    {
        assert(node == this->getRoot()); // if node's parent is null, then node should the root
        this->setRoot(newParent);
    }
    if(newParent->right != nullptr) newParent->right->parent = node;
    node->left = newParent->right;
    node->parent = newParent;
    newParent->right = node;

    this->updateNodeInfo(node);
    this->updateNodeInfo(newParent);

    // int newParentLeftHeight = (newParent->left == nullptr)? -1 : newParent->left->height;
    // int newParentRightHeight = (newParent->right == nullptr)? -1 : newParent->right->height;
    // newParent->height = std::max(newParentLeftHeight, newParentRightHeight) + 1;
    // newParent->leftSize = (newParent->left == nullptr)? 0 : newParent->left->leftSize + newParent->left->rightSize + newParent->left->duplications;
    // newParent->rightSize = (newParent->right == nullptr)? 0 : newParent->right->leftSize + newParent->right->rightSize + newParent->right->duplications;
}

template<typename T>
void BinaryTree<T>::fastRotateLeft(Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    Node *newParent = node->right;
    assert(newParent != nullptr); // node->right can not be null if we want a left rotation!
    newParent->parent = node->parent;
    if(node->parent != nullptr)
    {
        if(node == node->parent->left) node->parent->left = newParent;
        else node->parent->right = newParent;
    }
    else
    {
        assert(node == this->getRoot()); // if node's parent is null, then node should the root
        this->setRoot(newParent);
    }
    if(newParent->left != nullptr) newParent->left->parent = node;
    node->right = newParent->left;
    node->parent = newParent;
    newParent->left = node;

    this->updateNodeInfo(node);
    this->updateNodeInfo(newParent);

    // int myLeftHeight = (node->left == nullptr)? -1 : node->left->height;
    // int myRightHeight = (node->right == nullptr)? -1 : node->right->height;
    // node->height = std::max(myLeftHeight, myRightHeight) + 1;
    // node->leftSize = (node->left == nullptr)? 0 : node->left->leftSize + node->left->rightSize + node->left->duplications;
    // node->rightSize = (node->right == nullptr)? 0 : node->right->leftSize + node->right->rightSize + node->right->duplications;

    // int newParentLeftHeight = (newParent->left == nullptr)? -1 : newParent->left->height;
    // int newParentRightHeight = (newParent->right == nullptr)? -1 : newParent->right->height;
    // newParent->height = std::max(newParentLeftHeight, newParentRightHeight) + 1;
    // newParent->leftSize = (newParent->left == nullptr)? 0 : newParent->left->leftSize + newParent->left->rightSize + newParent->left->duplications;
    // newParent->rightSize = (newParent->right == nullptr)? 0 : newParent->right->leftSize + newParent->right->rightSize + newParent->right->duplications;
}

template<typename T>
void BinaryTree<T>::getAllDecendantsHelper(Node *node, std::vector<T> &vec)
{
    if(node == nullptr)
    {
        return;
    }
    this->getAllDecendantsHelper(node->left, vec);
    vec.push_back(node->value);
    this->getAllDecendantsHelper(node->right, vec);
}

#ifdef DEBUG_MODE
template<typename T>
void BinaryTree<T>::printHelper(const Node *node, int tabs) const
{
    if(node == nullptr)
    {
        for(int i = 0; i < tabs; i++){std::cout << "\t";}
        std::cout << "nullptr" << std::endl;
        return;
    }
    for(int i = 0; i < tabs; i++){std::cout << "\t";}
    std::cout << node->value << " (height: " << node->height << ", duplications: " << node->duplications << ", minSubtree: " << node->minInSubtree << ", maxSubtree: " <<
            node->maxInSubtree << ", right-height: " << ((node->right == nullptr)? 0 : node->right->height) << ", left-height: " << ((node->left == nullptr)? 0 : node->left->height) << ", height: " << node->height << ", left-size: " << node->leftSize << ", right-size: " << node->rightSize << ")" << std::endl;
    for(int i = 0; i < tabs; i++){std::cout << "\t";}
    std::cout << "[L]" << std::endl;
    printHelper(node->left, tabs + 1);
    for(int i = 0; i < tabs; i++){std::cout << "\t";}
    std::cout << "[R]" << std::endl;
    printHelper(node->right, tabs + 1);
}

template<typename T>
long int BinaryTree<T>::validateHelper(const Node *node) const
{
    if(node == nullptr)
    {
        return 0;
    }
    long int leftSize = validateHelper(node->left);
    long int rightSize = validateHelper(node->right);
    int leftHeight = (node->left == nullptr)? -1 : node->left->height;
    int rightHeight = (node->right == nullptr)? -1 : node->right->height;
    if(node->left != nullptr)
    {
        assert(node->left->parent == node);
        assert(node->left != node);
        assert(node->left != node->right);
        assert(node->parent != node->left->parent); // cycle!!
    }
    if(node->right != nullptr)
    {
        assert(node->right->parent == node);
        assert(node->right != node);
        assert(node->left != node->right);
        assert(node->parent != node->right->parent); // cycle!!
    }
    if(node->height != std::max(leftHeight, rightHeight) +1) std::cout << node->height << " left=" << leftHeight << ", right=" << rightHeight << ", expected:" << std::max(leftHeight, rightHeight) +1 << std::endl;
    assert(node->height == std::max(leftHeight, rightHeight) + 1);
    assert(node->leftSize == leftSize);
    assert(node->rightSize == rightSize);
    return leftSize + rightSize + node->duplications;
}
#endif 


