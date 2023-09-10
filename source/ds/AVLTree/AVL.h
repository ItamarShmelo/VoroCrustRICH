#ifndef _AVL_H
#define _AVL_H

#include "../BinaryTree.h"

template<typename T>
class AVLTree : public BinaryTree<T>
{
public:
    inline AVLTree(typename BinaryTree<T>::Node *root, const typename BinaryTree<T>::Compare &compare): BinaryTree<T>(root, compare){};
    inline AVLTree(const typename BinaryTree<T>::Compare &compare): BinaryTree<T>(compare){};
    inline AVLTree(): BinaryTree<T>(){};
    template<typename InputIterator>
    inline AVLTree(const InputIterator &begin, const InputIterator &end): AVLTree<T>(){for(InputIterator it = begin; it < end; it++) this->insert(*it);};
    virtual ~AVLTree() = default;

    virtual bool insert(const T &value) override;
    virtual bool remove(const T &value) override;

protected:
    inline int getBalance(const typename BinaryTree<T>::Node *node) const
    {if(node == nullptr) return 0; int leftHeight = (node->left == nullptr)? -1 : node->left->height; int rightHeight = (node->right == nullptr)? -1 : node->right->height; return leftHeight - rightHeight;}
    void recursiveUpdateNodeInfo(typename BinaryTree<T>::Node *node) override;
    long int validateHelper(const typename AVLTree<T>::Node *node) const override;
};

#endif // _AVL_H