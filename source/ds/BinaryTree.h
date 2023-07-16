/**
 * Implements a d-dimensional data structure for fast range queries, called range trees.
 * See more here: https://www.cs.umd.edu/class/fall2020/cmsc420-0201/Lects/lect17-range-tree.pdf
 * 
 * @author maor miz
*/

#ifndef _RICH_BST_H
#define _RICH_BST_H

#include <functional>
#include <vector>
#include <utility>
#include <assert.h>

// #define DEBUG_MODE
#undef DEBUG_MODE

template<typename T>
class BinaryTree
{
protected:
    using Compare=std::function<bool(const T&, const T&)>; 

    class Node
    {
    public:
        T value;
        Node *left, *right;
        Node *parent;
        long int leftSize;
        long int rightSize;
        int height;
        int duplications;
        T maxInSubtree;
        T minInSubtree;

        Node(const T &value): value(value), left(nullptr), right(nullptr), parent(nullptr), leftSize(0), rightSize(0), height(0), duplications(1), maxInSubtree(value), minInSubtree(value){};
        Node(): Node(T()){};
        virtual ~Node() = default;
    };

    const Node *findNode(const T &value) const;
    const Node *findParent(const T &value) const;
    const Node *findClosestAbove(const T &value) const;
    const Node *findClosestBelow(const T &value) const;
    const Node *findSuccessor(const T &value) const;
    const Node *findPredecessor(const T &value) const;
    const Node *findMax() const{const Node *node = this->getRoot(); while(node != nullptr and node->right != nullptr) node = node->right; return node;};
    const Node *findMin() const{const Node *node = this->getRoot(); while(node != nullptr and node->left != nullptr) node = node->left; return node;};
    inline Node *findNode(const T &value){return const_cast<Node*>((std::as_const(*this)).findNode(value));};
    inline Node *findParent(const T &value){return const_cast<Node*>((std::as_const(*this)).findParent(value));};
    inline Node *findClosestAbove(const T &value){return const_cast<Node*>((std::as_const(*this)).findClosestAbove(value));};
    inline Node *findClosestBelow(const T &value){return const_cast<Node*>((std::as_const(*this)).findClosestBelow(value));};
    inline Node *findSuccessor(const T &value){return const_cast<Node*>((std::as_const(*this)).findSuccessor(value));};
    inline Node *findPredecessor(const T &value){return const_cast<Node*>((std::as_const(*this)).findPredecessor(value));};
    inline Node *findMax(){return const_cast<Node*>((std::as_const(*this)).findMax());};
    inline Node *findMin(){return const_cast<Node*>((std::as_const(*this)).findMin());};
    void fastRotateRight(Node *node);
    void fastRotateLeft(Node *node);
    inline void rotateRight(Node *node){this->fastRotateRight(node); this->updateNodeInfo(node);};
    inline void rotateLeft(Node *node){this->fastRotateLeft(node); this->updateNodeInfo(node);};
    Node *_insert(const T &value);
    Node *_remove(Node *node);

    inline virtual Node *getRoot(){return this->root;};
    inline virtual const Node *getRoot() const{return this->root;};
    inline virtual void setRoot(Node *other){this->root = other;};

    virtual inline std::vector<T> getAllDecendants(Node *node){std::vector<T> vec; this->getAllDecendantsHelper(node, vec); return vec;};

    virtual void updateNodeInfo(Node *node);
    inline virtual void recursiveUpdateNodeInfo(Node *node){while(node != nullptr){this->updateNodeInfo(node); node = node->parent;}};

    #ifdef DEBUG_MODE
    void printHelper(const Node *node, int tabs) const;
    virtual long int validateHelper(const Node *node) const;
    #endif 

    
    Node *root;
    Compare compare;
    long int treeSize;


public:
    inline BinaryTree(BinaryTree<T>::Node *root, const Compare &compare): root(root), compare(compare), treeSize(0){};
    inline BinaryTree(const Compare &compare): BinaryTree(nullptr, compare){};
    inline BinaryTree(): BinaryTree([](const T &a, const T &b){return a <= b;}){};
    template<typename InputIterator>
    inline BinaryTree(const InputIterator &begin, const InputIterator &end): BinaryTree(){for(InputIterator it = begin; it < end; it++) this->insert(*it);}
    inline virtual ~BinaryTree(){this->deleteSubtree(this->getRoot());};

    inline bool trySuccessor(const T &value, T &successor) const{const Node *successorNode = this->findSuccessor(value); if(successorNode != nullptr){successor = successorNode->value; return true;}else{successor = T(); return false;}};
    inline const T &successor(const T &value) const{return this->findSuccessor(value)->value;};
    inline bool tryCeil(const T &value, T &ceil) const{const Node *ceilNode = this->findClosestAbove(value); if(ceilNode != nullptr){ceil = ceilNode->value; return true;}else{ceil = T(); return false;}};
    inline const T &ceil(const T &value) const{return this->findClosestAbove(value)->value;};
    inline bool tryFloor(const T &value, T &floor) const{const Node *floorNode = this->findClosestBelow(value); if(floorNode != nullptr){floor = floorNode->value; return true;}else{floor = T(); return false;}};
    inline const T &floor(const T &value) const{return this->findClosestBelow(value)->value;};
    inline const T &getMax() const{assert(this->getRoot() != nullptr); return this->getRoot()->maxInSubtree;};
    inline const T &getMin() const{assert(this->getRoot() != nullptr); return this->getRoot()->minInSubtree;};
    inline const T &getRootValue() const{assert(this->getRoot() != nullptr); return this->getRoot()->value;}
    inline bool isEmpty() const{return this->getRoot() == nullptr;};
    inline bool find(const T &value) const{return this->findNode(value) != nullptr;};
    virtual bool insert(const T &value);
    virtual bool remove(const T &value);
    inline long int size() const{return this->treeSize;};

    void deleteSubtree(Node *node);

    #ifdef DEBUG_MODE
    inline void rotateRight(const T &value){this->rotateRight(this->findNode(value));};
    inline void rotateLeft(const T &value){this->rotateLeft(this->findNode(value));};
    #endif // _DEBUG_MODE

    #ifdef DEBUG_MODE
    inline void print() const{this->printHelper(this->getRoot(), 0);};
    inline bool validate() const{return this->treeSize == this->validateHelper(this->getRoot());};
    #endif // DEBUG_MODE

    struct Iterator
    {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        Iterator(Node *node, BinaryTree<T> *tree): node(node), tree(tree){};
        Iterator(const Iterator &other) = default;
        inline Iterator &operator=(const Iterator &other){this->tree = other.tree; this->node = other.node; return (*this);};
        inline const T &operator*() const{return node->value;};
        inline virtual Iterator &operator++(){if(this->node != nullptr) this->node = this->tree->findSuccessor(this->node->value); return (*this);};
        inline virtual Iterator &operator--(){if(this->node != nullptr) this->node = this->tree->findPredecessor(this->node->value); return (*this);};
        inline bool operator==(const Iterator &other) const{return this->node == other.node and this->tree == other.tree;};
        inline bool operator!=(const Iterator &other) const{return !this->operator==(other);};

        Node *node;
        BinaryTree<T> *tree;
    };

    struct ReverseIterator : Iterator
    {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        ReverseIterator(Node *node, BinaryTree<T> *tree): Iterator(node, tree){};
        ReverseIterator(const ReverseIterator &other) = default;

        inline Iterator &operator++() override{return this->Iterator::operator--();};
        inline Iterator &operator--() override{return this->Iterator::operator++();};
    };

    friend struct Iterator;

    typedef Iterator iterator;
    typedef const Iterator const_iterator;
    typedef ReverseIterator reverse_iterator;
    typedef const ReverseIterator const_reverse_iterator;

    iterator begin(){return iterator(this->findMin(), this);};
    const_iterator cbegin(){return const_iterator(this->findMin(), this);};
    reverse_iterator rbegin(){return reverse_iterator(this->findMax(), this);};
    const_reverse_iterator crbegin(){return const_reverse_iterator(this->findMax(), this);};
    iterator end(){return iterator(nullptr, this);};
    const_iterator cend(){return const_iterator(nullptr, this);};
    reverse_iterator rend(){return reverse_iterator(nullptr, this);};
    const_reverse_iterator crend(){return const_reverse_iterator(nullptr, this);};

private:
    void getAllDecendantsHelper(Node *node, std::vector<T> &vec);
};

#endif // _RICH_BST_H