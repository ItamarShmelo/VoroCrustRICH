#ifndef _DISTRIBUTED_OCTTREE_HPP
#define _DISTRIBUTED_OCTTREE_HPP

#ifdef RICH_MPI

#include <iostream> // todo remove
#include <vector>
#include <assert.h>
#include <utility>
#include <array>
#include <bitset>
#include <boost/container/flat_set.hpp>
#include <mpi.h>
#include "ds/OctTree/OctTree.hpp"

#define UNDEFINED_OWNER -1

template<typename T>
class DistributedOctTree
{

public:
    DistributedOctTree(const MPI_Comm &comm, const OctTree<T> *tree);
    inline DistributedOctTree(const OctTree<T> *tree): DistributedOctTree(MPI_COMM_WORLD, tree){};
    ~DistributedOctTree(){delete this->octTree;};

    void print() const{this->octTree->print();};
    boost::container::flat_set<int> getIntersectingRanks(const _Sphere<T> &sphere) const;
    inline boost::container::flat_set<int> getIntersectingRanks(const T &center, const typename T::coord_type radius) const{return this->getIntersectingRanks(_Sphere(center, radius));};

    #ifdef DEBUG_MODE
    inline bool validate() const{if(this->octTree != nullptr) return this->validateHelper(this->octTree->getRoot()); return true;};
    #endif // DEBUG_MODE

private:
    class _Wrapper
    {
        friend class DistributedOctTree;
    public:
        using coord_type = typename T::coord_type;

        T value;
        int owner; // rank of the owner

        _Wrapper(const T &value, int owner): value(value), owner(owner){};
        _Wrapper(): _Wrapper(T(), UNDEFINED_OWNER){};
        _Wrapper(const _Wrapper &other): _Wrapper(other.value, other.owner){};
        _Wrapper &operator=(const _Wrapper &other){this->value = other.value; this->owner = other.owner; return (*this);};
        inline typename T::coord_type operator[](size_t idx) const{return this->value[idx];};
        inline typename T::coord_type &operator[](size_t idx){return this->value[idx];};
        inline _Wrapper operator+(const _Wrapper &other) const{int owner = (this->owner == other.owner)? this->owner : UNDEFINED_OWNER; return _Wrapper(this->value + other.value, owner);};
        inline _Wrapper operator-(const _Wrapper &other) const{int owner = (this->owner == other.owner)? this->owner : UNDEFINED_OWNER; return _Wrapper(this->value - other.value, owner);};
        inline _Wrapper operator*(double constant) const{return _Wrapper(this->value * constant, owner);};
        inline _Wrapper operator/(double constant) const{return this->operator*(1/constant);};

        friend std::ostream &operator<<(std::ostream &stream, const _Wrapper &wrapper)
        {
            return stream << wrapper.value << " [owner: " << wrapper.owner << "]";
        }

    };

    OctTree<_Wrapper> *octTree = nullptr;
    MPI_Comm comm;
    int rank, size;

    void buildTreeHelper(typename OctTree<_Wrapper>::OctTreeNode *newNode, const typename OctTree<T>::OctTreeNode *node);
    void buildTree(const OctTree<T> *tree);

    #ifdef DEBUG_MODE
    bool validateHelper(const typename OctTree<_Wrapper>::OctTreeNode *node) const;
    #endif // DEBUG_MODE
};

template<typename T>
DistributedOctTree<T>::DistributedOctTree(const MPI_Comm &comm, const OctTree<T> *tree): comm(comm)
{
    MPI_Comm_rank(this->comm, &this->rank);
    MPI_Comm_size(this->comm, &this->size);
    this->buildTree(tree);
}

template<typename T>
void DistributedOctTree<T>::buildTreeHelper(typename OctTree<_Wrapper>::OctTreeNode *newNode, const typename OctTree<T>::OctTreeNode *node)
{
    assert(newNode != nullptr);

    // MPI_Barrier(this->comm);
        
    unsigned char valueToSend = 0; // assumes `CHILDREN` is 8. this variable contains 1 in the `i`th bit iff child `i` exists
    if(node != nullptr)
    {
        for(int i = 0; i < CHILDREN; i++)
        {
            bool bit = (node->children[i] != nullptr || (node->isValue and newNode->getChildNumberContaining(_Wrapper(node->value, UNDEFINED_OWNER)) == i));
            valueToSend |= (bit << i);
        }
    }
    std::vector<char> childBuff(this->size);
    MPI_Allgather(&valueToSend, 1, MPI_BYTE, &childBuff[0], 1, MPI_BYTE, this->comm);

    bool somebodyHolds = false;
    for(int i = 0; i < CHILDREN; i++)
    {
        bool recursiveBuild = false;
        int containingValue = UNDEFINED_OWNER; // who has the `i`th child
        for(int _rank = 0; _rank < this->size; _rank++)
        {
            if(((childBuff[_rank] >> i) & 0x1) == 1)
            {
                if(containingValue == UNDEFINED_OWNER)
                {
                    somebodyHolds = true;
                    containingValue = _rank;
                }
                else
                {
                    recursiveBuild = true;
                    break;
                }
            }
        }
        // if somebody holds the value
        if(containingValue != UNDEFINED_OWNER)
        {
            // someone holds the `i`th child
            newNode->createChild(i); // creates the child in my node
            if(recursiveBuild)
            {
                // there are several holders, call recursive build
                const typename OctTree<T>::OctTreeNode *nextNode = nullptr;
                if(node != nullptr)
                {
                    if(node->isValue)
                    {
                        nextNode = newNode->children[i]->boundingBox.contains(_Wrapper(node->value, UNDEFINED_OWNER))? node : nullptr;
                    }
                    else
                    {
                        nextNode = node->children[i];
                    }
                }
                else
                {
                    nextNode = nullptr;
                }
                this->buildTreeHelper(newNode->children[i], nextNode);
                newNode->children[i]->value.owner = UNDEFINED_OWNER;
                newNode->children[i]->boundingBox.ll.owner = newNode->children[i]->boundingBox.ur.owner = UNDEFINED_OWNER;
            }
            else
            {
                // there is only one holder, set its owner
                newNode->children[i]->value.owner = containingValue;
                newNode->children[i]->boundingBox.ll.owner = newNode->children[i]->boundingBox.ur.owner = containingValue;
                newNode->children[i]->isValue = true;
            }
        }
        else
        {
            newNode->children[i] = nullptr;
        }
    }
}

#ifdef DEBUG_MODE
template<typename T>
bool DistributedOctTree<T>::validateHelper(const typename OctTree<_Wrapper>::OctTreeNode *node) const
{
    if(node == nullptr)
    {
        return true;
    }
    bool hasChildren = false;
    for(int i = 0; i < CHILDREN; i++)
    {
        if(node->children[i] != nullptr)
        {
            hasChildren = true;
            break;
        }
    }
    if(hasChildren)
    {
        assert(node->isValue);
        assert(node->boundingBox.ll.owner != UNDEFINED_OWNER);
        assert(node->boundingBox.ur.owner != UNDEFINED_OWNER);
        assert(node->boundingBox.ll.owner == node->boundingBox.ur.owner);
    }
    for(int i = 0; i < CHILDREN; i++)
    {
        assert(this->validateHelper(node->children[i]));
    }
    return true;
}
#endif // DEBUG_MODE

template<typename T>
void DistributedOctTree<T>::buildTree(const OctTree<T> *tree)
{
    assert(this->octTree == nullptr);
    if(tree == nullptr or tree->getRoot() == nullptr)
    {
        return;
    }
    this->octTree = new OctTree<_Wrapper>(_Wrapper(tree->getRoot()->boundingBox.ll, UNDEFINED_OWNER), _Wrapper(tree->getRoot()->boundingBox.ur, UNDEFINED_OWNER));
    this->buildTreeHelper(this->octTree->getRoot(), tree->getRoot());
}

template<typename T>
boost::container::flat_set<int> DistributedOctTree<T>::getIntersectingRanks(const _Sphere<T> &sphere) const
{
    boost::container::flat_set<int> ranks;
    for(const _Wrapper &point : this->octTree->range(_Sphere<_Wrapper>(_Wrapper(sphere.center, UNDEFINED_OWNER), sphere.radius)))
    {
        ranks.insert(point.owner);
    }
    return ranks;
}

#endif // RICH_MPI

#endif // _DISTRIBUTED_OCTTREE_HPP

