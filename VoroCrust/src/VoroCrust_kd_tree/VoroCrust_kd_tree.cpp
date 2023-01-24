#include "VoroCrust_kd_tree.hpp"
#include <functional>
#include <numeric>



std::shared_ptr<Node> newNode(int index, int axis){
    std::shared_ptr<Node> node = std::make_shared<Node>();

    node->index = index;
    node->axis = axis;

    return node;
}
