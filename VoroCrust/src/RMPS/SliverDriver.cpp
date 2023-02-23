#include "SliverDriver.hpp"
#include <iostream>
#include "../../../source/misc/utils.hpp"

SliverDriver::SliverDriver(double const L_Lipschitz_) : L_Lipschitz(L_Lipschitz_), r_new_corner_balls(), r_new_edge_balls(), r_new_face_balls(), number_of_slivers_eliminated(0) {}
