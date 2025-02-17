#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include "trees/trees.hpp"
#include "RMPS/CornersRMPS.hpp"
#include "RMPS/EdgesRMPS.hpp"
#include "RMPS/FacesRMPS.hpp"
#include "RMPS/SliverDriver.hpp"
#include <sstream>
#include <filesystem>
#include <memory>

std::pair<std::vector<Seed>, std::vector<Seed>> randomSampleVolumeSeeds(PL_ComplexPtr const& plc, std::pair<std::vector<Seed>, std::vector<Seed>> const& zones_boundary_seeds, double const maxSize, Trees const& trees, double const L_Lipschitz);

bool enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree, double const L_Lipschitz);

VoroCrust_KD_Tree_Ball makeSeedBallTree(std::vector<Seed> const& seeds);

std::pair<std::vector<Seed>, std::vector<Seed>> determineSeedsInOut(std::vector<Seed> const& seeds, PL_ComplexPtr const& plc);

std::vector<Seed> getSeedsFromBallTree(VoroCrust_KD_Tree_Ball const& ball_tree);

void dumpSeeds(std::string const& dirname, std::vector<Seed> const& seeds);

std::vector<Seed> load_dumpSeeds(std::filesystem::path const& dirname);

#endif /* VOROCRUST_ALGORITHM */