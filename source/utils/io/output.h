#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "mpi.h"

void printArray(int rank, std::vector<int> &data);

void printArraySize(int rank, std::vector<int> &data);

void printAllArrays(std::vector<int> &myData);

void printAllArraysSizes(std::vector<int> &myData);

void printToFile(std::string fileName, std::vector<int> &vector);
