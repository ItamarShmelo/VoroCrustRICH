#include <iostream> // todo remove
#include <vector>
#include <algorithm>
#include "mpi.h"

#define MAX(x, y) (((x) > (y))? (x) : (y))
#define ROOT 0
#define LENGTH_TAG 0
#define DATA_TAG 1
#define SAMPLE_DEPTH size

/**
 * Using the sort function, each one of the processes sends its data, and an efficient parallel sorting is being performed.
 * @param data the data of the process, as a vector
*/
void sort(std::vector<int> &data);