#include "output.h"

void printArray(int rank, std::vector<int> &data)
{
    std::cout << "Final Process " << rank << ": " << data.size() << " Elements, as follows: " << std::endl;
    std::for_each(data.begin(), data.end(), [](int &value){std::cout << value << " ";});
    std::cout << std::endl;
}

void printAllArrays(std::vector<int> &myData)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int wait = 0;
    if(rank != 0)
    {
        MPI_Recv(&wait, 1, MPI_INT, rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    printArray(rank, myData);
    if(rank != size - 1)
    {
        MPI_Send(&wait, 1, MPI_INT, rank+1, 3, MPI_COMM_WORLD);
    }
}

inline void printArraySize(int rank, std::vector<int> &data)
{
    std::cout << "Final Process " << rank << ": " << data.size() << " Elements." << std::endl;
}

void printAllArraysSizes(std::vector<int> &myData)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int wait = 0;
    if(rank != 0)
    {
        MPI_Recv(&wait, 1, MPI_INT, rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    printArraySize(rank, myData);
    if(rank != size - 1)
    {
        MPI_Send(&wait, 1, MPI_INT, rank+1, 3, MPI_COMM_WORLD);
    }
}


void printToFile(std::string fileName, std::vector<int> &vector)
{
	std::ofstream file;
    file.open(fileName, std::ios::out | std::ios::app);
	if(!file.good())
    {
		std::cerr << "Error in writing to " << fileName << std::endl;
		exit(EXIT_FAILURE);
	}
	else
    {
		for(int &value : vector)
        {
            file << value << std::endl;
        }
	}
	file.close();
}