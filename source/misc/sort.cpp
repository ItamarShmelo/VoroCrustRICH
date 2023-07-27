#if 0
/*
Recommended reading about the samplesort algorithm here:
https://brunomaga.github.io/Distributed-Sort
*/

#include "sort.hpp"

/**
 * Gets a data vector, a number 'n', and an additional container, and samples 'n' numbers from the vector,
 * roughly equally distant, appending them to the additional container.
 * @param data the data to sample from
 * @param n the number of elements to sample
 * @param sampled the container to append the sampled elements into
*/
void _sample(std::vector<int> &data, size_t n, std::vector<int> &sampled)
{
    if(data.size() < n)
    {
        // not enough data to sample 'n' values, just append all
        sampled.insert(sampled.end(), data.cbegin(), data.cend());
        return;
    }
    
    size_t gap = data.size() / n; // gap between samples
    size_t i = data.size() - 1; // the index chosen to be sampled

    while(i >= n-1 and n > 0)
    {
        sampled.push_back(data[i]);
        i -= gap;
        n--;
    }
}

/**
 * gets a data container, and an array of vectors, one for each process, containing elements to send to the process.
 * In the end of this function, for each one of the processes,
 * the given container contains all the elements intended to be sent to it.
 * @param data the data container - will be overridden
 * @param toRanks an array containing a vector for each rank, in which the elements that should be sent to the rank are stored
*/
static void _distribute(std::vector<int> &data, std::vector<int> *&toRanks)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int i = 0; i < size; i++)
    {
        if(i == rank)
        {
            continue;
        }
        int toSend = toRanks[i].size();
        MPI_Send(&toSend, 1, MPI_INT, i, LENGTH_TAG, MPI_COMM_WORLD);
        if(toSend > 0)
        {
            MPI_Send(&(*(toRanks[i].cbegin())), toSend, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD);
        }
    }
    data.clear();
    data.insert(data.end(), toRanks[rank].cbegin(), toRanks[rank].cend());
}

/**
 * Gets a value and the ranges array, and determines to which processor the value should be sent.
 * @param value the value
 * @param ranges the ranges array
*/
static int _getRankForValue(int &value, int *&ranges)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(value <= ranges[0])
    {
        return 0; // the value matches to the first processor
    }

    // apply a binary search
    int begin = 0, end = size;
    while(true)
    {
        int mid = (begin + end) / 2;
        if(ranges[mid-1] < value and value <= ranges[mid])
        {
            return mid; // found the correct processor
        }
        if(ranges[mid-1] >= value)
        {
            end = mid; // go left, the value is smaller
        } 
        else if(ranges[mid] < value)
        {
            begin = mid; // go right, the value is bigger
        }
    }
}

/**
 * Expects to receive elements from each one of the processes, and appending them to 'data' (in an arbitary order,
 * according to which are the first processes messages to arrive)
 * @param data the container to hold the data
 * @param size the communicator size
*/
static void _receive(std::vector<int> &data)
{
    int size, rank; // todo remove rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int received = 0;
    MPI_Status status;

    while(received != size-1)
    {
        int length;
        MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, LENGTH_TAG, MPI_COMM_WORLD, &status);
        if(length > 0)
        {
            data.resize(data.size() + length);
            MPI_Recv(&(*(data.end() - length)), length, MPI_INT, status.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        received++;
    }
}

/**
 * In the end of this function, the relevant data is sent to the relevant processor, according to the ranges array.
 * @param data the data of the current processor
 * @param ranges the ranges array
*/
static void _rearrangeByRanges(std::vector<int> &data, int *&ranges)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> *toRanks = new std::vector<int>[size];
    for(int value : data)
    {
        toRanks[_getRankForValue(value, ranges)].push_back(value);
    }
    _distribute(data, toRanks); // send data to others
    delete[] toRanks;
    _receive(data); // receive my own data
}


/**
 * Using this function, each one of the processes sends its samples, and updates its ranges list,
 * the list that determines which numbers will be sent to which processors.
 * @param mySamples the processor's samples vector
 * @param ranges the ranges array, to be updated
*/
static void _getRanges(std::vector<int> &mySamples, int *&ranges)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank != ROOT)
    {
        // send the root my samples
        int length = mySamples.size();
        
        MPI_Send(&length, 1, MPI_INT, ROOT, LENGTH_TAG, MPI_COMM_WORLD);
        if(length > 0)
        {
            MPI_Send(&(*(mySamples.cbegin())), mySamples.size(), MPI_INT, ROOT, DATA_TAG, MPI_COMM_WORLD);
        }
    }
    else
    {
        // the root first appends all the samples to its own samples vector, then samples again
        _receive(mySamples);
        std::sort(mySamples.begin(), mySamples.end());
        std::vector<int> sampleOfSamples;
        _sample(mySamples, size, sampleOfSamples);
        std::sort(sampleOfSamples.begin(), sampleOfSamples.end());
        
        // now update the ranges array
        for(int i = 0; i < size; i++)
        {
            ranges[i] = sampleOfSamples[i];
        }
    }

    // the root sends the ranges array using broadcast
    MPI_Bcast(ranges, size, MPI_INT, ROOT, MPI_COMM_WORLD);
}

/**
 * Using the sort function, each one of the processes sends its data, and an efficient parallel sorting is being performed.
 * @param data the data of the process, as a vector
*/
void sort(std::vector<int> &data)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // phase 1: sort my own data
    std::sort(data.begin(), data.end());
    
    // phase 2: send the root 'size' elements, equally distant
    std::vector<int> samples;

    _sample(data, SAMPLE_DEPTH * size, samples);
    
    // phase 3: the root determines the ranges each one of the processes should hold
    int *ranges = new int[size];
    _getRanges(samples, ranges);

    // phase 4: after receiving the ranges, we now rearrange the data between processes
    _rearrangeByRanges(data, ranges);
    delete[] ranges;

    // phase 5: re-sort my new array
    std::sort(data.begin(), data.end());
}
#endif