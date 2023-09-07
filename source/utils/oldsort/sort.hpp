/**
 * DEPRECATED
 * \
 * \author Maor Mizrachi
 * \brief Implements parallel sort.
*/
#ifndef _RICH_PARALLEL_SORT
#define _RICH_PARALLEL_SORT

#include <iostream> // todo remove
#include <vector>
#include <algorithm>
#include <mpi.h>

#define SAMPLE_DEPTH size
#define DATA_TAG 305
#define SAMPLES_TAG 306
#define MAX_PER_MESSAGE 150
#define ROOT 0

namespace
{
    /**
     * Gets a data vector, a number 'n', and an additional container, and samples 'n' numbers from the vector,
     * roughly equally distant, appending them to the additional container.
     * @param data the data to sample from
     * @param n the number of elements to sample
     * @return the vector of the sampled elements
    */
    template<typename T>
    std::vector<T> _sample(std::vector<T> &data, size_t n)
    {
        std::vector<T> sampled;
        if(data.size() < n)
        {
            // not enough data to sample 'n' values, just append all
            sampled.insert(sampled.end(), data.cbegin(), data.cend());
            return sampled;
        }
        
        size_t gap = data.size() / n; // gap between samples
        size_t i = data.size() - 1; // the index chosen to be sampled

        while(i >= n-1 and n > 0)
        {
            sampled.push_back(data[i]);
            i -= gap;
            n--;
        }

        return std::vector<T>(sampled.rbegin(), sampled.rend());
    }

    /**
     * Expects to receive elements from each one of the processes, and appending them to 'data' (in an arbitary order,
     * according to which are the first processes messages to arrive)
     * @param data the container to hold the data
    */
    template<typename T>
    static void _receive(std::vector<T> &data)
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        std::vector<MPI_Request> requests;

        int received = 0;

        while(received != size - 1)
        {
            MPI_Status status;
            int length;
            MPI_Probe(MPI_ANY_SOURCE, SAMPLES_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &length);
            if(length > 0)
            {
                data.resize(data.size() + length / sizeof(T));
                /*
                TODO: maybe parallelizing the receive, making them non-blocking, might help in optimization.
                */
                requests.push_back(MPI_REQUEST_NULL);
                MPI_Irecv(&(*(data.end() - length / sizeof(T))), length * sizeof(T), MPI_BYTE, status.MPI_SOURCE, SAMPLES_TAG, MPI_COMM_WORLD, &requests[requests.size() - 1]);
            }
            received++;
        }
        MPI_Waitall(requests.size() ,&requests[0], MPI_STATUSES_IGNORE);
    }

    /**
     * Using this function, each one of the processes sends its samples, and updates its ranges list,
     * the list that determines which numbers will be sent to which processors.
     * @param mySamples the processor's samples vector
    */
    template<typename T>
    std::vector<T> _getRanges(std::vector<T> &mySamples)
    {
        int size, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        std::vector<T> sampleOfSamples(size);

        if(rank != ROOT)
        {
            MPI_Send(&(*(mySamples.cbegin())), mySamples.size() * sizeof(T), MPI_BYTE, ROOT, SAMPLES_TAG, MPI_COMM_WORLD);
        }
        else
        {
            // the root first appends all the samples to its own samples vector, then samples again
            _receive(mySamples);
            std::sort(mySamples.begin(), mySamples.end());
            sampleOfSamples = _sample(mySamples, size);
            std::sort(sampleOfSamples.begin(), sampleOfSamples.end());
        }

        // the root sends the ranges array to everyone using broadcast
        MPI_Bcast(&(*(sampleOfSamples.begin())), size * sizeof(T), MPI_BYTE, ROOT, MPI_COMM_WORLD);
        return sampleOfSamples;
    }

    /**
     * In the end of this function, the relevant data is sent to the relevant processor, according to the ranges array.
     * @param data the data of the current processor
     * @param ranges the ranges array
    */
    template<typename T>
    void _rearrangeByRanges(std::vector<T> &data, std::vector<T> &ranges)
    {
        int size, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        std::vector<MPI_Request> requests;
        int currentProcessor = 0;
        size_t index = 0;

        std::vector<T> myValues;

        // send stage:
        while(index < data.size())
        {
            std::vector<T> valuesToProcessor;
            while(index < data.size() and data[index] <= ranges[currentProcessor])
            {
                // get to the furthest value (`data[index]`) matching to this processor
                valuesToProcessor.push_back(data[index]);
                if(valuesToProcessor.size() == MAX_PER_MESSAGE)
                {
                    if(currentProcessor == rank)
                    {
                        myValues.insert(myValues.end(), valuesToProcessor.begin(), valuesToProcessor.end());
                    }
                    else
                    {
                        requests.push_back(MPI_REQUEST_NULL);
                        MPI_Send(&valuesToProcessor[0], sizeof(T) * valuesToProcessor.size(), MPI_BYTE, currentProcessor, DATA_TAG, MPI_COMM_WORLD);
                    }
                    valuesToProcessor.clear();
                }
                index++;

            }
            if(valuesToProcessor.size() > 0)
            {
                if(currentProcessor == rank)
                {
                    myValues.insert(myValues.end(), valuesToProcessor.begin(), valuesToProcessor.end());
                }
                else
                {
                    requests.push_back(MPI_REQUEST_NULL);
                    MPI_Send(&valuesToProcessor[0], sizeof(T) * valuesToProcessor.size(), MPI_BYTE, currentProcessor, DATA_TAG, MPI_COMM_WORLD);
                }
            }
            // now, if we haven't finished, get to the processor that matches the current data
            if(index != data.size())
            {
                while(!(data[index] <= ranges[currentProcessor]))
                {
                    currentProcessor += 1;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        data.clear();
        data = std::move(myValues);

        // receive stage:
        int newMessage = 0;
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE, DATA_TAG, MPI_COMM_WORLD, &newMessage, &status);
        while(newMessage != 0)
        {
            int length;
            MPI_Get_count(&status, MPI_BYTE, &length);
            size_t sizeBefore = data.size();
            data.resize(sizeBefore + (length / sizeof(T)));
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(&data[sizeBefore], length, MPI_BYTE, status.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD, &requests[requests.size() - 1]);
            MPI_Iprobe(MPI_ANY_SOURCE, DATA_TAG, MPI_COMM_WORLD, &newMessage, &status);
        }

        MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE); // make sure everyone received the message
    }
};

template<typename T>
void parallelSort(std::vector<T> &data)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // phase 1: sort my own data
    std::sort(data.begin(), data.end());
    
    // phase 2: send the root 'size' elements, equally distant
    std::vector<T> samples = _sample(data, SAMPLE_DEPTH * size);
    
    // phase 3: the root determines the ranges each one of the processes should hold
    std::vector<T> ranges = _getRanges(samples);

    // phase 4: after receiving the ranges, we now rearrange the data between processes
    _rearrangeByRanges(data, ranges);
    
    // phase 5: re-sort my new array
    std::sort(data.begin(), data.end());
}

#endif // _RICH_PARALLEL_SORT
//#endif