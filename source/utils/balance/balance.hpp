#include <iostream> // todo remove
#include <vector>
#include <boost/heap/priority_queue.hpp>
#include <mpi.h>
#include <functional>
#include <utility>
#include <assert.h>

// #define DEBUG_MODE
#undef DEBUG_MODE

#define BALANCE_ROOT 0
#define BALANCE_MINIMAL_REQUEST_TAG 901
#define BALANCE_MINIMAL_RESPONSE_TAG 902
#define BALANCE_STOP_ITERATION_TAG 903
#define BALANCE_FINISHED_TAG 904
#define BALANCE_NEW_BOUND_TAG 907

template<typename T>
class BalanceJob
{
private:
    using Compare = std::function<bool(const T&, const T&)>; 

    struct _T_Wrapper
    {
        T value;
        int rank;

        _T_Wrapper(const T &value, int rank): value(value), rank(rank){};
        /*
        bool operator==(const _T_Wrapper &other) const{return this->value == other.value;};
        bool operator<(const _T_Wrapper &other) const{return this->value < other.value;};
        bool operator<=(const _T_Wrapper &other) const{return this->value <= other.value;};
        bool operator>(const _T_Wrapper &other) const{return this->value > other.value;};
        bool operator>=(const _T_Wrapper &other) const{return this->value >= other.value;};
        */
    };

    template<typename _Compare>
    struct _T_Wrapper_Comp
    {
        _Compare compare;
        inline _T_Wrapper_Comp(const Compare &compare): compare(compare){};
        inline bool operator()(const _T_Wrapper &lhs, const _T_Wrapper &rhs) const{return !this->compare(lhs.value, rhs.value);};
    };

    using Heap = boost::heap::priority_queue<_T_Wrapper, boost::heap::compare<_T_Wrapper_Comp<Compare>>>;

    void helpBuildHeap();
    Heap buildHeap();
    void heapReplaceMin(Heap &heap);
    bool answerToRequests();
    void announceBound(const T &bound);
    void kthOrderStatistics(std::vector<size_t> &orders, std::vector<T> &bounds);
    void getNewBounds(std::vector<T> &bounds);

    std::vector<T> values;
    std::vector<T> newValues;
    std::vector<MPI_Request> requests;
    // std::vector<std::vector<T>> buffers;
    Compare comparator;
    int rank, size;
    size_t heapPos;

public:
    inline BalanceJob(const std::vector<T> &values, Compare comparator): values(values), comparator(comparator), heapPos(0)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
        MPI_Comm_size(MPI_COMM_WORLD, &this->size);
        std::sort(this->values.begin(), this->values.end(), this->comparator);
    }
    inline explicit BalanceJob(const std::vector<T> &values): BalanceJob(values, [](const T &a, const T &b){return a < b;}){};
    std::vector<T> getBorders();
};

template<typename T>
void BalanceJob<T>::helpBuildHeap()
{
    if(this->rank == BALANCE_ROOT)
    {
        return;
    }
    if(this->heapPos >= this->values.size())
    {
        int dummy = 0;
        MPI_Send(&dummy, 1, MPI_INT, BALANCE_ROOT, BALANCE_STOP_ITERATION_TAG, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&this->values[this->heapPos++], sizeof(T), MPI_BYTE, BALANCE_ROOT, BALANCE_MINIMAL_RESPONSE_TAG, MPI_COMM_WORLD);
    }
}

template<typename T>
typename BalanceJob<T>::Heap BalanceJob<T>::buildHeap()
{
    Heap heap(_T_Wrapper_Comp<Compare>(this->comparator));
    if(this->values.size() > this->heapPos)
    {
        heap.push(_T_Wrapper(this->values[this->heapPos++], this->rank));
    }
    for(int _rank = 0; _rank < this->size; _rank++)
    {
        if(_rank == this->rank)
        {
            continue;
        }
        MPI_Status status;
        MPI_Probe(_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch(status.MPI_TAG)
        {
            case BALANCE_STOP_ITERATION_TAG:
            {   
                int dummy;
                MPI_Recv(&dummy, 1, MPI_INT, _rank, BALANCE_STOP_ITERATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::cerr << "reached here, should not" << std::endl;
                break;
            }
            case BALANCE_MINIMAL_RESPONSE_TAG:
            {
                T value;
                MPI_Recv(&value, sizeof(T), MPI_BYTE, _rank, BALANCE_MINIMAL_RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                heap.push(_T_Wrapper(value, _rank));
                break;
            }
            default:
                // should not reach here
                std::cerr << "Invalid tag: " << status.MPI_TAG << " received to rank " << this->rank << " from rank " << status.MPI_SOURCE << std::endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                break;
        }
    }
    return heap;
}

template<typename T>
void BalanceJob<T>::heapReplaceMin(BalanceJob<T>::Heap &heap)
{
    if(heap.empty())
    {
        std::cerr << "Error: heap is empty" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    _T_Wrapper full_value = heap.top();
    heap.pop();
    if(full_value.rank == this->rank)
    {
        if(this->values.size() <= this->heapPos)
        {
            return;
        }
        heap.push(_T_Wrapper(this->values[this->heapPos++], this->rank));
    }
    else
    {
        int dummy = 0;
        MPI_Send(&dummy, 1, MPI_INT, full_value.rank, BALANCE_MINIMAL_REQUEST_TAG, MPI_COMM_WORLD);
        MPI_Status status;
        int received = 0;
        while(!received)
        {
            MPI_Iprobe(full_value.rank, BALANCE_MINIMAL_RESPONSE_TAG, MPI_COMM_WORLD, &received, &status);
            if(received) break;
            MPI_Iprobe(full_value.rank, BALANCE_STOP_ITERATION_TAG, MPI_COMM_WORLD, &received, &status);
        }
        switch(status.MPI_TAG)
        {
            case BALANCE_MINIMAL_RESPONSE_TAG:
            {
                T value;
                MPI_Recv(&value, sizeof(T), MPI_BYTE, full_value.rank, BALANCE_MINIMAL_RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                heap.push(_T_Wrapper(value, full_value.rank));
                break;
            }
            case BALANCE_STOP_ITERATION_TAG:
            {
                // no value for rank, just leave it
                MPI_Recv(&dummy, 1, MPI_INT, full_value.rank, BALANCE_STOP_ITERATION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break;
            }
            default:
                // should not reach here
                std::cerr << "Invalid tag: " << status.MPI_TAG << " received to rank " << this->rank << " from rank " << status.MPI_SOURCE << std::endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                break;
        }
    }
}

template<typename T>
bool BalanceJob<T>::answerToRequests()
{
    MPI_Status status;
    int dummy;

    int finished = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, BALANCE_FINISHED_TAG, MPI_COMM_WORLD, &finished, &status);
    if(finished)
    {
        MPI_Recv(&dummy, 1, MPI_INT, status.MPI_SOURCE, BALANCE_FINISHED_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return true;
    }

    int received = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, BALANCE_MINIMAL_REQUEST_TAG, MPI_COMM_WORLD, &received, &status);
    while(received)
    {
        MPI_Recv(&dummy, 1, MPI_INT, status.MPI_SOURCE, BALANCE_MINIMAL_REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if(this->heapPos >= this->values.size())
        {
            int dummy = 0;
            MPI_Send(&dummy, 1, MPI_INT, status.MPI_SOURCE, BALANCE_STOP_ITERATION_TAG, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&this->values[this->heapPos++], sizeof(T), MPI_BYTE, status.MPI_SOURCE, BALANCE_MINIMAL_RESPONSE_TAG, MPI_COMM_WORLD);
        }
        break;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &received, &status);
    }
    return false;
}

template<typename T>
void BalanceJob<T>::announceBound(const T &bound)
{
    for(int _rank = 0; _rank < this->size; _rank++)
    {
        MPI_Send(&bound, sizeof(T), MPI_BYTE, _rank, BALANCE_NEW_BOUND_TAG, MPI_COMM_WORLD);
    }
}

template<typename T>
void BalanceJob<T>::kthOrderStatistics(std::vector<size_t> &orders, std::vector<T> &bounds)
{
    this->heapPos = 0;
    this->helpBuildHeap();

    if(this->rank == BALANCE_ROOT)
    {
        std::sort(orders.begin(), orders.end());
        auto heap = this->buildHeap();
        size_t i = 0, j = 0;
        while(i != orders.size())
        {
            // replace minimum in heap, until we reached orders[i]-th element
            while(j != orders[i])
            {
                this->heapReplaceMin(heap);
                j++;
            }
            const _T_Wrapper &value = heap.top();
            this->announceBound(value.value);
            i++;
        }
        // done all! send finish messages
        for(int _rank = 0; _rank < this->size; _rank++)
        {
            if(_rank != this->rank)
            {
                int dummy = 0;
                MPI_Send(&dummy, 1, MPI_INT, _rank, BALANCE_FINISHED_TAG, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        // answer to queries
        bool finished = false;
        while(!finished)
        {
            this->getNewBounds(bounds);
            finished = this->answerToRequests();
        }
    }
}

template<typename T>
void BalanceJob<T>::getNewBounds(std::vector<T> &bounds)
{
    // checks if a new bound has arrived. If yes, adds it
    if(bounds.size() == this->size)
    {
        return;
    }
    MPI_Status status;
    int arrived = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, BALANCE_NEW_BOUND_TAG, MPI_COMM_WORLD, &arrived, &status);
    while(arrived)
    {
        T bound;
        MPI_Recv(&bound, sizeof(T), MPI_BYTE, status.MPI_SOURCE, BALANCE_NEW_BOUND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bounds.push_back(bound);
        MPI_Iprobe(MPI_ANY_SOURCE, BALANCE_NEW_BOUND_TAG, MPI_COMM_WORLD, &arrived, &status);
    }
}

template<typename T>
std::vector<T> BalanceJob<T>::getBorders()
{
    std::vector<size_t> currLengths(this->size);
    size_t mySize = this->values.size();
    MPI_Allgather(&mySize, sizeof(size_t), MPI_BYTE, &currLengths[0], sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);
    size_t totalSize = 0;
    for(int i = 0; i < this->size; i++) totalSize += currLengths[i];
    // size_t requestedSize = (totalSize / this->size) + ((this->rank == this->size - 1)? (totalSize - (totalSize / this->size) * this->size) : 0); // if `totalSize` is not divisible, the last one gets the rest
    std::vector<size_t> stats(this->size);
    for(int i = 0; i < this->size - 1; i++)
    {
        stats[i] = (totalSize / this->size) * (i + 1);
    }
    stats[this->size - 1] = totalSize - 1;

    std::vector<T> bounds;

    this->kthOrderStatistics(stats, bounds);

    while(bounds.size() != this->size)
    {
        this->getNewBounds(bounds);
        usleep(10);
    }
    return bounds;
}

