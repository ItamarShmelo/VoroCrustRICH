#ifndef _RICH_RANGE_AGENT_H_
#define _RICH_RANGE_AGENT_H_

#ifdef RICH_MPI

#include <iostream> // todo remove
#include <algorithm>
#include <cmath>
#include <set>
#include <queue>
#include <vector>
#include <mpi.h>

// set data structure:
#include <boost/container/flat_set.hpp>
#include <unordered_set>

#include "ds/DistributedOctTree/DistributedOctTree.hpp"
#include "../hilbert/HilbertAgent.h"
#include "finders/RangeFinder.hpp"

#define TAG_REQUEST 200
#define TAG_RESPONSE 201
#define TAG_FINISHED 202

#define UNDEFINED_BUFFER_IDX -1
#define FLUSH_QUERIES_NUM 80

#define QUERY_AUTOFLUSH_NUM 100
#define RECEIVE_AUTOFLUSH_NUM 15
#define FINISH_AUTOFLUSH_NUM 50
#define MAX_RECEIVE_IN_CYCLE 20
#define MAX_ANSWER_IN_CYCLE 20

typedef struct RangeQueryData
{
    _3DPoint center;
    coord_t radius;
} RangeQueryData;

typedef struct SubQueryData
{
    RangeQueryData data;
    size_t parent_id;
} SubQueryData;

typedef struct QueryInfo
{
    RangeQueryData data;
    size_t id;
    // int subQueriesNum;
    std::vector<_3DPoint> finalResults;
} QueryInfo;

typedef struct QueryBatchInfo
{
    std::vector<QueryInfo> queriesAnswers;
    std::vector<Vector3D> newPoints;
} QueryBatchInfo;


class RangeAgent
{
public:
    template<typename T>
    using _set = boost::container::flat_set<T>;

    RangeAgent(MPI_Comm comm, const HilbertAgent &hilbertAgent, RangeFinder *rangeFinder);
    inline RangeAgent(const HilbertAgent &hilbertAgent, RangeFinder *rangeFinder): RangeAgent(MPI_COMM_WORLD, hilbertAgent, rangeFinder){};
    inline RangeAgent(MPI_Comm comm, const Vector3D &origin, const Vector3D &corner, int order, RangeFinder *rangeFinder): RangeAgent(comm, HilbertAgent(origin, corner, order), rangeFinder){};
    inline RangeAgent(const Vector3D &origin, const Vector3D &corner, int order, RangeFinder *rangeFinder): RangeAgent(MPI_COMM_WORLD, origin, corner, order, rangeFinder){};
    inline ~RangeAgent(){/*delete this->hilbertTree;*/};

    void receiveQueries(QueryBatchInfo &batch);
    void answerQueries();
    void sendQuery(const QueryInfo &query);
    QueryBatchInfo runBatch(std::queue<RangeQueryData> &queries);

    inline size_t getNumPoints(){return this->rangeFinder->size();};
    inline void buildHilbertTree(const OctTree<Vector3D> *tree){this->hilbertTree = new DistributedOctTree(tree);};

    std::vector<std::vector<size_t>> &getSentPoints(){return this->sentPoints;};
    std::vector<std::vector<size_t>> &getRecvPoints(){return this->recvPoints;};
    std::vector<int> &getSentProc(){return this->sentProcessorsRanks;};
    std::vector<int> &getRecvProc(){return this->recvProcessorsRanks;};

private:
    MPI_Comm comm;
    int rank, size;
    std::vector<MPI_Request> requests;
    std::vector<std::vector<char>> buffers;
    size_t receivedUntilNow;
    size_t shouldReceiveInTotal;
    HilbertAgent hilbertAgent;
    RangeFinder *rangeFinder;
    DistributedOctTree<Vector3D> *hilbertTree;

    std::vector<int> sentProcessorsRanks;
    std::vector<int> recvProcessorsRanks;
    std::vector<std::vector<size_t>> recvPoints; 
    std::vector<std::vector<size_t>> sentPoints; 
    std::vector<_set<size_t>> sentPointsSet; 
    std::vector<int> ranksBufferIdx;
    
    std::vector<Vector3D> getRangeResult(const SubQueryData &query, int rank);
    _set<int> getIntersectingRanks(const Vector3D &center, coord_t radius) const;
    void sendFinish();
    int checkForFinishMessages() const;
    inline void flushBuffer(int node)
    {
        int bufferIdx = this->ranksBufferIdx[node];
        if(bufferIdx == UNDEFINED_BUFFER_IDX)
        {
            return;
        }
        if(this->buffers[bufferIdx].size() > 0)
        {
            this->requests.push_back(MPI_REQUEST_NULL);
            MPI_Isend(&this->buffers[bufferIdx][0], this->buffers[bufferIdx].size(), MPI_BYTE, node, TAG_REQUEST, this->comm, &this->requests[this->requests.size() - 1]);
        }
        this->ranksBufferIdx[node] = UNDEFINED_BUFFER_IDX;
    }
    inline void flushAll(){for(int i = 0; i < this->size; i++) this->flushBuffer(i);};
};

#endif // RICH_MPI

#endif // _RICH_RANGE_AGENT_H_
