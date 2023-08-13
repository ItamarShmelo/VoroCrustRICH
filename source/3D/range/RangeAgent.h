#ifndef _RICH_RANGE_AGENT_H_
#define _RICH_RANGE_AGENT_H_

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

#define QUERY_AUTOFLUSH_NUM 5
#define RECEIVE_AUTOFLUSH_NUM 10
#define FINISH_AUTOFLUSH_NUM 20
#define MAX_RECEIVE_IN_CYCLE 1000
#define MAX_ANSWER_IN_CYCLE 1000

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
    std::vector<std::vector<Vector3D>> pointsFromRanks;
    std::vector<Vector3D> newPoints;
} QueryBatchInfo;


class RangeAgent
{
    template<typename T>
    using _set = boost::container::flat_set<T>;

public:
    RangeAgent(MPI_Comm comm, const HilbertAgent &hilbertAgent, RangeFinder *rangeFinder);
    inline RangeAgent(const HilbertAgent &hilbertAgent, RangeFinder *rangeFinder): RangeAgent(MPI_COMM_WORLD, hilbertAgent, rangeFinder){};
    inline RangeAgent(MPI_Comm comm, const Vector3D &origin, const Vector3D &corner, int order, RangeFinder *rangeFinder): RangeAgent(comm, HilbertAgent(origin, corner, order), rangeFinder){};
    inline RangeAgent(const Vector3D &origin, const Vector3D &corner, int order, RangeFinder *rangeFinder): RangeAgent(MPI_COMM_WORLD, origin, corner, order, rangeFinder){};
    inline ~RangeAgent(){delete this->hilbertTree;};

    void receiveQueries(QueryBatchInfo &batch, bool blocking);
    void answerQueries(bool finishAnswering);
    void sendQuery(const QueryInfo &query);
    QueryBatchInfo runBatch(std::queue<RangeQueryData> &queries);
    void createArtificialQueries(coord_t radius);

    inline size_t getNumPoints(){return this->rangeFinder->size();};
    inline void buildHilbertTree(const OctTree<Vector3D> *tree){this->hilbertTree = new DistributedOctTree(tree);};

private:
    MPI_Comm comm;
    int rank, size;
    int order;
    std::vector<MPI_Request> requests;
    std::vector<std::vector<char>> buffers;
    size_t receivedUntilNow;
    size_t shouldReceiveInTotal;
    HilbertAgent hilbertAgent;
    RangeFinder *rangeFinder;
    DistributedOctTree<Vector3D> *hilbertTree;
    int myCurrentRound;
    int finishedForNextRound;

    std::vector<int> sentProcessorsRanks;
    std::vector<_set<size_t>> sentPoints; 
    std::vector<int> recvProcessorsRank;
    std::vector<std::vector<size_t>> recvPoints; 
    
    std::vector<Vector3D> getRangeResult(const SubQueryData &query, int node);
    _set<int> getIntersectingRanks(const Vector3D &center, coord_t radius) const;
    void sendFinish(const int &dummy);
    int checkForFinishMessages() const;
};

#endif // _RICH_RANGE_AGENT_H_
