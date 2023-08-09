#include "RangeAgent.h"

RangeAgent::RangeAgent(MPI_Comm comm, const HilbertAgent &hilbertAgent, RangeFinder *rangeFinder):
        globalComm(comm), hilbertAgent(hilbertAgent), rangeFinder(rangeFinder), hilbertTree(nullptr)
{
    MPI_Comm_rank(this->globalComm, &this->rank);
    MPI_Comm_size(this->globalComm, &this->size);
    MPI_Comm_dup(this->globalComm, &this->answersComm);
    MPI_Comm_dup(this->globalComm, &this->requestsComm);
    MPI_Comm_dup(this->globalComm, &this->finishedComm);
}

void RangeAgent::receiveQueries(QueryBatchInfo &batch, bool blocking)
{
    if(this->receivedUntilNow >= this->shouldReceiveInTotal)
    {
        return;
    }
    MPI_Status status;
    int receivedAnswer = 0;
    int received = 0;

    std::vector<QueryInfo> &queries = batch.queriesAnswers;
    std::vector<std::vector<Vector3D>> &pointsFromRanks = batch.pointsFromRanks;
    pointsFromRanks.resize(this->size);

    if(blocking)
    {
        MPI_Probe(MPI_ANY_SOURCE, /*TAG_RESPONSE*/ MPI_ANY_TAG, this->answersComm, &status);
    }
    else
    {
        MPI_Iprobe(MPI_ANY_SOURCE, /*TAG_RESPONSE*/ MPI_ANY_TAG, this->answersComm, &receivedAnswer, &status);
    }

    std::vector<char> buffer;

    while((!blocking and receivedAnswer and received < MAX_RECEIVE_IN_CYCLE) or (blocking and this->receivedUntilNow < this->shouldReceiveInTotal))
    {
        ++this->receivedUntilNow;
        ++received;

        int count;
        MPI_Get_count(&status, MPI_BYTE, &count);
        if(buffer.size() < static_cast<size_t>(count))
        {
            buffer.resize(count);
        }

        MPI_Recv(&(*(buffer.begin())), count, MPI_BYTE, status.MPI_SOURCE, /*TAG_RESPONSE*/ MPI_ANY_TAG, this->answersComm, MPI_STATUS_IGNORE);
        size_t id;
        long int length;
        int pos = 0;
        MPI_Unpack(&(*(buffer.begin())), count, &pos, &id, 1, MPI_UNSIGNED_LONG, this->answersComm);
        MPI_Unpack(&(*(buffer.begin())), count, &pos, &length, 1, MPI_LONG, this->answersComm);
        if(length > 0)
        {
            // std::cout << "id is " << id << std::endl;
            // insert the results to the points received by rank `status.MPI_SOURCE` and to the queries result
            queries[id].finalResults.resize(queries[id].finalResults.size() + length);
            MPI_Unpack(&(*(buffer.begin())), count, &pos, &(*(queries[id].finalResults.end() - length)), length * sizeof(_3DPoint), MPI_BYTE, this->answersComm);
            pointsFromRanks[status.MPI_SOURCE].resize(pointsFromRanks[status.MPI_SOURCE].size() + length);
            for(size_t i = 0; i < length; i++)
            {
                _3DPoint &point = *(queries[id].finalResults.end() - length + i);
                Vector3D vector(point.x, point.y, point.z);
                pointsFromRanks[status.MPI_SOURCE].push_back(vector);
                batch.newPoints.push_back(vector);
            }
        }
        else
        {
            assert(length >= 0);
        }
        if(blocking)
        {
            if(this->receivedUntilNow < this->shouldReceiveInTotal)
            {
                MPI_Probe(MPI_ANY_SOURCE, /*TAG_RESPONSE*/ MPI_ANY_TAG, this->answersComm, &status);
            }
            else
            {
                break; // finished all
            }
        }
        else
        {
            MPI_Iprobe(MPI_ANY_SOURCE, /*TAG_RESPONSE*/ MPI_ANY_TAG, this->answersComm, &receivedAnswer, &status);
        }
    }
}

/**
 * gets a query, and who requests it, and returns the answer that should be sent.
*/
std::vector<Vector3D> RangeAgent::getRangeResult(const SubQueryData &query, int node)
{
    // get what is the right index of `node` inside the sentProcessors vector. If it isn't there, create it
    std::vector<Vector3D> result;
    size_t rankIndex = std::find(this->sentProcessorsRanks.begin(), this->sentProcessorsRanks.end(), node) - this->sentProcessorsRanks.begin();
    if(rankIndex == this->sentProcessorsRanks.size())
    {
        // rank is not inside the sentProcessors rank, add it
        this->sentProcessorsRanks.push_back(node);
        this->sentPoints.push_back(_set<size_t>());
    }

    // get the real results, and filter it (do not send points you sent before)
    std::vector<size_t> nonFilteredResult = this->rangeFinder->range(Vector3D(query.data.center.x, query.data.center.y, query.data.center.z), query.data.radius);
    for(const size_t &pointIdx : nonFilteredResult)
    {
        // check if the point wasn't already sent to node, only after that, add it to the result
        //if(std::find(this->sentPoints[rankIndex].begin(), this->sentPoints[rankIndex].end(), pointIdx) == this->sentPoints[rankIndex].end())
        if(this->sentPoints[rankIndex].find(pointIdx) == this->sentPoints[rankIndex].end())
        {
            // point haven't been sent, send it
            result.push_back(this->rangeFinder->getPoint(pointIdx));
            this->sentPoints[rankIndex].insert(pointIdx);
        }
    }
    return result;
}

void RangeAgent::answerQueries(bool finishAnswering)
{
    MPI_Status status;
    int arrivedNew = 0;
    int answered = 0;

    MPI_Iprobe(MPI_ANY_SOURCE,  /*TAG_REQUEST*/ MPI_ANY_TAG, this->requestsComm, &arrivedNew, &status);

    while(arrivedNew != 0 and (finishAnswering or (!finishAnswering and answered < MAX_ANSWER_IN_CYCLE)))
    {
        SubQueryData query;
        MPI_Recv(&query, sizeof(SubQueryData), MPI_BYTE, status.MPI_SOURCE,  /*TAG_REQUEST*/ MPI_ANY_TAG, this->requestsComm, MPI_STATUS_IGNORE);
        answered++;
        std::vector<Vector3D> result = this->getRangeResult(query, status.MPI_SOURCE);
        long int resultSize = static_cast<long int>(result.size());

        int pos = 0;
        this->buffers.push_back(std::vector<char>());
        std::vector<char> &to_send = this->buffers[this->buffers.size() - 1];
        size_t msg_size = 2 * sizeof(size_t) + resultSize * sizeof(_3DPoint);
        to_send.resize(msg_size);

        MPI_Pack(&query.parent_id, 1, MPI_UNSIGNED_LONG, &to_send[0], msg_size, &pos, this->requestsComm);
        MPI_Pack(&resultSize, 1, MPI_LONG, &to_send[0], msg_size, &pos, this->requestsComm);
        if(resultSize > 0)
        {
            _3DPoint *points = reinterpret_cast<_3DPoint*>((&to_send[0]) + pos);
            for(size_t i = 0; i < resultSize; i++)
            {
                points[i] = {result[i].x, result[i].y, result[i].z};
            }
        }
        this->requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&to_send[0], msg_size, MPI_BYTE, status.MPI_SOURCE, TAG_RESPONSE, this->answersComm, &this->requests[requests.size() - 1]);

        MPI_Iprobe(MPI_ANY_SOURCE,  /*TAG_REQUEST*/ MPI_ANY_TAG, this->requestsComm, &arrivedNew, &status);
    }
}

typename RangeAgent::_set<int> RangeAgent::getIntersectingRanks(const Vector3D &center, coord_t radius) const
{
    _set<int> possibleNodes;

    if(this->hilbertTree == nullptr)
    {
        auto intersectionHilbertCells = this->hilbertAgent.getIntersectingCircle(center, radius);
        for(const hilbert_index_t &index : intersectionHilbertCells)
        {
            possibleNodes.insert(this->hilbertAgent.getCellOwner(index));
        }
    }
    else
    {
        for(const int &rank : this->hilbertTree->getIntersectingRanks(center, radius))
        {
            possibleNodes.insert(rank);
        }
    }

    return possibleNodes;
}

void RangeAgent::sendQuery(const QueryInfo &query)
{
    _set<int> possibleNodes = this->getIntersectingRanks(Vector3D(query.data.center.x, query.data.center.y, query.data.center.z), query.data.radius);
    for(const int &node : possibleNodes)
    {
        if(node == this->rank)
        {
            continue; // unnecessary to send
        }
        ++this->shouldReceiveInTotal;
        SubQueryData subQuery = {query.data, query.id};
        // this->requests.push_back(MPI_REQUEST_NULL);
        MPI_Send(&subQuery, sizeof(SubQueryData), MPI_BYTE, node, TAG_REQUEST, this->requestsComm/*, &this->requests[this->requests.size() - 1]*/);
    }
}

// todo: change the parallel model: each one sends a message when it finishes sending requests. Waits to `size` such messages, then answers all again
QueryBatchInfo RangeAgent::runBatch(std::queue<RangeQueryData> &queries)
{
    this->receivedUntilNow = 0; // reset the receive counter
    this->shouldReceiveInTotal = 0; // reset the should-be-received counter
    this->buffers.clear();
    this->requests.clear();

    QueryBatchInfo queriesBatch;
    std::vector<QueryInfo> &queriesInfo = queriesBatch.queriesAnswers;

    size_t i = 0;
    while(!queries.empty())
    {
        RangeQueryData queryData = queries.front();
        queries.pop();
        queriesInfo.push_back({queryData, i, std::vector<_3DPoint>()});
        QueryInfo &query = queriesInfo[queriesInfo.size() - 1];
        this->sendQuery(query);
        ++i;
        if(i % QUERY_AUTOFLUSH_NUM == 0)
        {
            this->answerQueries(false);
        }
        if(i % RECEIVE_AUTOFLUSH_NUM == 0)
        {
            this->receiveQueries(queriesBatch, false);
        }
    }
    MPI_Barrier(this->globalComm); // ensure everyone stopped sending
    this->answerQueries(true); // answer the arrived queries
    this->receiveQueries(queriesBatch, true); // receive the remain results
    // MPI_Barrier(this->globalComm);
    if(this->requests.size() > 0)
    {
        MPI_Waitall(this->requests.size(), &(*(this->requests.begin())), MPI_STATUSES_IGNORE); // make sure any query was indeed received
    }
    // this->answerQueries(true);
    return queriesBatch;

}
