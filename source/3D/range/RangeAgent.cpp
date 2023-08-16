#include "RangeAgent.h"

RangeAgent::RangeAgent(MPI_Comm comm, const HilbertAgent &hilbertAgent, RangeFinder *rangeFinder):
        comm(comm), hilbertAgent(hilbertAgent), rangeFinder(rangeFinder), hilbertTree(nullptr), myCurrentRound(0), finishedForNextRound(0)
{
    MPI_Comm_rank(this->comm, &this->rank);
    MPI_Comm_size(this->comm, &this->size);
}

void RangeAgent::receiveQueries(QueryBatchInfo &batch)
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

    MPI_Iprobe(MPI_ANY_SOURCE, TAG_RESPONSE, this->comm, &receivedAnswer, &status);

    std::vector<char> buffer;

    while(receivedAnswer and received < MAX_RECEIVE_IN_CYCLE)
    {
        ++this->receivedUntilNow;
        ++received;

        int count;
        MPI_Get_count(&status, MPI_BYTE, &count);
        if(buffer.size() < static_cast<size_t>(count))
        {
            buffer.resize(count);
        }

        MPI_Recv(&(*(buffer.begin())), count, MPI_BYTE, status.MPI_SOURCE, TAG_RESPONSE, this->comm, MPI_STATUS_IGNORE);
        size_t id;
        long int length;
        int pos = 0;
        MPI_Unpack(&(*(buffer.begin())), count, &pos, &id, 1, MPI_UNSIGNED_LONG, this->comm);
        MPI_Unpack(&(*(buffer.begin())), count, &pos, &length, 1, MPI_LONG, this->comm);
        if(length > 0)
        {
            // insert the results to the points received by rank `status.MPI_SOURCE` and to the queries result
            queries[id].finalResults.resize(queries[id].finalResults.size() + length);
            MPI_Unpack(&(*(buffer.begin())), count, &pos, &(*(queries[id].finalResults.end() - length)), length * sizeof(_3DPoint), MPI_BYTE, this->comm);
            pointsFromRanks[status.MPI_SOURCE].resize(pointsFromRanks[status.MPI_SOURCE].size() + length);
            for(size_t i = 0; i < static_cast<size_t>(length); i++)
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
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_RESPONSE, this->comm, &receivedAnswer, &status);
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

    MPI_Iprobe(MPI_ANY_SOURCE, TAG_REQUEST, this->comm, &arrivedNew, &status);

    // while arrived new messages, and we should answer until the end, or answer until a bound we haven't reached to
    while(arrivedNew != 0 and (finishAnswering or (!finishAnswering and answered < MAX_ANSWER_IN_CYCLE)))
    {
        SubQueryData query;
        MPI_Recv(&query, sizeof(SubQueryData), MPI_BYTE, status.MPI_SOURCE,  TAG_REQUEST, this->comm, MPI_STATUS_IGNORE);
        answered++;
        // calculate the result
        std::vector<Vector3D> result = this->getRangeResult(query, status.MPI_SOURCE);
        long int resultSize = static_cast<long int>(result.size());

        int pos = 0;
        this->buffers.push_back(std::vector<char>());
        std::vector<char> &to_send = this->buffers[this->buffers.size() - 1];
        size_t msg_size = 2 * sizeof(size_t) + resultSize * sizeof(_3DPoint);
        to_send.resize(msg_size);

        MPI_Pack(&query.parent_id, 1, MPI_UNSIGNED_LONG, &to_send[0], msg_size, &pos, this->comm);
        MPI_Pack(&resultSize, 1, MPI_LONG, &to_send[0], msg_size, &pos, this->comm);
        if(resultSize > 0)
        {
            _3DPoint *points = reinterpret_cast<_3DPoint*>((&to_send[0]) + pos);
            for(size_t i = 0; i < static_cast<size_t>(resultSize); i++)
            {
                points[i] = {result[i].x, result[i].y, result[i].z};
            }
        }
        this->requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&to_send[0], msg_size, MPI_BYTE, status.MPI_SOURCE, TAG_RESPONSE, this->comm, &this->requests[requests.size() - 1]);

        MPI_Iprobe(MPI_ANY_SOURCE,  TAG_REQUEST, this->comm, &arrivedNew, &status);
    }
}

/**
 * gets a center and a radius, and returns all the ranks which the sphere with cener `center` and radius `radius`, intersects with their areas.
*/
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
        this->buffers.push_back(std::vector<char>(sizeof(SubQueryData)));
        SubQueryData &subQuery = *reinterpret_cast<SubQueryData*>(&this->buffers[this->buffers.size() - 1][0]);
        subQuery.data = query.data;
        subQuery.parent_id = query.id;
        this->requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&subQuery, sizeof(SubQueryData), MPI_BYTE, node, TAG_REQUEST, this->comm, &this->requests[this->requests.size() - 1]);
    }
}

void RangeAgent::sendFinish()
{
    int dummy;
    for(int _rank = 0; _rank < this->size; _rank++)
    {
        this->requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&dummy, 1, MPI_BYTE, _rank, TAG_FINISHED, this->comm, &this->requests[this->requests.size() - 1]);
    }
}

int RangeAgent::checkForFinishMessages() const
{
    int arrived = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, TAG_FINISHED, this->comm, &arrived, &status);
    if(arrived)
    {
        int dummy = 0;
        MPI_Recv(&dummy, 1, MPI_BYTE, MPI_ANY_SOURCE, TAG_FINISHED, this->comm, MPI_STATUS_IGNORE);
        return 1;
    }
    return 0;
}

// todo: change the parallel model: each one sends a message when it finishes sending requests. Waits to `size` such messages, then answers all again
QueryBatchInfo RangeAgent::runBatch(std::queue<RangeQueryData> &queries)
{
    this->receivedUntilNow = 0; // reset the receive counter
    this->shouldReceiveInTotal = 0; // reset the should-be-received counter
    this->buffers.clear();
    this->buffers.reserve(3 * queries.size()); // heuristic
    this->requests.clear();
    QueryBatchInfo queriesBatch;
    std::vector<QueryInfo> &queriesInfo = queriesBatch.queriesAnswers;

    int finishedReceived = 0;
    bool sentFinished = false;
    size_t i = 0;
    while((finishedReceived < this->size) or (!queries.empty()))
    {
        if(!queries.empty())
        {
            RangeQueryData queryData = queries.front();
            queries.pop();
            queriesInfo.push_back({queryData, i, std::vector<_3DPoint>()});
            QueryInfo &query = queriesInfo[queriesInfo.size() - 1];
            this->sendQuery(query);
        }
        else
        {
            if(i % FINISH_AUTOFLUSH_NUM == 0)
            {
                if((this->shouldReceiveInTotal <= this->receivedUntilNow) and !sentFinished)
                {
                    sentFinished = true;
                    this->sendFinish();
                }
                finishedReceived += this->checkForFinishMessages();
            }
        }
        if(i % QUERY_AUTOFLUSH_NUM == 0)
        {
            this->answerQueries(false);
        }
        if(i % RECEIVE_AUTOFLUSH_NUM == 0)
        {
            this->receiveQueries(queriesBatch);
        }
        ++i;
    }
    if(this->requests.size() > 0)
    {
        MPI_Waitall(this->requests.size(), &(*(this->requests.begin())), MPI_STATUSES_IGNORE); // make sure any query was indeed received
    }
    return queriesBatch;

}
