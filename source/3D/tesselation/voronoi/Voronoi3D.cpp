#include "Voronoi3D.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif
#include <algorithm>
#include <cfloat>
#include <stack>
#include "../../elementary/Mat33.hpp"
#include "../utils/Predicates3D.hpp"
#include "misc/utils.hpp"
#include "misc/io3D.hpp"
#include <fstream>
#include <iostream>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include "3D/GeometryCommon/Intersections.hpp"
#ifdef RICH_MPI
#include "3D/range/finders/BruteForce.hpp"
#include "3D/range/finders/RangeTree.hpp"
#include "3D/range/finders/OctTree.hpp"
#include "3D/range/finders/SmartBruteForce.hpp"
#include "3D/range/finders/HashBruteForce.hpp"
#include "3D/range/finders/GroupRangeTree.hpp"
#endif // RICH_MPI
#include "misc/int2str.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/container/static_vector.hpp>
#include <omp.h>

// #define VORONOI_DEBUG

bool PointInPoly(Tessellation3D const &tess, Vector3D const &point, std::size_t index)
{
    face_vec const &faces = tess.GetCellFaces(index);
    vector<Vector3D> const &points = tess.GetFacePoints();
    std::size_t N = faces.size();
    std::array<Vector3D, 4> vec;
    for (std::size_t i = 0; i < N; ++i)
    {
        double R = fastsqrt(tess.GetArea(faces[i]));
        size_t N1 = 0;
        size_t N2 = 0;
        Vector3D V1, V2;
        size_t counter = 0;
        point_vec const &InFace = tess.GetPointsInFace(faces[i]);
        size_t NinFace = InFace.size();
        N1 = 1;
        V1 = points[InFace[(counter + 1) % NinFace]] - points[InFace[0]];
        while (fastabs(V1) < 0.01 * R)
        {
            ++counter;
            assert(counter < NinFace);
            V1 = points[InFace[(counter + 1) % NinFace]] - points[InFace[0]];
            ++N1;
        }
        V2 = points[InFace[(counter + 2) % NinFace]] - points[InFace[N1]];
        N2 = (counter + 2) % NinFace;
        while (fastabs(V2) < 0.01 * R || fastabs(CrossProduct(V1, V2)) < 0.0001 * tess.GetArea(faces[i]))
        {
            ++counter;
            if (counter > 2 * NinFace)
                break;
            V2 = points[InFace[(counter + 2) % NinFace]] - points[InFace[N1]];
            N2 = (counter + 2) % NinFace;
        }
        if (counter > 2 * NinFace)
        {
            std::cout << "Weird face in PointInPoly, cell " << index << " face " << faces[i] << " i " << i << " face area " << tess.GetArea(faces[i]) << std::endl;
            for (size_t j = 0; j < NinFace; ++j)
                std::cout << "Point j " << points[InFace[j]].x << "," << points[InFace[j]].y << "," << points[InFace[j]].z << std::endl;
            Vector3D normal = tess.GetFaceNeighbors(faces[i]).second == index ? tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).second) - tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).first) : tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).first) - tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).second);
            if (ScalarProd(normal, point - points[InFace[0]]) < 0)
                return false;
        }
        else
        {
            vec[0] = points[InFace[0]];
            vec[1] = points[InFace.at(N1)];
            vec[2] = points[InFace.at(N2)];
            vec[3] = tess.GetMeshPoint(index);
            double s1 = orient3d(vec);
            vec[3] = point;
            double s2 = orient3d(vec);
            if (s1 * s2 < -0)
                return false;
        }
    }
    return true;
}

bool PointInPoly(std::vector<Face> const& faces, Vector3D const &point)
{
    std::size_t const N = faces.size();
    std::array<Vector3D, 4> vec;
    vec[3] = point;
    for (std::size_t i = 0; i < N; ++i)
    {
        vec[0] = faces[i].vertices[0];
        vec[1] = faces[i].vertices[1];
        vec[2] = faces[i].vertices[2];
        double const s = orient3d(vec);
        if (s > 0)
            return false;
    }
    return true;
}

namespace
{
#ifdef RICH_MPI
    void GetPastDuplicate(size_t point, vector<size_t> &res, vector<vector<size_t>> const &sorted_to_duplicate,
                                                vector<size_t> const &procs)
    {
        res.clear();
        for (size_t i = 0; i < procs.size(); ++i)
        {
            if (std::binary_search(sorted_to_duplicate[i].begin(), sorted_to_duplicate[i].end(), point))
                res.push_back(procs[i]);
        }
    }
#endif
    boost::multiprecision::cpp_dec_float_50 Calc33Det(std::array<boost::multiprecision::cpp_dec_float_50, 9> const &points)
    {
        return points[0] * (points[4] * points[8] - points[5] * points[7]) + points[1] * (points[5] * points[6] - points[3] * points[8]) + points[2] * (points[3] * points[7] - points[4] * points[6]);
    }
}

namespace
{
    bool ShouldCalcTetraRadius(Tetrahedron const &T, size_t Norg)
    {
        for (size_t i = 0; i < 4; ++i)
            if (T.points[i] < Norg)
                return true;
        return false;
    }

    void FirstCheckList(std::stack<std::size_t> &check_stack, vector<unsigned char> &future_check, size_t Norg,
                                            Delaunay3D const &del, vector<tetra_vec> const &PointsInTetra)
    {
        //        check_stack.empty();
        future_check.resize(Norg, 0);
        size_t Ntetra = del.tetras_.size();
        vector<unsigned char> tetra_check(Ntetra, 0);

        for (size_t i = 0; i < Ntetra; ++i)
        {
            Tetrahedron const &tetra = del.tetras_[i];
            for (size_t j = 0; j < 4; ++j)
            {
                if (tetra.points[j] >= Norg)
                {
                    for (size_t k = 0; k < 4; ++k)
                    {
                        size_t tetcheck = tetra.points[k];
                        if (tetra.points[k] < Norg)
                        {
                            size_t ntet = PointsInTetra[tetcheck].size();
                            for (size_t z = 0; z < ntet; ++z)
                                tetra_check[PointsInTetra[tetcheck][z]] = 1;
                        }
                    }
                    break;
                }
            }
        }
        for (size_t i = 0; i < Ntetra; ++i)
        {
            if (tetra_check[i] == 1)
            {
                Tetrahedron const &tetra = del.tetras_[i];
                for (size_t j = 0; j < 4; ++j)
                {
                    if (tetra.points[j] < Norg)
                        future_check[tetra.points[j]] = 1;
                }
            }
        }
        for (size_t i = 0; i < Norg; ++i)
            if (future_check[i] == 1)
                check_stack.push(i);
    }

    vector<Face> BuildBox(Vector3D const &ll, Vector3D const &ur)
    {
        double dx = ur.x - ll.x;
        double dy = ur.y - ll.y;
        double dz = ur.z - ll.z;
        vector<Face> res(6);
        vector<Vector3D> points;
        points.push_back(ll);
        points.push_back(ll + Vector3D(dx, 0, 0));
        points.push_back(ll + Vector3D(dx, dy, 0));
        points.push_back(ll + Vector3D(0, dy, 0));
        points.push_back(ll + Vector3D(0, 0, dz));
        points.push_back(ll + Vector3D(dx, 0, dz));
        points.push_back(ll + Vector3D(dx, dy, dz));
        points.push_back(ll + Vector3D(0, dy, dz));
        points.push_back(ur);
        res[0].vertices.push_back(points[0]);
        res[0].vertices.push_back(points[1]);
        res[0].vertices.push_back(points[2]);
        res[0].vertices.push_back(points[3]);
        res[1].vertices.push_back(points[0]);
        res[1].vertices.push_back(points[4]);
        res[1].vertices.push_back(points[5]);
        res[1].vertices.push_back(points[1]);
        res[2].vertices.push_back(points[3]);
        res[2].vertices.push_back(points[7]);
        res[2].vertices.push_back(points[4]);
        res[2].vertices.push_back(points[0]);
        res[3].vertices.push_back(points[2]);
        res[3].vertices.push_back(points[6]);
        res[3].vertices.push_back(points[7]);
        res[3].vertices.push_back(points[3]);
        res[4].vertices.push_back(points[1]);
        res[4].vertices.push_back(points[5]);
        res[4].vertices.push_back(points[6]);
        res[4].vertices.push_back(points[2]);
        res[5].vertices.push_back(points[5]);
        res[5].vertices.push_back(points[4]);
        res[5].vertices.push_back(points[7]);
        res[5].vertices.push_back(points[6]);
        return res;
    }
#ifdef RICH_MPI
    vector<Vector3D> GetBoxNormals(Vector3D const &ll, Vector3D const &ur, vector<Face> const& box_faces_)
    {
        const vector<Face> faces = box_faces_.empty() ? BuildBox(ll, ur) : box_faces_;
        vector<Vector3D> res(faces.size());
        size_t N = res.size();
        for (size_t i = 0; i < N; ++i)
        {
            CrossProduct(faces[i].vertices[2] - faces[i].vertices[0], faces[i].vertices[1] - faces[i].vertices[0], res[i]);
            res[i] *= 1.0 /abs(res[i]);
        }
        return res;
    }

    size_t BoxIndex(vector<Vector3D> const &fnormals, const Vector3D &normal)
    {
        double max_angle = ScalarProd(fnormals[0], normal);
        size_t loc = 0;
        size_t N = fnormals.size();
        for (size_t i = 1; i < N; i++)
        {
            double temp = ScalarProd(fnormals[i], normal);
            if (temp > max_angle)
            {
                max_angle = temp;
                loc = i;
            }
        }
        return loc;
    }
#endif

    double CleanDuplicates(std::array<size_t, 128> const &indeces, const vector<Vector3D> &points,
                                                 boost::container::small_vector<size_t, 8> &res, double R,
                                                 std::array<double, 128> &diffs,
                                                 std::array<Vector3D, 128> &vtemp, const size_t N)
    {
        res.clear();
        for (size_t i = 0; i < N; ++i)
            vtemp[i] = points[indeces[i]];
        for (size_t i = N - 1; i > 0; --i)
        {
            vtemp[i].x -= vtemp[i - 1].x;
            vtemp[i].y -= vtemp[i - 1].y;
            vtemp[i].z -= vtemp[i - 1].z;
        }
        vtemp[0] -= points[indeces[N - 1]];
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(max \
                                                     : R)
#endif
        for (size_t i = 0; i < N; ++i)
        {
            diffs[i] = ScalarProd(vtemp[i], vtemp[i]);
            R = std::max(R, diffs[i]);
        }
        for (size_t i = 0; i < N; ++i)
            if (diffs[i] > R * 1e-16)
                res.push_back(indeces[i]);
        return R;
    }

    size_t SetPointTetras(vector<tetra_vec> &PointTetras, size_t Norg, vector<Tetrahedron> &tetras,
                                                boost::container::flat_set<size_t> const &empty_tetras)
    {
        PointTetras.clear();
        PointTetras.resize(Norg);
        size_t Ntetra = tetras.size();
        size_t bigtet(0);
        bool has_good, has_big;
        // change empty tetras to be not relevant
        for (boost::container::flat_set<size_t>::const_iterator it = empty_tetras.begin(); it !=
                                                                                                                                                                             empty_tetras.end();
                 ++it)
        {
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
#endif
            for (size_t i = 0; i < 4; ++i)
            {
                tetras[*it].points[i] = std::numeric_limits<std::size_t>::max();
                tetras[*it].neighbors[i] = std::numeric_limits<std::size_t>::max();
            }
        }

        for (size_t i = 0; i < Ntetra; ++i)
        {
            has_good = false;
            has_big = false;
            for (size_t j = 0; j < 4; ++j)
            {
                size_t temp = tetras[i].points[j];
                if (temp < Norg)
                {
                    PointTetras[temp].push_back(i);
                    has_good = true;
                }
                else
                    has_big = true;
            }
            if (has_big && has_good)
                bigtet = i;
        }
        return bigtet;
    }

    bool CleanSameLine(boost::container::small_vector<size_t, 8> &indeces, vector<Vector3D> const& face_points, std::array<double, 128> &area_vec_temp)
    {
        point_vec old;
        size_t const N = indeces.size();
        double const small_fraction = 1e-14;
        double const medium_fraction = 3e-1;
        // Find correct normal
        Vector3D good_normal;
        for(size_t i = 0; i < N; ++i)
        {
        area_vec_temp[i] = fastabs(CrossProduct(face_points[indeces[i]] - face_points[indeces[(N + i - 1) % N]], face_points[indeces[(i + 1) % N]] 
        - face_points[indeces[(N + i - 1) % N]]));
        old.push_back(indeces[i]);
        }
        double const area_scale = *std::max_element(area_vec_temp.begin(), area_vec_temp.begin() + N);
        good_normal = CrossProduct(face_points[indeces[0]] - face_points[indeces[N - 1]], face_points[indeces[1]] 
        - face_points[indeces[N - 1]]);
        good_normal *= 1.0 / (100 * std::numeric_limits<double>::min() + fastabs(good_normal));
        for(size_t i = 1; i < N; ++i)
        {
        if(area_vec_temp[i] > area_scale * medium_fraction)
        {
            good_normal = CrossProduct(face_points[indeces[i]] - face_points[indeces[i - 1]], face_points[indeces[(i + 1) % N]] 
            - face_points[indeces[i - 1]]);
            good_normal *= 1.0 / fastabs(good_normal);
            break;
        }
        }
        size_t Nindeces = indeces.size();
        for(size_t i = 0; i < Nindeces; ++i)
        {
        Vector3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
        double const area = fastabs(normal_temp);
        normal_temp *= 1.0 / (100 * std::numeric_limits<double>::min() + area);
        if(area < area_scale * small_fraction || ScalarProd(normal_temp, good_normal) < 0.99998)
        {
            indeces.erase(indeces.begin() + i);
            if(i == indeces.size() - 1)
            break;
            --i;
            Nindeces = indeces.size();
        }
        }
        if(Nindeces < 3)
        return false;
        else
        return true;
    }


    void MakeRightHandFace(boost::container::small_vector<size_t, 8> &indeces, Vector3D const &point, vector<Vector3D> const &face_points,
                                                 std::array<size_t, 128> &temp, double areascale)
    {
        Vector3D V1, V2;
        size_t counter = 0;
        const size_t N = indeces.size();
        V1 = face_points[indeces[counter + 1]];
        V1 -= face_points[indeces[counter]];
        double AScale = 1e-14 * areascale;
        while (ScalarProd(V1, V1) < AScale)
        {
            ++counter;
            assert(counter < N);
            V1 = face_points[indeces[(counter + 1) % N]];
            V1 -= face_points[indeces[counter]];
        }
        V2 = face_points[indeces[(counter + 2) % N]];
        V2 -= face_points[indeces[(counter + 1) % N]];
        while (ScalarProd(V2, V2) < AScale)
        {
            ++counter;
            assert(counter < 2 * N);
            V2 = face_points[indeces[(counter + 2) % N]];
            V2 -= face_points[indeces[(counter + 1) % N]];
        }
        // Do we need to flip handness?
        if (ScalarProd(CrossProduct(V1, V2), point - face_points[indeces[0]]) > 0)
        {
            const size_t Ninner = indeces.size();
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
#endif
            for (size_t j = 0; j < Ninner; ++j)
                temp[j] = indeces[j];
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
#endif
            for (size_t i = 0; i < N; ++i)
                indeces[i] = temp[(N - i - 1)];
        }
    }

    size_t NextLoopTetra(Tetrahedron const &cur_tetra, size_t last_tetra, size_t N0, size_t N1)
    {
        size_t i = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (; i < 4; i++)
        {
            size_t point = cur_tetra.points[i];
            if (point != N0 && point != N1 && cur_tetra.neighbors[i] != last_tetra)
                break;
        }
        if(i >= 4)
            throw UniversalError("Bad NextLoopTetra");
        return cur_tetra.neighbors[i];
    }

#ifdef RICH_MPI
    std::pair<Vector3D, Vector3D> GetBoundingBox(Tessellation3D const &tproc, int rank)
    {
        vector<Vector3D> const &face_points = tproc.GetFacePoints();
        face_vec faces = tproc.GetCellFaces(static_cast<size_t>(rank));
        Vector3D ll = face_points[tproc.GetPointsInFace(faces[0])[0]];
        Vector3D ur(ll);
        size_t Nface = faces.size();
        for (size_t i = 0; i < Nface; ++i)
        {
            point_vec const &findex = tproc.GetPointsInFace(faces[i]);
            size_t Nindex = findex.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
            for (size_t j = 0; j < Nindex; j++)
            {
                ll.x = std::min(ll.x, face_points[findex[j]].x);
                ll.y = std::min(ll.y, face_points[findex[j]].y);
                ll.z = std::min(ll.z, face_points[findex[j]].z);
                ur.x = std::max(ur.x, face_points[findex[j]].x);
                ur.y = std::max(ur.y, face_points[findex[j]].y);
                ur.z = std::max(ur.z, face_points[findex[j]].z);
            }
        }
        return std::pair<Vector3D, Vector3D>(ll, ur);
    }

    void TalkSymmetry(vector<int> &to_talk_with)
    {
        int wsize;
        MPI_Comm_size(MPI_COMM_WORLD, &wsize);
        vector<int> totalk(static_cast<std::size_t>(wsize), 0);
        vector<int> scounts(totalk.size(), 1);
        for (std::size_t i = 0; i < to_talk_with.size(); ++i)
            totalk[to_talk_with[i]] = 1;
        int nrecv = 0;
        MPI_Reduce_scatter(&totalk[0], &nrecv, &scounts[0], MPI_INT, MPI_SUM,
                                             MPI_COMM_WORLD);

        vector<MPI_Request> req(to_talk_with.size());
        std::vector<int> dummy_send(to_talk_with.size());
        for (std::size_t i = 0; i < to_talk_with.size(); ++i)
            MPI_Isend(&dummy_send[i], 1, MPI_INT, to_talk_with[i], 3, MPI_COMM_WORLD, &req[i]);
        vector<int> talkwithme;
        for (int i = 0; i < nrecv; ++i)
        {
            MPI_Status status;
            MPI_Recv(&wsize, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
            talkwithme.push_back(status.MPI_SOURCE);
        }
        if (!to_talk_with.empty())
            MPI_Waitall(static_cast<int>(to_talk_with.size()), &req[0], MPI_STATUSES_IGNORE);
        vector<int> new_talk_with_me;
        for (std::size_t i = 0; i < to_talk_with.size(); ++i)
            if (std::find(talkwithme.begin(), talkwithme.end(), to_talk_with[i]) != talkwithme.end())
                new_talk_with_me.push_back(to_talk_with[i]);
        to_talk_with = new_talk_with_me;
    }
#endif //RICH_MPI

    void CalcFaceAreaCM(boost::container::small_vector<size_t, 8> const &indeces, std::vector<Vector3D> const &allpoints,
                                            std::array<Vector3D, 128> &points, double &Area, Vector3D &CM,
                                            std::array<double, 128> &Atemp)
    {
        //CM.Set(0.0, 0.0, 0.0);
        size_t Nloop = indeces.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (size_t i = 0; i < Nloop; i++)
            points[i] = allpoints[indeces[i]];
        Nloop -= 2;
        Area = 0;
        //Vector3D temp3, temp4, temp5;
        for (size_t i = 0; i < Nloop; i++)
        {
            //temp4.Set(points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z);
            Vector3D temp4(points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z);
            //temp5.Set(points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z);
            Vector3D temp5(points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z);
            Vector3D temp3;
            CrossProduct(temp4, temp5, temp3);
            Atemp[i] = 0.3333333333333333 * 0.5 * fastsqrt(ScalarProd(temp3, temp3));
        }
        double x = 0, y = 0, z = 0;
#ifdef __INTEL_COMPILER
#pragma vector aligned
        //#pragma omp simd reduction(+:x, y, z, Area)
#endif
        for (size_t i = 0; i < Nloop; i++)
        {
            double A = Atemp[i];
            x += A * points[0].x;
            y += A * points[0].y;
            z += A * points[0].z;
            x += A * points[i + 1].x;
            y += A * points[i + 1].y;
            z += A * points[i + 1].z;
            x += A * points[i + 2].x;
            y += A * points[i + 2].y;
            z += A * points[i + 2].z;
            Area += 3.0 * A;
        }
        CM.Set(x, y, z);
        CM *= (1.0 / (Area + DBL_MIN * 100)); //prevent overflow
    }

    bool PointInDomain(Vector3D const &ll, Vector3D const &ur, Vector3D const &point)
    {
        if (point.x > ll.x && point.x < ur.x && point.y > ll.y && point.y < ur.y && point.z > ll.z && point.z < ur.z)
            return true;
        else
            return false;
    }

    Vector3D MirrorPoint(Face const &face, Vector3D const &point)
    {
        Vector3D normal = CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
        normal = normal / abs(normal);
        return point - (2 * ScalarProd(point - face.vertices[0], normal)) * normal;
    }
}

#ifdef RICH_MPI
vector<Vector3D> Voronoi3D::UpdateMPIPoints(Tessellation3D const &vproc, int rank,
                                                                                        vector<Vector3D> const &points, vector<std::size_t> &selfindex, vector<int> &sentproc,
                                                                                        vector<vector<std::size_t>> &sentpoints)
{
    vector<Vector3D> res;
    res.reserve(points.size());
    selfindex.clear();
    std::size_t npoints = points.size();
    std::size_t nproc = vproc.GetPointNo();
    vector<std::size_t> neighbors = vproc.GetNeighbors(static_cast<std::size_t>(rank));
    vector<std::size_t> realneigh;
    sentpoints.clear();
    sentproc.clear();
    for (std::size_t i = 0; i < neighbors.size(); ++i)
        if (static_cast<std::size_t>(neighbors[i]) < nproc)
        {
            realneigh.push_back(neighbors[i]);
            sentproc.push_back(static_cast<int>(neighbors[i]));
        }
    std::size_t Nreal = realneigh.size();
    sentpoints.resize(sentproc.size());

    for (std::size_t i = 0; i < npoints; ++i)
    {
        Vector3D temp = points[i];
        if (PointInPoly(vproc, temp, static_cast<std::size_t>(rank)))
        {
            res.push_back(temp);
            selfindex.push_back(i);
            continue;
        }
        bool good = false;
        for (std::size_t j = 0; j < Nreal; ++j)
        {
            if (PointInPoly(vproc, temp, realneigh[j]))
            {
                sentpoints[j].push_back(i);
                good = true;
                break;
            }
        }
        if (good)
            continue;
        for (std::size_t j = 0; j < nproc; ++j)
        {
            if (std::find(realneigh.begin(), realneigh.end(), j) != realneigh.end() || j == static_cast<std::size_t>(rank))
                continue;
            if (PointInPoly(vproc, temp, j))
            {
                good = true;
                std::size_t index = std::find(sentproc.begin(), sentproc.end(), j) - sentproc.begin();
                if (index >= sentproc.size())
                {
                    sentproc.push_back(static_cast<int>(j));
                    sentpoints.push_back(vector<std::size_t>(1, i));
                }
                else
                    sentpoints[index].push_back(i);
                break;
            }
        }
        if (good)
            continue;
        UniversalError eo("Point is not inside any processor");
        eo.addEntry("CPU rank", rank);
        eo.addEntry("Point number", static_cast<double>(i));
        eo.addEntry("Point x cor", points[i].x);
        eo.addEntry("Point y cor", points[i].y);
        eo.addEntry("Point z cor", points[i].z);
        vproc.output("vproc_" + int2str(rank) + ".bin");
        std::array<Vector3D, 4> vec;
        face_vec faces_error = vproc.GetCellFaces(static_cast<size_t>(rank));
        for (size_t j = 0; j < faces_error.size(); ++j)
        {
            point_vec_v f_points = VectorValues(vproc.GetFacePoints(), vproc.GetPointsInFace(faces_error[j]));
            for (size_t k = 0; k < f_points.size(); ++k)
            {
                std::cout << "Rank " << rank << " face " << faces_error[j] << " point " << k << " cor " << f_points[k].x
                                    << " " << f_points[k].y << " " << f_points[k].z << std::endl;
            }
            vec[0] = f_points[0];
            vec[1] = f_points[1];
            vec[2] = f_points[2];
            vec[3] = vproc.GetMeshPoint(rank);
            double s1 = orient3d(vec);
            vec[3] = points[i];
            double s2 = orient3d(vec);
            std::cout << "s1 = " << s1 << " s2 = " << s2 << std::endl;
        }
        for (std::size_t l = 0; l < Nreal; ++l)
        {
            faces_error = vproc.GetCellFaces(static_cast<size_t>(realneigh[l]));
            for (size_t j = 0; j < faces_error.size(); ++j)
            {
                point_vec_v f_points = VectorValues(vproc.GetFacePoints(), vproc.GetPointsInFace(faces_error[j]));
                for (size_t k = 0; k < f_points.size(); ++k)
                {
                    std::cout << "Rank " << realneigh[l] << " face " << faces_error[j] << " point " << k << " cor " << f_points[k].x
                                        << " " << f_points[k].y << " " << f_points[k].z << std::endl;
                }
                vec[0] = f_points[0];
                vec[1] = f_points[1];
                vec[2] = f_points[2];
                vec[3] = vproc.GetMeshPoint(rank);
                double s1 = orient3d(vec);
                vec[3] = points[i];
                double s2 = orient3d(vec);
                std::cout << "s1 = " << s1 << " s2 = " << s2 << std::endl;
            }
        }
        throw eo;
    }
    // Send/Recv the points
    // Communication
    int wsize;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    vector<int> totalk(static_cast<std::size_t>(wsize), 0);
    vector<int> scounts(totalk.size(), 1);
    for (std::size_t i = 0; i < sentproc.size(); ++i)
        totalk[sentproc[i]] = 1;
    int nrecv;
    MPI_Reduce_scatter(&totalk[0], &nrecv, &scounts[0], MPI_INT, MPI_SUM,
                                         MPI_COMM_WORLD);

    vector<MPI_Request> req(sentproc.size());
    std::vector<int> dummy_send(req.size());
    for (std::size_t i = 0; i < sentproc.size(); ++i)
        MPI_Isend(&dummy_send[i], 1, MPI_INT, sentproc[i], 3, MPI_COMM_WORLD, &req[i]);
    vector<int> talkwithme;
    for (int i = 0; i < nrecv; ++i)
    {
        MPI_Status status;
        MPI_Recv(&wsize, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
        talkwithme.push_back(status.MPI_SOURCE);
    }
    MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    for (std::size_t i = 0; i < talkwithme.size(); ++i)
    {
        if (std::find(sentproc.begin(), sentproc.end(), talkwithme[i]) == sentproc.end())
        {
            sentproc.push_back(talkwithme[i]);
            sentpoints.push_back(vector<std::size_t>());
        }
    }
    // Point exchange
    vector<vector<Vector3D>> incoming;
    if (points.empty())
    {
        vector<Vector3D> dummy(1);
        incoming = MPI_exchange_data(sentproc, sentpoints, dummy);
    }
    else
        incoming = MPI_exchange_data(sentproc, sentpoints, points);
    // Combine the vectors
    for (std::size_t i = 0; i < incoming.size(); ++i)
        for (std::size_t j = 0; j < incoming[i].size(); ++j)
            res.push_back(incoming[i][j]);
    return res;
}
#endif //RICH_MPI

Voronoi3D::Voronoi3D() : ll_(Vector3D()), ur_(Vector3D()), Norg_(0), bigtet_(0), set_temp_(std::set<int>()), stack_temp_(std::stack<int>()),
                        del_(Delaunay3D()), PointTetras_(vector<tetra_vec>()), R_(vector<double>()), tetra_centers_(vector<Vector3D>()),
                        FacesInCell_(vector<face_vec>()),
                        PointsInFace_(vector<point_vec>()),
                        FaceNeighbors_(vector<std::pair<std::size_t, std::size_t>>()),
                        CM_(vector<Vector3D>()), Face_CM_(vector<Vector3D>()),
                        volume_(vector<double>()), area_(vector<double>()), duplicated_points_(vector<vector<std::size_t>>()),
                        sentprocs_(vector<int>()), duplicatedprocs_(vector<int>()), sentpoints_(vector<vector<std::size_t>>()), Nghost_(vector<vector<std::size_t>>()),
                        self_index_(vector<std::size_t>()), temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>())
                        #ifdef RICH_MPI
<<<<<<< HEAD
                        , hilbertAgent(HilbertAgent(this->ll_, this->ur_, DEFAULT_HILBERT_ACCURACY)), initialRadius(0.0), firstCall(true)
=======
                        ,hilbertAgent(HilbertAgent(this->ll_, this->ur_, DEFAULT_HILBERT_ACCURACY)), initialRadius(0.0)
>>>>>>> 342f02048afef504c85b1f72ffc45f1cefee8c24
                        #endif // RICH_MPI
{
}

Voronoi3D::Voronoi3D(std::vector<Face> const& box_faces) : Norg_(0), bigtet_(0), set_temp_(std::set<int>()), stack_temp_(std::stack<int>()),
<<<<<<< HEAD
                                                        del_(Delaunay3D()), PointTetras_(vector<tetra_vec>()), R_(vector<double>()), tetra_centers_(vector<Vector3D>()),
                                                        FacesInCell_(vector<face_vec>()),
                                                        PointsInFace_(vector<point_vec>()),
                                                        FaceNeighbors_(vector<std::pair<std::size_t, std::size_t>>()),
                                                        CM_(vector<Vector3D>()), Face_CM_(vector<Vector3D>()),
                                                        volume_(vector<double>()), area_(vector<double>()), duplicated_points_(vector<vector<std::size_t>>()),
                                                        sentprocs_(vector<int>()), duplicatedprocs_(vector<int>()), sentpoints_(vector<vector<std::size_t>>()), Nghost_(vector<vector<std::size_t>>()),
                                                        self_index_(vector<std::size_t>()), temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>()), box_faces_(box_faces)
                                                        #ifdef RICH_MPI
                                                        , hilbertAgent(HilbertAgent(this->ll_, this->ur_, DEFAULT_HILBERT_ACCURACY)), initialRadius(0.0), firstCall(true)
                                                        #endif // RICH_MPI
=======
                                                          del_(Delaunay3D()), PointTetras_(vector<tetra_vec>()), R_(vector<double>()), tetra_centers_(vector<Vector3D>()),
                                                          FacesInCell_(vector<face_vec>()),
                                                          PointsInFace_(vector<point_vec>()),
                                                          FaceNeighbors_(vector<std::pair<std::size_t, std::size_t>>()),
                                                          CM_(vector<Vector3D>()), Face_CM_(vector<Vector3D>()),
                                                          volume_(vector<double>()), area_(vector<double>()), duplicated_points_(vector<vector<std::size_t>>()),
                                                          sentprocs_(vector<int>()), duplicatedprocs_(vector<int>()), sentpoints_(vector<vector<std::size_t>>()), Nghost_(vector<vector<std::size_t>>()),
                                                          self_index_(vector<std::size_t>()), temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>()), box_faces_(box_faces)
                                                          #ifdef RICH_MPI
                                                          ,hilbertAgent(HilbertAgent(this->ll_, this->ur_, DEFAULT_HILBERT_ACCURACY)), initialRadius(0.0)
                                                          #endif // RICH_MPI
>>>>>>> 342f02048afef504c85b1f72ffc45f1cefee8c24
{
    size_t const Nfaces = box_faces.size();
    if(Nfaces < 4)
        throw UniversalError("Zero face vector in Voronoi3D constructor");
    ll_ = box_faces[0].vertices[0];
    ur_ = ll_;
    for(size_t i = 0; i < Nfaces; ++i)
    {
        size_t const Nvertices = box_faces[i].vertices.size();
        for(size_t j = 0; j < Nvertices; ++j)
        {
            ll_.x = std::min(ll_.x, box_faces[i].vertices[j].x);
            ll_.y = std::min(ll_.y, box_faces[i].vertices[j].y);
            ll_.z = std::min(ll_.z, box_faces[i].vertices[j].z);
            ur_.x = std::max(ur_.x, box_faces[i].vertices[j].x);
            ur_.y = std::max(ur_.y, box_faces[i].vertices[j].y);
            ur_.z = std::max(ur_.z, box_faces[i].vertices[j].z);
        }
    }
}

Voronoi3D::Voronoi3D(Vector3D const &ll, Vector3D const &ur) : ll_(ll), ur_(ur), Norg_(0), bigtet_(0), set_temp_(std::set<int>()), stack_temp_(std::stack<int>()),
                                                              del_(Delaunay3D()), PointTetras_(vector<tetra_vec>()), R_(vector<double>()), tetra_centers_(vector<Vector3D>()),
                                                              FacesInCell_(vector<face_vec>()),
                                                              PointsInFace_(vector<point_vec>()),
                                                              FaceNeighbors_(vector<std::pair<std::size_t, std::size_t>>()),
                                                              CM_(vector<Vector3D>()), Face_CM_(vector<Vector3D>()),
                                                              volume_(vector<double>()), area_(vector<double>()), duplicated_points_(vector<vector<std::size_t>>()),
                                                              sentprocs_(vector<int>()), duplicatedprocs_(vector<int>()), sentpoints_(vector<vector<std::size_t>>()), Nghost_(vector<vector<std::size_t>>()),
                                                              self_index_(vector<std::size_t>()), temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>()), box_faces_(std::vector<Face> ())
                                                              #ifdef RICH_MPI
<<<<<<< HEAD
                                                              , hilbertAgent(HilbertAgent(this->ll_, this->ur_, DEFAULT_HILBERT_ACCURACY)), initialRadius(0.0), firstCall(true)
=======
                                                              ,hilbertAgent(HilbertAgent(this->ll_, this->ur_, DEFAULT_HILBERT_ACCURACY))
>>>>>>> 342f02048afef504c85b1f72ffc45f1cefee8c24
                                                              #endif // RICH_MPI
                                                              {}

void Voronoi3D::CalcRigidCM(std::size_t face_index)
{
    Vector3D normal = normalize(del_.points_[FaceNeighbors_[face_index].first] - del_.points_[FaceNeighbors_[face_index].second]);
    std::size_t real, other;
    if (FaceNeighbors_[face_index].first >= Norg_)
    {
        real = FaceNeighbors_[face_index].second;
        other = FaceNeighbors_[face_index].first;
    }
    else
    {
        real = FaceNeighbors_[face_index].first;
        other = FaceNeighbors_[face_index].second;
    }
    CM_[other] = CM_[real] - 2 * normal * ScalarProd(normal, CM_[real] - tetra_centers_[PointsInFace_[face_index][0]]);
}

vector<Vector3D> Voronoi3D::CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t>> const &to_duplicate,
                                                 vector<vector<size_t>> &past_duplicate)
{
    size_t Ncheck = to_duplicate.size();
    vector<std::pair<std::size_t, std::size_t>> to_add;
    to_add.reserve(Ncheck);
    vector<Face> faces = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<Vector3D> res;
    bool first_time = past_duplicate.empty();
    if (first_time)
        past_duplicate.resize(faces.size());
    for (std::size_t i = 0; i < Ncheck; ++i)
    {
        if (first_time || !std::binary_search(past_duplicate[to_duplicate[i].first].begin(),
                                              past_duplicate[to_duplicate[i].first].end(), to_duplicate[i].second))
        {
            res.push_back(MirrorPoint(faces[to_duplicate[i].first], del_.points_[to_duplicate[i].second]));
            to_add.push_back(to_duplicate[i]);
        }
    }
    for (size_t i = 0; i < to_add.size(); ++i)
        past_duplicate[to_add[i].first].push_back(to_add[i].second);
    for (size_t i = 0; i < past_duplicate.size(); ++i)
        std::sort(past_duplicate[i].begin(), past_duplicate[i].end());
    return res;
}

#ifdef RICH_MPI
vector<Vector3D> Voronoi3D::CreateBoundaryPointsMPI(vector<std::pair<std::size_t, std::size_t>> const &to_duplicate,
                                                    Tessellation3D const &tproc, vector<vector<size_t>> &self_duplicate)
{
    vector<vector<size_t>> to_send;

    vector<Face> box_faces = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<Vector3D> box_normals = GetBoxNormals(ll_, ur_, box_faces_);
    vector<vector<size_t>> box_candidates(box_normals.size());
    self_duplicate.resize(box_faces.size());

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::size_t Nproc = tproc.GetPointNo();
    size_t Ndup = to_duplicate.size();
    for (std::size_t i = 0; i < Ndup; i++)
    {
        std::pair<std::size_t, std::size_t> const &neigh = tproc.GetFaceNeighbors(to_duplicate[i].first);
        if (neigh.first < Nproc && static_cast<int>(neigh.first) != rank)
        {
            if (std::find(duplicatedprocs_.begin(), duplicatedprocs_.end(), neigh.first) == duplicatedprocs_.end())
                duplicatedprocs_.push_back(static_cast<int>(neigh.first));
        }
        if (neigh.second < Nproc && static_cast<int>(neigh.second) != rank)
        {
            if (std::find(duplicatedprocs_.begin(), duplicatedprocs_.end(), neigh.second) == duplicatedprocs_.end())
                duplicatedprocs_.push_back(static_cast<int>(neigh.second));
        }
    }
    TalkSymmetry(duplicatedprocs_);
    to_send.resize(duplicatedprocs_.size());
    duplicated_points_.resize(duplicatedprocs_.size());
    vector<Vector3D> res;
    // Get the indeces and deal with selfboundary
    for (std::size_t i = 0; i < Ndup; ++i)
    {
        std::pair<std::size_t, std::size_t> const &neigh = tproc.GetFaceNeighbors(to_duplicate[i].first);
        if (neigh.first != static_cast<std::size_t>(rank))
        {
            if (neigh.first < Nproc)
            {
                std::size_t index = std::find(duplicatedprocs_.begin(), duplicatedprocs_.end(), static_cast<int>(neigh.first)) - duplicatedprocs_.begin();
                if (index < duplicatedprocs_.size())
                    to_send[index].push_back(to_duplicate[i].second);
            }
            else
            {
                Vector3D face_normal = tproc.GetMeshPoint(neigh.first) - tproc.GetMeshPoint(neigh.second);
                face_normal *= 1.0 / abs(face_normal);
                size_t index = BoxIndex(box_normals, face_normal);
                box_candidates[index].push_back(to_duplicate[i].second);
            }
        }
        if (neigh.second != static_cast<std::size_t>(rank))
        {
            if (neigh.second < Nproc)
            {
                std::size_t index = std::find(duplicatedprocs_.begin(), duplicatedprocs_.end(), static_cast<int>(neigh.second)) - duplicatedprocs_.begin();
                if (index < duplicatedprocs_.size())
                    to_send[index].push_back(to_duplicate[i].second);
            }
            else
            {
                Vector3D face_normal = tproc.GetMeshPoint(neigh.second) - tproc.GetMeshPoint(neigh.first);
                face_normal *= 1.0 / abs(face_normal);
                size_t index = BoxIndex(box_normals, face_normal);
                box_candidates[index].push_back(to_duplicate[i].second);
            }
        }
    }
    // Clean
    for (size_t i = 0; i < duplicated_points_.size(); ++i)
    {
        std::sort(to_send[i].begin(), to_send[i].end());
        to_send[i] = unique(to_send[i]);
        vector<size_t> temp;
        size_t Nsend = to_send[i].size();
        vector<size_t> indeces = sort_index(duplicated_points_[i]);
        sort(duplicated_points_[i].begin(), duplicated_points_[i].end());

        for (size_t j = 0; j < Nsend; ++j)
        {
            if (duplicated_points_[i].empty() ||
                    !std::binary_search(duplicated_points_[i].begin(), duplicated_points_[i].end(), to_send[i][j]))
                temp.push_back(to_send[i][j]);
        }
        to_send[i] = temp;
        duplicated_points_[i].insert(duplicated_points_[i].end(), temp.begin(), temp.end());
        Nsend = indeces.size();
        temp = duplicated_points_[i];
        for (size_t j = 0; j < Nsend; ++j)
            duplicated_points_[i][indeces[j]] = temp[j];
    }
    for (size_t i = 0; i < box_candidates.size(); ++i)
    {
        std::sort(box_candidates[i].begin(), box_candidates[i].end());
        box_candidates[i] = unique(box_candidates[i]);
        vector<size_t> indeces = sort_index(self_duplicate[i]);
        sort(self_duplicate[i].begin(), self_duplicate[i].end());

        vector<size_t> temp;

        for (size_t j = 0; j < box_candidates[i].size(); ++j)
        {
            if (self_duplicate[i].empty() ||
                    !std::binary_search(self_duplicate[i].begin(), self_duplicate[i].end(), box_candidates[i][j]))
            {
                temp.push_back(box_candidates[i][j]);
                res.push_back(MirrorPoint(box_faces[i], del_.points_[temp.back()]));
            }
        }
        self_duplicate[i].insert(self_duplicate[i].end(), temp.begin(), temp.end());
        box_candidates[i] = temp;
        size_t Nsend = indeces.size();
        temp = self_duplicate[i];
        for (size_t j = 0; j < Nsend; ++j)
            self_duplicate[i][indeces[j]] = temp[j];
    }
    // Communicate
    vector<vector<Vector3D>> toadd = MPI_exchange_data(duplicatedprocs_, to_send, del_.points_);
    // Add points
    Nghost_.resize(toadd.size());
    size_t temp_add = del_.points_.size();
    for (std::size_t i = 0; i < toadd.size(); ++i)
        for (std::size_t j = 0; j < toadd[i].size(); ++j)
        {
            Nghost_[i].push_back(temp_add + res.size());
            res.push_back(toadd[i][j]);
        }
    return res;
}
#endif //RICH_MPI

vector<vector<std::size_t>> const &Voronoi3D::GetGhostIndeces(void) const
{
    return Nghost_;
}

#ifdef RICH_MPI
void Voronoi3D::Build(vector<Vector3D> const &points, Tessellation3D const &tproc)
{
  //assert(points.size() > 0);
  // Clear data
  PointTetras_.clear();
  R_.clear();
  R_.reserve(points.size() * 11);
  tetra_centers_.clear();
  tetra_centers_.reserve(points.size() * 11);
  del_.Clean();
  // Voronoi Data
  FacesInCell_.clear();
  PointsInFace_.clear();
  FaceNeighbors_.clear();
  CM_.clear();
  Face_CM_.clear();
  volume_.clear();
  area_.clear();
  Nghost_.clear();
  duplicatedprocs_.clear();
  duplicated_points_.clear();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vector<Vector3D> new_points = UpdateMPIPoints(tproc, rank, points, self_index_, sentprocs_, sentpoints_);
  Norg_ = new_points.size();
  /*	if (Norg_ == 0)
	{
	std::cout << "Zero Norg in rank " << rank << std::endl;
	std::cout << "Rank CM " << tproc.GetCellCM(static_cast<size_t>(rank)).x << ","
	<< tproc.GetCellCM(static_cast<size_t>(rank)).y << "," << tproc.GetCellCM(static_cast<size_t>(rank)).z << std::endl;
	std::cout << "Rank point " << tproc.GetMeshPoint(static_cast<size_t>(rank)).x << ","
	<< tproc.GetMeshPoint(static_cast<size_t>(rank)).y << "," << tproc.GetMeshPoint(static_cast<size_t>(rank)).z << std::endl;
	std::cout << "Rank R " << tproc.GetWidth(static_cast<size_t>(rank)) << std::endl;
	}
	assert(Norg_ > 0);*/
  std::pair<Vector3D, Vector3D> bounding_box = GetBoundingBox(tproc, rank);

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = MPI_Wtime();
#endif
  std::vector<size_t> order = HilbertOrder3D(new_points);

#ifdef vdebug
  std::vector<Vector3D> bbox;
  bbox.push_back(bounding_box.first);
  bbox.push_back(bounding_box.second);
  write_vecst(order, "order_" + int2str(rank) + ".bin");
  write_vec3d(new_points, "points0_" + int2str(rank) + ".bin");
  write_vec3d(bbox, "bb_" + int2str(rank) + ".bin");
#endif

  del_.Build(new_points, bounding_box.second, bounding_box.first, order);

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  if (rank == 0)
    std::cout << "First build time " << t1 - t0 << std::endl;
#endif

  R_.resize(del_.tetras_.size());
  std::fill(R_.begin(), R_.end(), -1);
  tetra_centers_.resize(R_.size());
  bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

  vector<vector<size_t>> self_duplicate;
  vector<std::pair<std::size_t, std::size_t>> ghost_index;
  MPIFirstIntersections(tproc, ghost_index);

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  if (rank == 0)
    std::cout << "First ghost time " << t0 - t1 << std::endl;
#endif

  vector<Vector3D> extra_points = CreateBoundaryPointsMPI(ghost_index, tproc, self_duplicate);

#ifdef vdebug
  write_vec3d(extra_points, "points1_" + int2str(rank) + ".bin");
#endif

  try
  {
    del_.BuildExtra(extra_points);
  }
  catch (UniversalError &eo)
  {
    std::cout << "Error in first extra rank " << rank << std::endl;
    string fname("extra_" + int2str(rank) + ".bin");
    output_buildextra(fname);
    tproc.output("vproc_" + int2str(rank) + ".bin");
    throw;
  }

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (rank == 0)
    std::cout << "First ghost build time " << t1 - t0 << std::endl;
#endif

  R_.resize(del_.tetras_.size());
  std::fill(R_.begin(), R_.end(), -1);
  tetra_centers_.resize(R_.size());
  bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);
  vector<unsigned char> checked_clear(Norg_, 0);
  ghost_index = FindIntersections(tproc, 1, checked_clear); // intersecting tproc face, point index

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  if (rank == 0)
    std::cout << "Second ghost time " << t0 - t1 << std::endl;
#endif

  extra_points = CreateBoundaryPointsMPI(ghost_index, tproc, self_duplicate);

#ifdef vdebug
  write_vec3d(extra_points, "points2_" + int2str(rank) + ".bin");
#endif

  try
  {
    del_.BuildExtra(extra_points);
  }
  catch (UniversalError &eo)
  {
    std::cout << "Error in second extra rank " << rank << std::endl;
    string fname("extra_" + int2str(rank) + ".bin");
    output_buildextra(fname);
    tproc.output("vproc_" + int2str(rank) + ".bin");
    throw;
  }

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (rank == 0)
    std::cout << "Second ghost build time " << t1 - t0 << std::endl;
#endif

  R_.resize(del_.tetras_.size());
  std::fill(R_.begin(), R_.end(), -1);
  tetra_centers_.resize(R_.size());
  bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

  ghost_index = FindIntersections(tproc, 2, checked_clear);

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  if (rank == 0)
    std::cout << "Third ghost time " << t0 - t1 << std::endl;
#endif

  extra_points = CreateBoundaryPointsMPI(ghost_index, tproc, self_duplicate);

#ifdef vdebug
  write_vec3d(extra_points, "points3_" + int2str(rank) + ".bin");
#endif

  try
  {
    del_.BuildExtra(extra_points);
  }
  catch (UniversalError &eo)
  {
    std::cout << "Error in third extra rank " << rank << std::endl;
    string fname("extra_" + int2str(rank) + ".bin");
    output_buildextra(fname);
    tproc.output("vproc_" + int2str(rank) + ".bin");
    throw;
  }

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (rank == 0)
    std::cout << "Third ghost build time " << t1 - t0 << std::endl;
#endif

  R_.resize(del_.tetras_.size());
  std::fill(R_.begin(), R_.end(), -1);
  tetra_centers_.resize(R_.size());
  bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

  ghost_index = FindIntersections(tproc, 3, checked_clear);

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  if (rank == 0)
    std::cout << "Fourth ghost time " << t0 - t1 << std::endl;
#endif

  extra_points = CreateBoundaryPointsMPI(ghost_index, tproc, self_duplicate);

#ifdef vdebug
  write_vec3d(extra_points, "points4_" + int2str(rank) + ".bin");
#endif

  try
  {
    del_.BuildExtra(extra_points);
  }
  catch (UniversalError &eo)
  {
    std::cout << "Error in fourth extra rank " << rank << std::endl;
    string fname("extra_" + int2str(rank) + ".bin");
    output_buildextra(fname);
    tproc.output("vproc_" + int2str(rank) + ".bin");
    throw;
  }

#ifdef timing
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (rank == 0)
    std::cout << "Fourth ghost build time " << t1 - t0 << std::endl;
#endif
  bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);
  std::vector<std::pair<size_t, size_t>>().swap(ghost_index);
  std::vector<Vector3D>().swap(extra_points);

  R_.resize(del_.tetras_.size());
  std::fill(R_.begin(), R_.end(), -1);
  tetra_centers_.resize(R_.size());

  CM_.resize(del_.points_.size());
  volume_.resize(Norg_);

  // Create Voronoi
  BuildVoronoi(order);

  std::vector<double>().swap(R_);
  std::vector<tetra_vec>().swap(PointTetras_);

  CalcAllCM();
  for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
    if (BoundaryFace(i))
      CalcRigidCM(i);
  // communicate the ghost CM
  vector<vector<Vector3D>> incoming = MPI_exchange_data(duplicatedprocs_, duplicated_points_, CM_);
  // Add the recieved CM
  for (size_t i = 0; i < incoming.size(); ++i)
    for (size_t j = 0; j < incoming.at(i).size(); ++j)
      CM_[Nghost_.at(i).at(j)] = incoming[i][j];
}
#endif

/**
 * gets a point index, and returns the maximal radius of the tetrahedra containing that point.
 * @param index the index of the point (within the points list)
*/
double Voronoi3D::GetMaxRadius(std::size_t index)
{
    std::size_t N = PointTetras_[index].size();
    double res = 0;
    #ifdef __INTEL_COMPILER
    #pragma ivdep
    #endif
    for(std::size_t i = 0; i < N; ++i)
    {
        res = std::max(res, GetRadius(PointTetras_[index][i]));
    }
    return res;
}

/**
 * if the initial box does not exist, builds its faces according to the leftmost and rightmost points.
 * If it does, does not build the faces again.
 * @return the normals to the faces
*/
void Voronoi3D::InitialBoxBuild(std::vector<Face> &box, std::vector<Vector3D> &normals)
{
    box = box_faces_.empty() ? BuildBox(this->ll_, this->ur_) : this->box_faces_;
    size_t Nfaces = box.size();
    normals.resize(Nfaces);

    // calculates the normals for each one of the box's faces
    for (size_t i = 0; i < Nfaces; ++i)
    {
        normals[i] = CrossProduct(box[i].vertices[1] - box[i].vertices[0], box[i].vertices[2] - box[i].vertices[0]);
        normals[i] *= (1.0 / fastsqrt(ScalarProd(normals[i], normals[i])));
    }
}

#ifdef RICH_MPI
namespace
{
    /**
     * returns the number of new <finish> messages to arrive.
    */
    int getNewFinished()
    {
        MPI_Status status;
        int newFinished = 0;
        int receivedFinished = 0;

        MPI_Iprobe(MPI_ANY_SOURCE, RICH_TESELLATION_FINISHED_TAG, MPI_COMM_WORLD, &receivedFinished, MPI_STATUS_IGNORE);
        while(receivedFinished)
        {
            int dummy;
            MPI_Recv(&dummy, 1, MPI_BYTE, MPI_ANY_SOURCE, RICH_TESELLATION_FINISHED_TAG, MPI_COMM_WORLD, &status);
            newFinished++;
            MPI_Iprobe(MPI_ANY_SOURCE, RICH_TESELLATION_FINISHED_TAG, MPI_COMM_WORLD, &receivedFinished, MPI_STATUS_IGNORE);
        }
        return newFinished;
    }

    void sendFinished()
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        int dummy = 0;

        for(int _rank = 0; _rank < size; _rank++)
        {
            MPI_Send(&dummy, 1, MPI_BYTE, _rank, RICH_TESELLATION_FINISHED_TAG, MPI_COMM_WORLD);
        }
    }

    template<typename T>
    void reportDuplications(const std::vector<T> &vector)
    {
        for(size_t i = 0; i < vector.size(); i++)
        {
            for(size_t j = 0; j < vector.size(); j++)
            {
                if(i == j) continue;
                if(vector[i] == vector[j])
                {
                    std::cout << "duplication found in indices " << i << " and " << j << ": " << vector[i] << std::endl;
                }
            }
        }
    }
}

void Voronoi3D::BuildInitialize(size_t num_points)
{
    // assert(num_points > 0);
    // Clear data
    PointTetras_.clear();
    R_.clear();
    if(num_points > 0) R_.reserve(num_points * 11);
    tetra_centers_.clear();
     if(num_points > 0) tetra_centers_.reserve(num_points * 11);
    // Voronoi Data
    del_.Clean();
    FacesInCell_.clear();
    PointsInFace_.clear();
    FaceNeighbors_.clear();
    CM_.clear();
    Face_CM_.clear();
    volume_.clear();
    area_.clear();
    Norg_ = num_points;
    duplicatedprocs_.clear();
    duplicated_points_.clear();
    self_index_.clear();
    sentprocs_.clear();
    sentpoints_.clear();
    Nghost_.clear();
}

void Voronoi3D::CheckToMirror(const Vector3D &point, double radius, std::vector<Face> &box, std::vector<Vector3D> &normals, std::vector<Vector3D> &points)
{
    for(size_t i = 0; i < box.size(); i++)
    {
        // check for intersecting the sphere with radius `radius` around `point`, with the `i`th face of `box`
        Sphere sphere(point, radius);
        if(FaceSphereIntersections(box[i], sphere, normals[i]))
        {
            // intersects! mirror the point
            points.push_back(MirrorPoint(box[i], point));
        }
    }
}

bool Voronoi3D::CheckForRebalance(const std::vector<Vector3D> &points) const
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    size_t mySize = points.size();
    size_t N;
    MPI_Allreduce(&mySize, &N, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    size_t ideal = N / size;
    int I_say = (mySize >= (BALANCE_FACTOR * static_cast<double>(ideal)))? 1 : 0; // if I say 'rebalance' or not
    int rebalance = 0;
    MPI_Allreduce(&I_say, &rebalance, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return (rebalance > 0);
}

/**
 * \author Maor Mizrachi
 * \brief The algorithm follows arepro paper (https://www.mpa-garching.mpg.de/~volker/arepo/arepo_paper.pdf), section 2.4.
*/
void Voronoi3D::PrepareToBuildHilbert(const std::vector<Vector3D> &points)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<Face> box;
    std::vector<Vector3D> normals;
    this->InitialBoxBuild(box, normals);

    bool sent_finished = false; // if I sent a finished message
    int finished = 0; // the number of finished ranks
    std::deque<int> current;    // the indexes of the current "bad" points (points with h_i <= s_i)
    for(size_t i = 0; i < points.size(); i++)
    {
        current.push_back(i);
    }

    //BruteForceFinder rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    //RangeTreeFinder rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    //SmartBruteForceFinder rangeFinder(&hilbertAgent, this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    //HashBruteForceFinder rangeFinder(&hilbertAgent, this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    OctTreeFinder rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_, this->ll_, this->ur_);
    //GroupRangeTreeFinder<256> rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    
    OctTree<Vector3D> pointsTree(this->ll_, this->ur_, points);
    RangeAgent rangeAgent(this->hilbertAgent, &rangeFinder);
    
    rangeAgent.buildHilbertTree(&pointsTree);

    std::vector<Vector3D> allMirrored;

    while(finished != size)
    {
        std::queue<RangeQueryData> queries;
        std::vector<Vector3D> mirroedPoints;

        for(std::size_t i = 0; i < current.size(); i++)
        {
            this->CheckToMirror(this->del_.points_[current[i]], this->radiuses[current[i]], box, normals, mirroedPoints);
            Vector3D &point = this->del_.points_[current[i]];
            queries.push({{point.x, point.y, point.z}, this->radiuses[current[i]]});
        }
    
        if(current.empty() and !sent_finished)
        {
            sendFinished();
            sent_finished = true;
        }

        QueryBatchInfo batchInfo = rangeAgent.runBatch(queries);

        finished += getNewFinished();

        std::vector<Vector3D> &newPoints = batchInfo.newPoints;

        for(int _rank = 0; _rank < size; _rank++)
        {
          std::vector<size_t> &receivedFromRank = rangeAgent.getRecvPoints()[_rank];
          if(receivedFromRank.empty())
          {
            continue;
          }
          size_t rankIdx = std::find(this->duplicatedprocs_.begin(), this->duplicatedprocs_.end(), _rank) - this->duplicatedprocs_.begin();
          if(rankIdx == this->duplicatedprocs_.size())
          {
            // new rank in this->duplicatedprocs_, initialize it
            this->duplicatedprocs_.push_back(_rank);
            this->duplicated_points_.emplace_back(std::vector<size_t>());
            this->Nghost_.emplace_back(std::vector<size_t>());
          }
          for(const size_t &RelativePointIdx : receivedFromRank)
          {
            // batchInfo.pointsFromRanks[_rank][i] holds an index of point, but this point will be added to my delaunay, so
            // its index there will be this->del_.points_.size() + batchInfo.pointsFromRanks[_rank][i]
            this->Nghost_[rankIdx].push_back(this->del_.points_.size() + RelativePointIdx);
          }
        }

        for(const Vector3D &point : mirroedPoints)
        {
            if(std::find(allMirrored.begin(), allMirrored.end(), point) == allMirrored.end())
            {
                allMirrored.push_back(point);
                newPoints.push_back(point);
            }
        }

        // performs internal tesselation:
        this->del_.BuildExtra(newPoints);

        // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
        this->R_.resize(this->del_.tetras_.size());
        std::fill(this->R_.begin(), this->R_.end(), -1);
        this->tetra_centers_.resize(this->R_.size());
        this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);

        for(std::size_t i = 0; i < current.size(); i++)
        {
            double s = 2 * this->GetMaxRadius(current[i]);

            if(this->radiuses[current[i]] >= s)
            {
                // h_i > s_i, remove point from current
                current.erase(current.cbegin() + i);
            }
            else
            {
                double radius = std::min(this->radiuses[current[i]] * 1.618, s); // todo define 1.618 as a constant
                this->radiuses[current[i]] = radius;
            }
        }
    }

    const std::vector<std::vector<size_t>> &sentPoints = rangeAgent.getSentPoints();
    const std::vector<int> &sentProc = rangeAgent.getSentProc();

    for(size_t i = 0; i < sentProc.size(); i++)
    {
      int _rank = sentProc[i];
      size_t rankIdx = std::find(this->duplicatedprocs_.begin(), this->duplicatedprocs_.end(), _rank) - this->duplicatedprocs_.begin();
      if(rankIdx == this->duplicatedprocs_.size())
      {
        // new rank in this->duplicatedprocs_, initialize it
        this->duplicatedprocs_.push_back(_rank);
        this->duplicated_points_.emplace_back(std::vector<size_t>());
        this->Nghost_.emplace_back(std::vector<size_t>());
      }
      for(const size_t &pointIdx : sentPoints[i])
      {
        this->duplicated_points_[rankIdx].push_back(pointIdx);
      }
    }
}

void Voronoi3D::CalculateInitialRadius(size_t pointsSize)
{
    this->radiuses.resize(pointsSize);
    if(this->initialRadius <= __DBL_EPSILON__)
    {
      double volume = (this->ur_[0] - this->ll_[0]) * (this->ur_[1] - this->ll_[1]) * (this->ur_[2] - this->ll_[2]);
      size_t N;
      MPI_Allreduce(&pointsSize, &N, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      this->initialRadius = 0.25 * std::pow(volume / N, 0.333333f); // heuristic
    } 
    std::fill(this->radiuses.begin(), this->radiuses.end(), this->initialRadius);
}

/**
 * \author Maor Mizrachi
 * \brief Build the voronoi, after rebalancing the points using a proper hilbert curve.
*/
void Voronoi3D::BuildHilbert(const std::vector<Vector3D> &points)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(this->firstCall)
    {
        this->CalculateInitialRadius(points.size());
    }
    if(this->radiuses.size() < points.size())
    {
        this->radiuses.resize(points.size(), this->initialRadius);
    }

    if(this->firstCall or this->CheckForRebalance(points))
    {
        this->firstCall = false;
        // first, calculate (by heuristic) the order we should sort by
        OctTree<Vector3D> tree(this->ll_, this->ur_, points);
        int depth = tree.getDepth();
        MPI_Allreduce(&depth, &depth, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        // notify the agent with the order
        this->hilbertAgent.setOrder(depth);
        // calculate the new borders
        this->hilbertAgent.determineBorders(points);
    }
    std::vector<Vector3D> new_points = this->hilbertAgent.pointsExchange(points, this->self_index_, this->sentprocs_, this->sentpoints_, this->radiuses);
    // hilbertAgent.calculateBoundingBox();
    // std::cout << "rank " << rank << ", I have " << new_points.size() << " points" << std::endl;
    this->BuildInitialize(new_points.size());

    std::vector<size_t> order;

    if(new_points.size() != 0)
    {
        std::pair<Vector3D, Vector3D> bounding_box = std::make_pair(new_points[0], new_points[0]);
        for(const Vector3D &point : new_points)
        {
            bounding_box.first.x = std::min(bounding_box.first.x, point.x);
            bounding_box.second.x = std::max(bounding_box.second.x, point.x);
            bounding_box.first.y = std::min(bounding_box.first.y, point.y);
            bounding_box.second.y = std::max(bounding_box.second.y, point.y);
            bounding_box.first.z = std::min(bounding_box.first.z, point.y);
            bounding_box.second.z = std::max(bounding_box.second.z, point.z);
        }

        // performs internal tesselation:
        // reportDuplications(new_points);
        order = HilbertOrder3D(new_points);

        // initial build for the points
        this->del_.Build(new_points, bounding_box.second, bounding_box.first, order);
    }

    // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
    this->R_.resize(this->del_.tetras_.size());
    std::fill(this->R_.begin(), this->R_.end(), -1);
    this->tetra_centers_.resize(this->R_.size());
    this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);

    if(this->radiuses.size() != new_points.size())
    {
        throw UniversalError("Rank " + std::to_string(rank) + ", wrong size of radiuses (in voronoi build) (given this->radiuses.size()=" + std::to_string(this->radiuses.size()) + " while should be new_points.size()=" + std::to_string(new_points.size()) + ")");
    }

    this->PrepareToBuildHilbert(new_points);

    CM_.resize(del_.points_.size());
    volume_.resize(Norg_);

    if(new_points.size() != 0)
    {
        // Create Voronoi
        BuildVoronoi(order);
    }
    std::vector<double>().swap(R_);
    std::vector<tetra_vec>().swap(PointTetras_);

    CalcAllCM();
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace(i))
            CalcRigidCM(i);

    // communicate the ghost CM
    
    // TODO: I suspect this part as being wrong (are the radiuses being sent correctly?)

    // vector<vector<Vector3D>> incoming = MPI_Exchange_serializable(CM_, duplicatedprocs_, duplicated_points_);
    vector<vector<Vector3D>> incoming = MPI_exchange_data(duplicatedprocs_, duplicated_points_, CM_);
    // Add the recieved CM
    for (size_t i = 0; i < incoming.size(); ++i)
        for (size_t j = 0; j < incoming.at(i).size(); ++j)
            CM_[Nghost_.at(i).at(j)] = incoming[i][j];
}

#endif // RICH_MPI

vector<vector<std::size_t>> &Voronoi3D::GetGhostIndeces(void)
{
    return Nghost_;
}

void Voronoi3D::CalcAllCM(void)
{
    std::array<Vector3D, 4> tetra;
    size_t Nfaces = FaceNeighbors_.size();
    Vector3D vtemp;
    std::vector<Vector3D> vectemp;
    double vol;
    for (size_t i = 0; i < Nfaces; ++i)
    {
        size_t N0 = FaceNeighbors_[i].first;
        size_t N1 = FaceNeighbors_[i].second;
        size_t Npoints = PointsInFace_[i].size();
        vectemp.resize(Npoints);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (size_t j = 0; j < Npoints; ++j)
            vectemp[j] = tetra_centers_[PointsInFace_[i][j]];
        Npoints -= 2;
        tetra[0] = vectemp[0];
        for (std::size_t j = 0; j < Npoints; ++j)
        {
            tetra[1] = vectemp[j + 1];
            tetra[2] = vectemp[j + 2];
            if (N1 < Norg_)
            {
                tetra[3] = del_.points_[N1];
                vol = std::abs(GetTetraVolume(tetra));
                GetTetraCM(tetra, vtemp);
                volume_[N1] += vol;
                vtemp *= vol;
                CM_[N1] += vtemp;
            }
            tetra[3] = del_.points_[N0];
            vol = std::abs(GetTetraVolume(tetra));
            GetTetraCM(tetra, vtemp);
            volume_[N0] += vol;
            vtemp *= vol;
            CM_[N0] += vtemp;
        }
    }
#ifdef __INTEL_COMPILER
    //#pragma vector aligned
#pragma omp simd
#endif
    for (size_t i = 0; i < Norg_; ++i)
        CM_[i] *= (1.0 / volume_[i]);
    // Recalc points with high aspect ratio
    for (size_t i = 0; i < Norg_; ++i)
    {
        if (fastabs(CM_[i] - del_.points_[i]) > 0.4 * GetWidth(i))
        {
            tetra[3] = CM_[i];
            CM_[i] = Vector3D();
            volume_[i] = 0;
            Nfaces = FacesInCell_[i].size();
            for (size_t k = 0; k < Nfaces; ++k)
            {
                size_t Face = FacesInCell_[i][k];
                size_t Npoints = PointsInFace_[Face].size();
                tetra[0] = tetra_centers_[PointsInFace_[Face][0]];
                for (std::size_t j = 0; j < Npoints - 2; ++j)
                {
                    tetra[1] = tetra_centers_[PointsInFace_[Face][j + 1]];
                    tetra[2] = tetra_centers_[PointsInFace_[Face][j + 2]];
                    double vol2 = std::abs(GetTetraVolume(tetra));
                    volume_[i] += vol2;
                    GetTetraCM(tetra, vtemp);
                    CM_[i] += vol2 * vtemp;
                }
            }
            CM_[i] *= (1.0 / volume_[i]);
        }
    }
}

std::pair<Vector3D, Vector3D> Voronoi3D::GetBoxCoordinates(void) const
{
    return std::pair<Vector3D, Vector3D>(ll_, ur_);
}

void Voronoi3D::BuildNoBox(vector<Vector3D> const &points, vector<vector<Vector3D>> const &ghosts, vector<size_t> toduplicate)
{
    assert(points.size() > 0);
    // Clear data
    PointTetras_.clear();
    R_.clear();
    R_.reserve(points.size());
    tetra_centers_.clear();
    tetra_centers_.reserve(points.size() * 7);
    del_.Clean();
    // Voronoi Data
    FacesInCell_.clear();
    PointsInFace_.clear();
    FaceNeighbors_.clear();
    CM_.clear();
    Face_CM_.clear();
    volume_.clear();
    area_.clear();
    Norg_ = points.size();
    duplicatedprocs_.clear();
    duplicated_points_.clear();
    Nghost_.clear();
    std::vector<size_t> order = HilbertOrder3D(points);
    del_.Build(points, ur_, ll_, order);
    for (size_t i = 0; i < ghosts.size(); ++i)
    {
        del_.BuildExtra(ghosts[i]);
    }
    vector<std::pair<size_t, size_t>> duplicate(6);
    for (size_t j = 0; j < toduplicate.size(); ++j)
    {
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
        for (size_t i = 0; i < 6; ++i)
            duplicate[i] = std::pair<size_t, size_t>(i, toduplicate[j]);
        vector<vector<size_t>> past_duplicates;
        vector<Vector3D> extra_points = CreateBoundaryPoints(duplicate, past_duplicates);
        del_.BuildExtra(extra_points);
    }

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), -1);
    tetra_centers_.resize(R_.size());
    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    CM_.resize(Norg_);
    volume_.resize(Norg_, 0);
    // Create Voronoi
    BuildVoronoi(order);

    CalcAllCM();
    CM_.resize(del_.points_.size());
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace(i))
            CalcRigidCM(i);
}

void Voronoi3D::BuildDebug(int rank)
{
    std::vector<size_t> order = read_vecst("order_" + int2str(rank) + ".bin");
    std::vector<Vector3D> points = read_vec3d("points0_" + int2str(rank) + ".bin");
    Norg_ = points.size();
    std::vector<Vector3D> bb = read_vec3d("bb_" + int2str(rank) + ".bin");
    del_.Build(points, bb[1], bb[0], order);
    points = read_vec3d("points1_" + int2str(rank) + ".bin");
    del_.BuildExtra(points);
    points = read_vec3d("points2_" + int2str(rank) + ".bin");
    del_.BuildExtra(points);
    points = read_vec3d("points3_" + int2str(rank) + ".bin");
    del_.BuildExtra(points);
    points = read_vec3d("points4_" + int2str(rank) + ".bin");
    del_.BuildExtra(points);

    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), -1);
    tetra_centers_.resize(R_.size());

    CM_.resize(del_.points_.size());
    volume_.resize(Norg_, 0);
    // Create Voronoi
    BuildVoronoi(order);

    std::vector<double>().swap(R_);
    std::vector<tetra_vec>().swap(PointTetras_);
    std::vector<Tetrahedron>().swap(del_.tetras_);

    CalcAllCM();
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace(i))
            CalcRigidCM(i);
}

void Voronoi3D::Build(vector<Vector3D> const &points)
{
    assert(points.size() > 0);
    // Clear data
    PointTetras_.clear();
    R_.clear();
    R_.reserve(points.size() * 11);
    tetra_centers_.clear();
    tetra_centers_.reserve(points.size() * 11);
    // Voronoi Data
    FacesInCell_.clear();
    PointsInFace_.clear();
    FaceNeighbors_.clear();
    CM_.clear();
    Face_CM_.clear();
    volume_.clear();
    area_.clear();
    Norg_ = points.size();
    duplicatedprocs_.clear();
    duplicated_points_.clear();
    Nghost_.clear();
    std::vector<size_t> order = HilbertOrder3D(points);
    del_.Build(points, ur_, ll_, order);

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), -1);
    tetra_centers_.resize(R_.size());
    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    vector<std::pair<std::size_t, std::size_t>> ghost_index = SerialFirstIntersections();
    vector<vector<size_t>> past_duplicates;
    vector<Vector3D> extra_points = CreateBoundaryPoints(ghost_index, past_duplicates);

    del_.BuildExtra(extra_points);

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), -1);
    tetra_centers_.resize(R_.size());
    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);
    ghost_index = SerialFindIntersections(true);
    extra_points = CreateBoundaryPoints(ghost_index, past_duplicates);
    del_.BuildExtra(extra_points);

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), -1);
    tetra_centers_.resize(R_.size());
    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);
    ghost_index = SerialFindIntersections(false);
    extra_points = CreateBoundaryPoints(ghost_index, past_duplicates);
    del_.BuildExtra(extra_points);
    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    std::vector<std::pair<size_t, size_t>>().swap(ghost_index);
    std::vector<std::vector<size_t>>().swap(past_duplicates);
    std::vector<Vector3D>().swap(extra_points);

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), -1);
    tetra_centers_.resize(R_.size());

    CM_.resize(del_.points_.size());
    volume_.resize(Norg_, 0);
    // Create Voronoi
    BuildVoronoi(order);

    std::vector<double>().swap(R_);
    std::vector<tetra_vec>().swap(PointTetras_);
    std::vector<Tetrahedron>().swap(del_.tetras_);

    CalcAllCM();
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace(i))
            CalcRigidCM(i);
}

void Voronoi3D::BuildVoronoi(std::vector<size_t> const &order)
{
    FacesInCell_.resize(Norg_);
    area_.resize(Norg_ * 10);
    Face_CM_.resize(Norg_ * 10);
    FaceNeighbors_.resize(Norg_ * 10);
    PointsInFace_.resize(Norg_ * 10);

    std::array<size_t, 128> temp, temp3;
    // Build all voronoi points
    std::size_t Ntetra = del_.tetras_.size();
    for (size_t i = 0; i < Ntetra; ++i)
        if (ShouldCalcTetraRadius(del_.tetras_[i], Norg_))
            CalcTetraRadiusCenter(i);
    // Organize the faces and assign them to cells
    std::array<double, 128> diffs, Atempvec;

    size_t FaceCounter = 0;
    boost::container::flat_set<size_t> neigh_set;
    point_vec *temp_points_in_face;
    std::array<Vector3D, 128> clean_vec;
    std::array<double, 128> area_vec_temp;

    //std::vector<Vector3D, boost::alignment::aligned_allocator<Vector3D, 32> > clean_vec;
    for (size_t i = 0; i < Norg_; ++i)
    {
        neigh_set.clear();
        neigh_set.reserve(20);
        size_t point = order[i];
        size_t ntet = PointTetras_[point].size();
        // for each point loop over its tetras
        for (size_t j = 0; j < ntet; ++j)
        {
            const size_t tetcheck = PointTetras_[point][j];
            for (size_t k = 0; k < 4; ++k)
            {
                size_t point_other = del_.tetras_[tetcheck].points[k];
                if (point_other != point && point_other > point)
                {
                    // Did we already build this face?
                    if (neigh_set.find(point_other) == neigh_set.end())
                    {
                        size_t temp_size = 0;
                        // Find all tetras for face
                        temp[0] = tetcheck;
                        ++temp_size;
                        size_t next_check = NextLoopTetra(del_.tetras_[tetcheck], tetcheck, point, point_other);
                        size_t cur_check = next_check;
                        size_t last_check = tetcheck;
                        while (next_check != tetcheck)
                        {
                            Tetrahedron const &tet_check = del_.tetras_[cur_check];
                            temp[temp_size] = cur_check;
                            ++temp_size;
                            next_check = NextLoopTetra(tet_check, last_check, point, point_other);
                            last_check = cur_check;
                            cur_check = next_check;
                        }
                        // Is face too small?
                        if (temp_size < 3)
                            continue;
                        temp_points_in_face = &PointsInFace_[FaceCounter];
                        //temp_points_in_face->reserve(8);
                        double Asize = CleanDuplicates(temp, tetra_centers_, *temp_points_in_face, ScalarProd(del_.points_[point] - del_.points_[point_other], del_.points_[point] - del_.points_[point_other]), diffs, clean_vec, temp_size);
                        if (temp_points_in_face->size() < 3)
                            continue;
                        CalcFaceAreaCM(*temp_points_in_face, tetra_centers_, clean_vec, area_[FaceCounter],
                                                     Face_CM_[FaceCounter], Atempvec);
                        if (area_[FaceCounter] < (Asize * (IsPointOutsideBox(point_other) ? 1e-14 : 1e-15)))
                            continue;
                        if (point_other >= Norg_ && point_other < (Norg_ + 4))
                        {
                            UniversalError eo("Neighboring big tet point");
                            throw eo;
                        }
                        // Make faces right handed
                        MakeRightHandFace(*temp_points_in_face, del_.points_[point], tetra_centers_, temp3, area_[FaceCounter]);
                        CleanSameLine(*temp_points_in_face, tetra_centers_, area_vec_temp);
                        FaceNeighbors_[FaceCounter].first = point;
                        FaceNeighbors_[FaceCounter].second = point_other;

                        FacesInCell_[point].push_back(FaceCounter);
                        if (point_other < Norg_)
                        {
                            FacesInCell_[point_other].push_back(FaceCounter);
                        }
                        neigh_set.insert(point_other);
                        ++FaceCounter;
                        // realloc memory if needed
                        if (FaceCounter == FaceNeighbors_.size())
                        {
                            area_.resize(static_cast<size_t>(static_cast<double>(area_.size()) * 1.25));
                            Face_CM_.resize(static_cast<size_t>(static_cast<double>(Face_CM_.size()) * 1.25));
                            FaceNeighbors_.resize(static_cast<size_t>(static_cast<double>(FaceNeighbors_.size()) * 1.25));
                            PointsInFace_.resize(static_cast<size_t>(static_cast<double>(PointsInFace_.size()) * 1.25));
                        }
                    }
                }
            }
        }
    }

    // Fix Face CM (this prevents large face velocities for close by points)
    size_t Nfaces = FaceCounter;
    Vector3D mid, norm;
    for (size_t i = 0; i < Nfaces; ++i)
    {
        mid = del_.points_[FaceNeighbors_[i].first];
        mid += del_.points_[FaceNeighbors_[i].second];
        mid *= 0.5;
        norm = del_.points_[FaceNeighbors_[i].second];
        norm -= del_.points_[FaceNeighbors_[i].first];
        Face_CM_[i] -= ScalarProd(Face_CM_[i] - mid, norm) * norm / ScalarProd(norm, norm);
    }

    area_.resize(FaceCounter);
    area_.shrink_to_fit();
    Face_CM_.resize(FaceCounter);
    Face_CM_.shrink_to_fit();
    FaceNeighbors_.resize(FaceCounter);
    FaceNeighbors_.shrink_to_fit();
    PointsInFace_.resize(FaceCounter);
    PointsInFace_.shrink_to_fit();
    for (size_t i = 0; i < Norg_; ++i)
        FacesInCell_[i].shrink_to_fit();
}

inline double Voronoi3D::GetRadius(std::size_t index)
{
    return (R_[index] = (R_[index] < 0)? CalcTetraRadiusCenter(index) : R_[index]);
}

void Voronoi3D::FindIntersectionsSingle(vector<Face> const &box, std::size_t point, Sphere &sphere,
                                                                                vector<size_t> &intersecting_faces, std::vector<double> &Rtemp, std::vector<Vector3D> &vtemp)
{
    intersecting_faces.clear();
    std::size_t N = PointTetras_[point].size();
    Rtemp.resize(N);
    vtemp.resize(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        Rtemp[i] = GetRadius(PointTetras_[point][i]);
        vtemp[i] = tetra_centers_[PointTetras_[point][i]];
    }
    size_t bsize = box.size();
    for (std::size_t j = 0; j < bsize; ++j)
    {
        Vector3D normal = CrossProduct(box[j].vertices[1] - box[j].vertices[0], box[j].vertices[2] - box[j].vertices[0]);
        normal *= (1.0 / fastsqrt(ScalarProd(normal, normal)));
        for (std::size_t i = 0; i < N; ++i)
        {
            sphere.radius = Rtemp[i];
            sphere.center = vtemp[i];
            if (FaceSphereIntersections(box[j], sphere, normal))
            {
                intersecting_faces.push_back(j);
                break;
            }
        }
    }
}

#ifdef RICH_MPI
void Voronoi3D::FindIntersectionsFirstMPI(vector<std::size_t> &res, std::size_t point,
                                                                                    Sphere &sphere, std::vector<Face> const &faces, bool &skipped, face_vec const &face_index)
{
    res.clear();
    std::size_t Ntetra = PointTetras_[point].size();
    size_t Nfaces = faces.size();
    skipped = true;
    double maxR = GetRadius(PointTetras_[point].at(0));
    for (std::size_t j = 1; j < Ntetra; ++j)
        maxR = std::max(maxR, GetRadius(PointTetras_[point][j]));
    for (size_t i = 0; i < Nfaces; ++i)
    {
        Face const &f = faces[i];
        Vector3D normal = CrossProduct(f.vertices[1] - f.vertices[0], f.vertices[2] - f.vertices[0]);
        normal *= (1.0 / fastsqrt(ScalarProd(normal, normal)));

        // Quick check if there is no intersection for sure 
        sphere.radius = 2 * maxR;
        sphere.center = GetMeshPoint(point);
        if (!FaceSphereIntersections(f, sphere, normal))
            continue;

        for (std::size_t j = 0; j < Ntetra; ++j)
        {
            sphere.radius = GetRadius(PointTetras_[point][j]);
            sphere.center = tetra_centers_[PointTetras_[point][j]];
            if (FaceSphereIntersections(f, sphere, normal))
            {
                res.push_back(face_index[i]);
                skipped = false;
                break;
            }
        }
    }
    std::sort(res.begin(), res.end());
    res = unique(res);
}

void Voronoi3D::FindIntersectionsRecursive(vector<std::size_t> &res, Tessellation3D const &tproc, std::size_t rank, std::size_t point,
                                                                                     Sphere &sphere, size_t mode, boost::container::flat_set<size_t> &visited, std::stack<std::size_t> &to_check,
                                                                                     bool &skipped, face_vec &faces, vector<size_t> &past_duplicate)
{
    res.clear();
    std::size_t N = tproc.GetPointNo();
    assert(to_check.empty());
    std::size_t Ntetra = PointTetras_[point].size();
    vector<size_t> neigh;
    if (mode == 1)
    {
        faces = tproc.GetCellFaces(rank);
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            size_t other = tproc.GetFaceNeighbors(faces[i]).first == rank ? tproc.GetFaceNeighbors(faces[i]).second : tproc.GetFaceNeighbors(faces[i]).first;
            if (other < tproc.GetPointNo() && fastabs(tproc.GetCellCM(other) - del_.points_[point]) > 50 * tproc.GetWidth(other))
                continue;
            to_check.push(faces[i]);
        }
    }
    else
    {
        if (mode == 2)
        {
            neigh.push_back(rank);
            if (past_duplicate.empty())
            {
                tproc.GetNeighbors(rank, past_duplicate);
                size_t Nn = past_duplicate.size();
                for (size_t j = 0; j < Nn; ++j)
                    past_duplicate[j] = std::min(past_duplicate[j], tproc.GetPointNo() - 1);
            }
        }
        else
        {
            tproc.GetNeighbors(rank, neigh);
            std::sort(neigh.begin(), neigh.end());
            {
                tproc.GetNeighbors(rank, past_duplicate);
                size_t Nn = past_duplicate.size();
                for (size_t j = 0; j < Nn; ++j)
                    past_duplicate[j] = std::min(past_duplicate[j], tproc.GetPointNo() - 1);
            }
        }
        for (size_t i = 0; i < past_duplicate.size(); ++i)
        {
            faces = tproc.GetCellFaces(past_duplicate[i]);
            for (size_t j = 0; j < faces.size(); ++j)
            {
                std::pair<size_t, size_t> const &fneigh = tproc.GetFaceNeighbors(faces[j]);
                if (fneigh.first != past_duplicate[i])
                {
                    if (!std::binary_search(neigh.begin(), neigh.end(), fneigh.first))
                    {
                        if (fneigh.first < tproc.GetPointNo() && fastabs(tproc.GetCellCM(fneigh.first) - del_.points_[point]) > 50 * tproc.GetWidth(fneigh.first))
                            continue;
                        to_check.push(faces[j]);
                        continue;
                    }
                }
                if (fneigh.second != past_duplicate[i])
                {
                    if (!std::binary_search(neigh.begin(), neigh.end(), fneigh.second))
                    {
                        if (fneigh.second < tproc.GetPointNo() && fastabs(tproc.GetCellCM(fneigh.second) - del_.points_[point]) > 50 * tproc.GetWidth(fneigh.second))
                            continue;
                        to_check.push(faces[j]);
                    }
                }
            }
        }
    }

    visited.clear();
    skipped = true;
    while (!to_check.empty())
    {
        std::size_t cur = to_check.top();
        to_check.pop();
        if (visited.find(cur) != visited.end())
            continue;
        visited.insert(cur);
        double maxR = GetRadius(PointTetras_[point].at(0));
        for (std::size_t j = 1; j < Ntetra; ++j)
            maxR = std::max(maxR, GetRadius(PointTetras_[point][j]));
        sphere.radius = 2 * maxR;
        sphere.center = GetMeshPoint(point);
        Face f(VectorValues(tproc.GetFacePoints(), tproc.GetPointsInFace(cur)), tproc.GetFaceNeighbors(cur).first,
                     tproc.GetFaceNeighbors(cur).second);
        Vector3D normal = CrossProduct(f.vertices[1] - f.vertices[0], f.vertices[2] - f.vertices[0]);
        normal *= (1.0 / fastsqrt(ScalarProd(normal, normal)));

        // Quick check if there is no intersection for sure
        if (!FaceSphereIntersections(f, sphere, normal))
            continue;

        for (std::size_t j = 0; j < Ntetra; ++j)
        {
            sphere.radius = GetRadius(PointTetras_[point][j]);
            sphere.center = tetra_centers_[PointTetras_[point][j]];
            if (FaceSphereIntersections(f, sphere, normal))
            {
                res.push_back(cur);
                if (mode == 1 || mode == 2)
                    skipped = false;
                if (mode == 3)
                {
                    if (f.neighbors.first < N && f.neighbors.first != rank)
                    {
                        face_vec const &faces_temp = tproc.GetCellFaces(f.neighbors.first);
                        for (std::size_t i = 0; i < faces_temp.size(); ++i)
                        {
                            if (visited.find(faces_temp[i]) == visited.end())
                            {
                                size_t other = tproc.GetFaceNeighbors(faces_temp[i]).first == f.neighbors.first ? tproc.GetFaceNeighbors(faces_temp[i]).second : tproc.GetFaceNeighbors(faces_temp[i]).first;
                                if (other < tproc.GetPointNo() && fastabs(tproc.GetCellCM(other) - del_.points_[point]) > 50 * tproc.GetWidth(other))
                                    continue;
                                to_check.push(faces_temp[i]);
                            }
                        }
                    }
                    if (f.neighbors.second < N && f.neighbors.second != rank)
                    {
                        face_vec const &faces_temp = tproc.GetCellFaces(f.neighbors.second);
                        for (std::size_t i = 0; i < faces_temp.size(); ++i)
                            if (visited.find(faces_temp[i]) == visited.end())
                            {
                                size_t other = tproc.GetFaceNeighbors(faces_temp[i]).first == f.neighbors.second ? tproc.GetFaceNeighbors(faces_temp[i]).second : tproc.GetFaceNeighbors(faces_temp[i]).first;
                                if (other < tproc.GetPointNo() && fastabs(tproc.GetCellCM(other) - del_.points_[point]) > 50 * tproc.GetWidth(other))
                                    continue;
                                to_check.push(faces_temp[i]);
                            }
                    }
                }
                break;
            }
        }
    }
    std::sort(res.begin(), res.end());
    res = unique(res);
}
#endif // RICH_MPI

void Voronoi3D::GetPointToCheck(std::size_t point, vector<unsigned char> const &checked, vector<std::size_t> &res)
{
    res.clear();
    std::size_t ntetra = PointTetras_[point].size();
    for (std::size_t i = 0; i < ntetra; ++i)
    {
        size_t tetra = PointTetras_[point][i];
        for (std::size_t j = 0; j < 4; ++j)
            if (del_.tetras_[tetra].points[j] < Norg_ && checked[del_.tetras_[tetra].points[j]] == 0)
                res.push_back(del_.tetras_[tetra].points[j]);
    }
    std::sort(res.begin(), res.end());
    res = unique(res);
}

std::size_t Voronoi3D::GetFirstPointToCheck(void) const
{
    std::size_t i;
    Tetrahedron const &tet = del_.tetras_[bigtet_];
    for (i = 0; i < 4; ++i)
        if (tet.points[i] < Norg_)
            break;
    if (i < 4)
        return tet.points[i];
    else
        throw UniversalError("Can't find first point to start boundary search");
}

#ifdef RICH_MPI
vector<std::pair<std::size_t, std::size_t>> Voronoi3D::FindIntersections(Tessellation3D const &tproc, size_t mode,
                                                                         vector<unsigned char> &checked_clear)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    vector<std::pair<std::size_t, std::size_t>> res;
    if (Norg_ == 0)
        return res;
    Sphere sphere;
    bool skipped = false;
    vector<std::size_t> intersecting_faces;
    // First attempt, only copy from neighboring cpus
    if (mode == 1)
    {
        // Create faces to check for intersection
        face_vec faces = tproc.GetCellFaces(static_cast<size_t>(rank));
        std::vector<Face> cell_faces(faces.size());
        for (size_t i = 0; i < cell_faces.size(); ++i)
            cell_faces[i] = Face(VectorValues(tproc.GetFacePoints(),
                                              tproc.GetPointsInFace(faces[i])),
                                              tproc.GetFaceNeighbors(faces[i]).first,
                                              tproc.GetFaceNeighbors(faces[i]).second);
        // Search for intersections
        for (size_t i = 0; i < Norg_; ++i)
        {
            FindIntersectionsFirstMPI(intersecting_faces, i, sphere, cell_faces, skipped, faces);
            for (std::size_t j = 0; j < intersecting_faces.size(); ++j)
                res.push_back(std::pair<std::size_t, std::size_t>(intersecting_faces[j], i));
            if (intersecting_faces.empty())
                checked_clear[i] = 1;
        }
    }
    else
    {
        std::stack<size_t> intersection_check;
        boost::container::flat_set<size_t> visited;
        vector<size_t> past_duplicate;
        face_vec vtemp;
        vector<vector<size_t>> sorted_to_duplicate = duplicated_points_;
        for (size_t i = 0; i < sorted_to_duplicate.size(); ++i)
            std::sort(sorted_to_duplicate[i].begin(), sorted_to_duplicate[i].end());
        vector<size_t> duplicatedprocs(duplicatedprocs_.size());
        for (size_t i = 0; i < duplicatedprocs_.size(); ++i)
            duplicatedprocs[i] = static_cast<size_t>(duplicatedprocs_[i]);
        for (size_t i = 0; i < Norg_; ++i)
        {
            if (checked_clear[i] == 1)
                continue;
            GetPastDuplicate(i, past_duplicate, sorted_to_duplicate, duplicatedprocs);
            FindIntersectionsRecursive(intersecting_faces, tproc, static_cast<std::size_t>(rank), i, sphere,
                                       mode, visited, intersection_check, skipped, vtemp, past_duplicate);
            for (std::size_t j = 0; j < intersecting_faces.size(); ++j)
                res.push_back(std::pair<std::size_t, std::size_t>(intersecting_faces[j], i));
            if (intersecting_faces.empty())
                checked_clear[i] = 1;
        }
    }
    return res;
}

void Voronoi3D::MPIFirstIntersections(Tessellation3D const &tproc, vector<std::pair<std::size_t, std::size_t>> &ghost_index)
{
    ghost_index.clear();
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    vector<size_t> neigh = tproc.GetNeighbors(static_cast<size_t>(rank));
    face_vec faces = tproc.GetCellFaces(static_cast<size_t>(rank));
    size_t Nneigh = neigh.size();
    vector<size_t> to_add;
    vector<Vector3D> neigh_points(Nneigh);
    vector<double> radii(Nneigh);
    for (size_t i = 0; i < Nneigh; ++i)
    {
        neigh_points[i] = tproc.GetMeshPoint(neigh[i]);
        radii[i] = neigh[i] >= tproc.GetPointNo() ? tproc.GetWidth(static_cast<size_t>(rank)) : tproc.GetWidth(neigh[i]);
    }
    size_t Ntetra = del_.tetras_.size();
    Vector3D proc_point = tproc.GetMeshPoint(static_cast<size_t>(rank));
    for (size_t i = 0; i < Ntetra; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            if (del_.tetras_[i].points[j] >= Norg_)
            {
                for (size_t k = 0; k < 4; ++k)
                {
                    if (del_.tetras_[i].points[k] < Norg_)
                    {
                        to_add.clear();
                        Vector3D const &point = del_.points_[del_.tetras_[i].points[k]];
                        double r0 = fastabs(point - proc_point);
                        size_t index = 0;
                        double mind_1 = 0;
                        for (size_t z = 0; z < Nneigh; ++z)
                        {
                            double r1 = fastabs(point - neigh_points[z]);
                            double temp = r0 > r1 ? r1 / r0 : r0 / r1;
                            if (temp > 0.75 && r1 < 2.5 * radii[z])
                                to_add.push_back(z);
                            if (r1 * mind_1 < 1)
                            {
                                mind_1 = 1.0 / r1;
                                index = z;
                            }
                        }
                        if (mind_1 * radii[index] > 0.25)
                            to_add.push_back(index);
                        std::sort(to_add.begin(), to_add.end());
                        to_add = unique(to_add);
                        for (size_t l = 0; l < to_add.size(); ++l)
                            ghost_index.push_back(std::pair<size_t, size_t>(faces[to_add[l]], del_.tetras_[i].points[k]));
                    }
                }
                break;
            }
        }
    }
}
#endif //RICH_MPI

vector<std::pair<std::size_t, std::size_t>> Voronoi3D::SerialFirstIntersections(void)
{
    vector<Face> box;
    vector<Vector3D> normals;
    this->InitialBoxBuild(box, normals);
    size_t Nfaces = box.size();

    //    vector<std::size_t> point_neigh;
    vector<std::pair<std::size_t, std::size_t>> res;
    Sphere sphere;
    vector<unsigned char> will_check(Norg_, 0);
    std::size_t cur_loc;
    std::stack<std::size_t> check_stack;
    FirstCheckList(check_stack, will_check, Norg_, del_, PointTetras_);
    std::vector<double> vdist(Nfaces);
    std::vector<Vector3D> vtemp(Nfaces);
    while (!check_stack.empty())
    {
        cur_loc = check_stack.top();
        check_stack.pop();
        double inv_max = 0;
        size_t max_loc = 0;
        size_t j = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (; j < Nfaces; ++j)
        {
            vtemp[j] = del_.points_[cur_loc];
            vtemp[j] -= box[j].vertices[0];
            double sprod = vtemp[j].x * normals[j].x + vtemp[j].y * normals[j].y + vtemp[j].z * normals[j].z;
            vdist[j] = 1.0 / std::abs(sprod);
        }
        j = 0;
        for (; j < Nfaces; ++j)
        {
            if (vdist[j] > inv_max)
            {
                inv_max = vdist[j];
                max_loc = j;
            }
        }
        res.push_back(std::pair<std::size_t, std::size_t>(max_loc, cur_loc));
    }
    return res;
}

vector<std::pair<std::size_t, std::size_t>> Voronoi3D::SerialFindIntersections(bool first_run)
{
    if (Norg_ < 50)
    {
        vector<std::pair<std::size_t, std::size_t>> res;
        size_t const Nfaces = box_faces_.empty() ? 6 : box_faces_.size();
        res.reserve(Norg_ * Nfaces);
        for (size_t i = 0; i < Norg_; ++i)
            for (size_t j = 0; j < Nfaces; ++j)
                res.push_back(std::pair<std::size_t, std::size_t>(j, i));
        return res;
    }
    std::stack<std::size_t> check_stack;
    vector<Face> box = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<std::size_t> point_neigh;
    vector<std::pair<std::size_t, std::size_t>> res;
    Sphere sphere;
    vector<unsigned char> checked(Norg_, 0), will_check(Norg_, 0);
    std::size_t cur_loc;
    if (first_run)
    {
        FirstCheckList(check_stack, will_check, Norg_, del_, PointTetras_);
        //cur_loc = check_stack.top();
        check_stack.pop();
    }
    else
    {
        cur_loc = GetFirstPointToCheck();
        check_stack.push(cur_loc);
        will_check[cur_loc] = true;
    }
    vector<size_t> intersecting_faces;
    std::vector<double> Rtemp;
    std::vector<Vector3D> vtemp;
    while (!check_stack.empty())
    {
        cur_loc = check_stack.top();
        check_stack.pop();
        checked[cur_loc] = true;
        // Does sphere have any intersections?
        bool added = false;
        FindIntersectionsSingle(box, cur_loc, sphere, intersecting_faces, Rtemp, vtemp);
        if (!intersecting_faces.empty())
        {
            added = true;
            for (std::size_t j = 0; j < intersecting_faces.size(); ++j)
                res.push_back(std::pair<std::size_t, std::size_t>(intersecting_faces[j], cur_loc));
        }
        if (added && !first_run)
        {
            GetPointToCheck(cur_loc, checked, point_neigh);
            std::size_t Nneigh = point_neigh.size();
            for (std::size_t j = 0; j < Nneigh; ++j)
                if (point_neigh[j] < Norg_ && !will_check[point_neigh[j]])
                {
                    check_stack.push(point_neigh[j]);
                    will_check[point_neigh[j]] = true;
                }
        }
    }
    return res;
}

double Voronoi3D::CalcTetraRadiusCenter(std::size_t index)
{
    Vector3D v2(del_.points_[del_.tetras_[index].points[1]]);
    v2 -= del_.points_[del_.tetras_[index].points[0]];
    Vector3D v3(del_.points_[del_.tetras_[index].points[2]]);
    v3 -= del_.points_[del_.tetras_[index].points[0]];
    Vector3D v4(del_.points_[del_.tetras_[index].points[3]]);
    v4 -= del_.points_[del_.tetras_[index].points[0]];

    Mat33<double> m_a(v2.x, v2.y, v2.z,
                                        v3.x, v3.y, v3.z,
                                        v4.x, v4.y, v4.z);
    double a = m_a.determinant();
    if(std::abs(a) < 100 * std::numeric_limits<double>::min())
        return CalcTetraRadiusCenterHiPrecision(index);
    Mat33<double> m_Dx(ScalarProd(v2, v2), v2.y, v2.z,
                                         ScalarProd(v3, v3), v3.y, v3.z,
                                         ScalarProd(v4, v4), v4.y, v4.z);
    double DDx = m_Dx.determinant();

    Mat33<double> m_Dy(ScalarProd(v2, v2), v2.x, v2.z,
                                         ScalarProd(v3, v3), v3.x, v3.z,
                                         ScalarProd(v4, v4), v4.x, v4.z);
    double DDy = -m_Dy.determinant();

    Mat33<double> m_Dz(ScalarProd(v2, v2), v2.x, v2.y,
                                         ScalarProd(v3, v3), v3.x, v3.y,
                                         ScalarProd(v4, v4), v4.x, v4.y);
    double DDz = m_Dz.determinant();
    Vector3D center = Vector3D(DDx / (2 * a), DDy / (2 * a), DDz / (2 * a)) + del_.points_[del_.tetras_[index].points[0]];
    tetra_centers_[index] = center;
    double Rres = 0.5 * std::sqrt(DDx * DDx + DDy * DDy + DDz * DDz) / std::abs(a);
    // Sanity check
    /*double Rcheck0 = fastabs(del_.points_[del_.tetras_[index].points[0]] - center);
        double Rcheck1 = fastabs(del_.points_[del_.tetras_[index].points[1]] - center);
        double Rcheck2 = fastabs(del_.points_[del_.tetras_[index].points[2]] - center);
        double Rcheck3 = fastabs(del_.points_[del_.tetras_[index].points[3]] - center);*/
    Vector3D v1(del_.points_[del_.tetras_[index].points[0]]);
    double Rcheck0 = fastabs(v1 - center);
    v2 += v1;
    v2 -= center;
    double Rcheck1 = fastabs(v2);
    v3 += v1;
    v3 -= center;
    double Rcheck2 = fastabs(v3);
    v4 += v1;
    v4 -= center;
    double Rcheck3 = fastabs(v4);
    double tol = 1 + 1e-6;
    if (((Rcheck0 + Rcheck1 + Rcheck2 + Rcheck3) * tol < (4 * Rcheck0)) || ((Rcheck0 + Rcheck1 + Rcheck2 + Rcheck3) > (tol * 4 * Rcheck0)))
        return CalcTetraRadiusCenterHiPrecision(index);
    if (Rcheck0 > tol * Rres || Rcheck0 * tol < Rres)
        return CalcTetraRadiusCenterHiPrecision(index);
    double const a_tol = 1e-8;
    if (std::abs(a) < Rres * Rres * Rres * a_tol)
        return CalcTetraRadiusCenterHiPrecision(index);
    return Rres;
}

double Voronoi3D::CalcTetraRadiusCenterHiPrecision(std::size_t index)
{
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V0;
    V0[0] = del_.points_[del_.tetras_[index].points[0]].x;
    V0[1] = del_.points_[del_.tetras_[index].points[0]].y;
    V0[2] = del_.points_[del_.tetras_[index].points[0]].z;
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V2;
    V2[0] = del_.points_[del_.tetras_[index].points[1]].x;
    V2[1] = del_.points_[del_.tetras_[index].points[1]].y;
    V2[2] = del_.points_[del_.tetras_[index].points[1]].z;
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V3;
    V3[0] = del_.points_[del_.tetras_[index].points[2]].x;
    V3[1] = del_.points_[del_.tetras_[index].points[2]].y;
    V3[2] = del_.points_[del_.tetras_[index].points[2]].z;
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V4;
    V4[0] = del_.points_[del_.tetras_[index].points[3]].x;
    V4[1] = del_.points_[del_.tetras_[index].points[3]].y;
    V4[2] = del_.points_[del_.tetras_[index].points[3]].z;
    V2[0] -= V0[0];
    V2[1] -= V0[1];
    V2[2] -= V0[2];
    V3[0] -= V0[0];
    V3[1] -= V0[1];
    V3[2] -= V0[2];
    V4[0] -= V0[0];
    V4[1] -= V0[1];
    V4[2] -= V0[2];
    std::array<boost::multiprecision::cpp_dec_float_50, 9> mat;
    mat[0] = V2[0];
    mat[1] = V2[1];
    mat[2] = V2[2];
    mat[3] = V3[0];
    mat[4] = V3[1];
    mat[5] = V3[2];
    mat[6] = V4[0];
    mat[7] = V4[1];
    mat[8] = V4[2];
    boost::multiprecision::cpp_dec_float_50 ba = Calc33Det(mat);
    mat[0] = V2[0] * V2[0] + V2[1] * V2[1] + V2[2] * V2[2];
    mat[1] = V2[1];
    mat[2] = V2[2];
    mat[3] = V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2];
    mat[4] = V3[1];
    mat[5] = V3[2];
    mat[6] = V4[0] * V4[0] + V4[1] * V4[1] + V4[2] * V4[2];
    mat[7] = V4[1];
    mat[8] = V4[2];
    boost::multiprecision::cpp_dec_float_50 bDx = Calc33Det(mat);
    mat[0] = V2[0] * V2[0] + V2[1] * V2[1] + V2[2] * V2[2];
    mat[1] = V2[0];
    mat[2] = V2[2];
    mat[3] = V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2];
    mat[4] = V3[0];
    mat[5] = V3[2];
    mat[6] = V4[0] * V4[0] + V4[1] * V4[1] + V4[2] * V4[2];
    mat[7] = V4[0];
    mat[8] = V4[2];
    boost::multiprecision::cpp_dec_float_50 bDy = -Calc33Det(mat);
    mat[0] = V2[0] * V2[0] + V2[1] * V2[1] + V2[2] * V2[2];
    mat[1] = V2[0];
    mat[2] = V2[1];
    mat[3] = V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2];
    mat[4] = V3[0];
    mat[5] = V3[1];
    mat[6] = V4[0] * V4[0] + V4[1] * V4[1] + V4[2] * V4[2];
    mat[7] = V4[0];
    mat[8] = V4[1];
    boost::multiprecision::cpp_dec_float_50 bDz = Calc33Det(mat);
    boost::multiprecision::cpp_dec_float_50 temp = (bDx / (2 * ba) + V0[0]);
    tetra_centers_[index].x = temp.convert_to<double>();
    temp = (bDy / (2 * ba) + V0[1]);
    tetra_centers_[index].y = temp.convert_to<double>();
    temp = (bDz / (2 * ba) + V0[2]);
    tetra_centers_[index].z = temp.convert_to<double>();
    temp = (boost::multiprecision::sqrt(bDx * bDx + bDy * bDy + bDz * bDz) / ba);
    return 0.5 * temp.convert_to<double>();
}

void Voronoi3D::GetTetraCM(std::array<Vector3D, 4> const &points, Vector3D &CM) const
{
    double x = 0, y = 0, z = 0;
    //CM.Set(0, 0, 0);
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+ \
                                                     : x, y, z)
#endif
    for (std::size_t i = 0; i < 4; i++)
    {
        x += points[i].x;
        y += points[i].y;
        z += points[i].z;
    }
    CM.Set(x, y, z);
    CM *= 0.25;
}

double Voronoi3D::GetTetraVolume(std::array<Vector3D, 4> const &points) const
{
    return std::abs(orient3d(points)) / 6.0;
}

/*
void Voronoi3D::CalcCellCMVolume(std::size_t index)
{
    volume_[index] = 0;
    CM_[index] = Vector3D();
    std::size_t Nfaces = FacesInCell_[index].size();
    std::array<Vector3D, 4> tetra;
    tetra[3] = del_.points_[index];
    Vector3D vtemp;
    for (std::size_t i = 0; i < Nfaces; ++i)
        {
            std::size_t face = FacesInCell_[index][i];
            std::size_t Npoints = PointsInFace_[face].size();
            tetra[0] = tetra_centers_[PointsInFace_[face][0]];
            double fvol = 0;
            for (std::size_t j = 0; j < Npoints - 2; ++j)
	{
	    tetra[1] = tetra_centers_[PointsInFace_[face][j + 1]];
	    tetra[2] = tetra_centers_[PointsInFace_[face][j + 2]];
	    double vol = GetTetraVolume(tetra);
	    fvol += std::abs(vol);
	    GetTetraCM(tetra, vtemp);
	    CM_[index] += std::abs(vol)*vtemp;
	}
            volume_[index] += fvol;
        }
    CM_[index] = CM_[index] / volume_[index];
}
*/

void Voronoi3D::output(std::string const &filename) const
{

    std::ofstream file_handle(filename.c_str(), std::ios::out | std::ios::binary);
    assert(file_handle.is_open());
    binary_write_single_int(static_cast<int>(Norg_), file_handle);

    // Points
    for (std::size_t i = 0; i < Norg_; ++i)
    {
        binary_write_single_double(del_.points_[i].x, file_handle);
        binary_write_single_double(del_.points_[i].y, file_handle);
        binary_write_single_double(del_.points_[i].z, file_handle);
    }

    binary_write_single_int(static_cast<int>(tetra_centers_.size()), file_handle);
    // Face Points
    for (std::size_t i = 0; i < tetra_centers_.size(); ++i)
    {
        binary_write_single_double(tetra_centers_[i].x, file_handle);
        binary_write_single_double(tetra_centers_[i].y, file_handle);
        binary_write_single_double(tetra_centers_[i].z, file_handle);
    }

    // Faces in cell
    for (std::size_t i = 0; i < Norg_; ++i)
    {
        binary_write_single_int(static_cast<int>(FacesInCell_[i].size()), file_handle);
        for (std::size_t j = 0; j < FacesInCell_[i].size(); ++j)
            binary_write_single_int(static_cast<int>(FacesInCell_[i][j]), file_handle);
    }

    // Points in Face
    binary_write_single_int(static_cast<int>(PointsInFace_.size()), file_handle);
    for (std::size_t i = 0; i < PointsInFace_.size(); ++i)
    {
        binary_write_single_int(static_cast<int>(PointsInFace_[i].size()), file_handle);
        for (std::size_t j = 0; j < PointsInFace_[i].size(); ++j)
            binary_write_single_int(static_cast<int>(PointsInFace_[i][j]), file_handle);
    }

    file_handle.close();
}

#ifdef RICH_MPI
void Voronoi3D::output_buildextra(std::string const &filename) const
{
    std::ofstream file_handle(filename.c_str(), std::ios::out | std::ios::binary);
    assert(file_handle.is_open());
    size_t stemp = Norg_;
    binary_write_single_int(static_cast<int>(stemp), file_handle);
    stemp = del_.points_.size();
    binary_write_single_int(static_cast<int>(stemp), file_handle);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Points
    for (std::size_t i = 0; i < stemp; ++i)
    {
        binary_write_single_double(del_.points_[i].x, file_handle);
        binary_write_single_double(del_.points_[i].y, file_handle);
        binary_write_single_double(del_.points_[i].z, file_handle);
    }

    binary_write_single_int(static_cast<int>(duplicatedprocs_.size()), file_handle);
    // Procs
    assert(duplicatedprocs_.size() == Nghost_.size());
    for (size_t i = 0; i < duplicatedprocs_.size(); ++i)
    {
        binary_write_single_int(static_cast<int>(duplicatedprocs_[i]), file_handle);
        binary_write_single_int(static_cast<int>(Nghost_[i].size()), file_handle);
        for (size_t j = 0; j < Nghost_[i].size(); ++j)
            binary_write_single_int(static_cast<int>(Nghost_[i][j]), file_handle);
    }
    file_handle.close();
}
#endif

std::size_t Voronoi3D::GetPointNo(void) const
{
    return Norg_;
}

Vector3D Voronoi3D::GetMeshPoint(std::size_t index) const
{
    return del_.points_[index];
}

double Voronoi3D::GetArea(std::size_t index) const
{
    return area_[index];
}

Vector3D const &Voronoi3D::GetCellCM(std::size_t index) const
{
    return CM_[index];
}

std::size_t Voronoi3D::GetTotalFacesNumber(void) const
{
    return FaceNeighbors_.size();
}

double Voronoi3D::GetWidth(std::size_t index) const
{
    return std::pow(3 * volume_[index] * 0.25 / M_PI, 0.3333333333);
}

double Voronoi3D::GetVolume(std::size_t index) const
{
    return volume_[index];
}

face_vec const &Voronoi3D::GetCellFaces(std::size_t index) const
{
    return FacesInCell_[index];
}

vector<Vector3D> &Voronoi3D::accessMeshPoints(void)
{
    return del_.points_;
}

const vector<Vector3D> &Voronoi3D::getMeshPoints(void) const
{
    return del_.points_;
}

vector<std::size_t> Voronoi3D::GetNeighbors(std::size_t index) const
{
    const size_t N = FacesInCell_[index].size();
    vector<size_t> res(N);
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
    for (size_t i = 0; i < N; ++i)
    {
        size_t face = FacesInCell_[index][i];
        res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second : FaceNeighbors_[face].first;
    }
    return res;
}

void Voronoi3D::GetNeighbors(size_t index, vector<size_t> &res) const
{
    std::size_t N = FacesInCell_[index].size();
    res.resize(N);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
    for (std::size_t i = 0; i < N; ++i)
    {
        std::size_t face = FacesInCell_[index][i];
        res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second : FaceNeighbors_[face].first;
    }
}

Tessellation3D *Voronoi3D::clone(void) const
{
    return new Voronoi3D(*this);
}

Voronoi3D::Voronoi3D(Voronoi3D const &other) : ll_(other.ll_), ur_(other.ur_), Norg_(other.Norg_), bigtet_(other.bigtet_),
                                                set_temp_(other.set_temp_), stack_temp_(other.stack_temp_), del_(other.del_), PointTetras_(other.PointTetras_), R_(other.R_),
                                                tetra_centers_(other.tetra_centers_), FacesInCell_(other.FacesInCell_), PointsInFace_(other.PointsInFace_),
                                                FaceNeighbors_(other.FaceNeighbors_), CM_(other.CM_), Face_CM_(other.Face_CM_), volume_(other.volume_), area_(other.area_),
                                                duplicated_points_(other.duplicated_points_), sentprocs_(other.sentprocs_), duplicatedprocs_(other.duplicatedprocs_), sentpoints_(other.sentpoints_),
                                                Nghost_(other.Nghost_), self_index_(other.self_index_), temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>()), box_faces_(other.box_faces_)
                                                #ifdef RICH_MPI
<<<<<<< HEAD
                                                , hilbertAgent(other.hilbertAgent), initialRadius(other.initialRadius), firstCall(true) 
=======
                                                ,hilbertAgent(other.hilbertAgent), initialRadius(other.initialRadius)
>>>>>>> 342f02048afef504c85b1f72ffc45f1cefee8c24
                                                #endif // RICH_MPI
                                                {}

bool Voronoi3D::NearBoundary(std::size_t index) const
{
    std::size_t N = FacesInCell_[index].size();
    for (std::size_t i = 0; i < N; ++i)
    {
        if (BoundaryFace(FacesInCell_[index][i]))
            return true;
    }
    return false;
}

bool Voronoi3D::IsPointOutsideBox(size_t index) const
{
    if(box_faces_.empty())
        return !PointInDomain(ll_, ur_, del_.points_[index]);
    else
        return !PointInPoly(box_faces_, del_.points_[index]);
}

bool Voronoi3D::BoundaryFace(std::size_t index) const
{
    if (FaceNeighbors_[index].first >= Norg_ || FaceNeighbors_[index].second >= Norg_)
    {
#ifdef RICH_MPI
        if(box_faces_.empty())
        {
            if (PointInDomain(ll_, ur_, del_.points_[std::max(FaceNeighbors_[index].first, FaceNeighbors_[index].second)]))
                return false;
            else
                return true;
        }
        else
            if(PointInPoly(box_faces_, del_.points_[std::max(FaceNeighbors_[index].first, FaceNeighbors_[index].second)]))
                return false;
            else
#endif
            return true;
    }
    else
        return false;
}

vector<vector<std::size_t>> &Voronoi3D::GetDuplicatedPoints(void)
{
    return duplicated_points_;
}

vector<vector<std::size_t>> const &Voronoi3D::GetDuplicatedPoints(void) const
{
    return duplicated_points_;
}

std::size_t Voronoi3D::GetTotalPointNumber(void) const
{
    return del_.points_.size();
}

vector<Vector3D> &Voronoi3D::GetAllCM(void)
{
    return CM_;
}

vector<Vector3D> Voronoi3D::GetAllCM(void) const
{
    return CM_;
}

void Voronoi3D::GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point) const
{
    result.clear();
    result.reserve(70);
    vector<std::size_t> neigh = GetNeighbors(point);
    result = neigh;
    std::size_t N = neigh.size();
    std::sort(neigh.begin(), neigh.end());
    vector<std::size_t> temp;
    for (std::size_t i = 0; i < N; ++i)
    {
        if (neigh[i] < Norg_)
        {
            temp = GetNeighbors(neigh[i]);
            result.insert(result.end(), temp.begin(), temp.end());
        }
    }
    std::sort(result.begin(), result.end());
    result = unique(result);
    result = RemoveList(result, neigh);
    RemoveVal(result, point);
}

vector<boost::container::small_vector<size_t, 8>> &Voronoi3D::GetAllPointsInFace(void)
{
    return PointsInFace_;
}

vector<boost::container::small_vector<size_t, 8>> const& Voronoi3D::GetAllPointsInFace(void)const
{
    return PointsInFace_;
}

size_t &Voronoi3D::GetPointNo(void)
{
    return Norg_;
}

std::vector<std::pair<size_t, size_t>> &Voronoi3D::GetAllFaceNeighbors(void)
{
    return FaceNeighbors_;
}

vector<double> &Voronoi3D::GetAllVolumes(void)
{
    return volume_;
}

vector<double> Voronoi3D::GetAllVolumes(void) const
{
    return volume_;
}

Vector3D Voronoi3D::Normal(std::size_t faceindex) const
{
    return del_.points_[FaceNeighbors_[faceindex].second] - del_.points_[FaceNeighbors_[faceindex].first];
}

bool Voronoi3D::IsGhostPoint(std::size_t index) const
{
    return index >= Norg_;
}

Vector3D Voronoi3D::FaceCM(std::size_t index) const
{
    return Face_CM_[index];
}

Vector3D Voronoi3D::CalcFaceVelocity(std::size_t index, Vector3D const &v0, Vector3D const &v1) const
{
    std::size_t p0 = FaceNeighbors_[index].first;
    std::size_t p1 = FaceNeighbors_[index].second;
    Vector3D r0 = GetMeshPoint(p0);
    Vector3D r1 = GetMeshPoint(p1);
    Vector3D r_diff = r1 - r0;
    double abs_r_diff = ScalarProd(r_diff, r_diff);

    Vector3D f = FaceCM(index);
    r1 += r0;
    r1 *= 0.5;
    f -= r1;
    Vector3D delta_w = ScalarProd((v0 - v1), f) * r_diff / abs_r_diff;
#ifdef RICH_DEBUG
    double dw_abs = fastabs(delta_w);
#endif // RICH_DEBUG
    Vector3D w = (v0 + v1) * 0.5;
#ifdef RICH_DEBUG
    //    double w_abs = std::max(fastabs(v0),fastabs(v1));
#endif // RICH_DEBUG
    //if (dw_abs > w_abs)
    //	delta_w *= (1 + (std::atan(dw_abs / w_abs) - 0.25 * M_PI)*2) * (w_abs / dw_abs);
#ifdef RICH_DEBUG
    if (!std::isfinite(dw_abs))
    {
        r0 = GetMeshPoint(p0);
        r1 = GetMeshPoint(p1);
        f = FaceCM(index);
        UniversalError eo("Bad Face velocity");
        eo.addEntry("Face index", index);
        eo.addEntry("Neigh 0", p0);
        eo.addEntry("Neigh 1", p1);
        eo.addEntry("Neigh 0 x", r0.x);
        eo.addEntry("Neigh 0 y", r0.y);
        eo.addEntry("Neigh 0 z", r0.z);
        eo.addEntry("Neigh 0 CMx", CM_[p0].x);
        eo.addEntry("Neigh 0 CMy", CM_[p0].y);
        eo.addEntry("Neigh 0 CMz", CM_[p0].z);
        eo.addEntry("Neigh 1 x", r1.x);
        eo.addEntry("Neigh 1 y", r1.y);
        eo.addEntry("Neigh 1 z", r1.z);
        eo.addEntry("Neigh 1 CMx", CM_[p1].x);
        eo.addEntry("Neigh 1 CMy", CM_[p1].y);
        eo.addEntry("Neigh 1 CMz", CM_[p1].z);
        eo.addEntry("Face CMx", f.x);
        eo.addEntry("Face CMy", f.y);
        eo.addEntry("Face CMz", f.z);
        eo.addEntry("V0x", v0.x);
        eo.addEntry("V0y", v0.y);
        eo.addEntry("V0z", v0.z);
        eo.addEntry("V1x", v1.x);
        eo.addEntry("V1y", v1.y);
        eo.addEntry("V1z", v1.z);
        throw eo;
    }
#endif
    w += delta_w;
    return w;
}

vector<double> &Voronoi3D::GetAllArea(void)
{
    return area_;
}

vector<Vector3D> &Voronoi3D::GetAllFaceCM(void)
{
    return Face_CM_;
}

vector<face_vec> &Voronoi3D::GetAllCellFaces(void)
{
    return FacesInCell_;
}

vector<face_vec> const& Voronoi3D::GetAllCellFaces(void)const
{
    return FacesInCell_;
}

vector<Vector3D> &Voronoi3D::GetFacePoints(void)
{
    return tetra_centers_;
}

vector<Vector3D> const &Voronoi3D::GetFacePoints(void) const
{
    return tetra_centers_;
}

point_vec const &Voronoi3D::GetPointsInFace(std::size_t index) const
{
    return PointsInFace_[index];
}

std::pair<std::size_t, std::size_t> Voronoi3D::GetFaceNeighbors(std::size_t face_index) const
{
    return std::pair<std::size_t, std::size_t>(FaceNeighbors_[face_index]);
}

vector<int> Voronoi3D::GetDuplicatedProcs(void) const
{
    return duplicatedprocs_;
}

vector<int> Voronoi3D::GetSentProcs(void) const
{
    return sentprocs_;
}

vector<vector<std::size_t>> const &Voronoi3D::GetSentPoints(void) const
{
    return sentpoints_;
}

vector<std::size_t> const &Voronoi3D::GetSelfIndex(void) const
{
    return self_index_;
}

vector<int> &Voronoi3D::GetSentProcs(void)
{
    return sentprocs_;
}

vector<vector<std::size_t>> &Voronoi3D::GetSentPoints(void)
{
    return sentpoints_;
}

vector<std::size_t> &Voronoi3D::GetSelfIndex(void)
{
    return self_index_;
}

void Voronoi3D::SetBox(Vector3D const &ll, Vector3D const &ur)
{
    ll_ = ll;
    ur_ = ur;
}
