
#include "trees.hpp"
#include <boost/random.hpp>
#include <algorithm>

Trees::Trees(): VC_kd_sharp_corners(),
                VC_kd_sharp_edges(),
                VC_kd_faces(),
                ball_kd_vertices(),
                ball_kd_edges(),
                ball_kd_faces() {}

void Trees::loadPLC(PL_Complex const& plc, std::size_t const Nsample_edges, std::size_t const Nsample_faces){
    
    auto const& [sharp_corners_points, i_corners] = pointsFromVertices(plc.sharp_corners);
    auto const& [p_edges, v_edges, i_feature_edges, i_plc_edges] = superSampleEdges(plc.sharp_edges, Nsample_edges);
    auto const& [p_faces, v_faces, i_feature_faces, i_plc_faces] = superSampleFaces(plc.faces, Nsample_faces);
    
    VC_kd_sharp_corners = VoroCrust_KD_Tree_Boundary(sharp_corners_points);
    VC_kd_sharp_edges = VoroCrust_KD_Tree_Boundary(p_edges, v_edges, i_feature_edges, i_plc_edges);
    VC_kd_faces = VoroCrust_KD_Tree_Boundary(p_faces, v_faces, i_feature_faces, i_plc_faces);
    
}

std::tuple<std::vector<Vector3D>, std::vector<std::size_t>> Trees::pointsFromVertices(std::vector<Vertex> const& vertices){
    std::size_t const Npoints = vertices.size();
    std::vector<Vector3D> points(Npoints, {0, 0, 0});
    std::vector<std::size_t> plc_index(Npoints, 0);

    for(std::size_t i = 0; i<Npoints; ++i){
        points[i] = vertices[i]->vertex;
        plc_index[i] = vertices[i]->index;
    }
    
    return std::tuple(points, plc_index);
}

std::tuple<std::vector<Vector3D>, std::vector<Vector3D>, std::vector<std::size_t>, std::vector<std::size_t>> Trees::superSampleEdges(std::vector<Edge> const& edges, std::size_t const Nsample){
    // generate a random number generator
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937, boost::uniform_01<>> rand_gen(rng, zeroone);

    std::vector<double> start_len(edges.size(), 0.0);

    // calculate the total length of all edges and the start length of each individual edge
    // i.e. the total length up to it. 
    double total_len = 0;
    for(std::size_t i=0; i<edges.size(); ++i){
        Edge const& edge = edges[i];
        start_len[i] = total_len;

        double const edge_len = abs(edge->vertex2->vertex - edge->vertex1->vertex);
        total_len += edge_len;
    }

    std::cout << "\nEdge Samples : \n---------------------\n";
    std::vector<Vector3D> points(Nsample, {0, 0, 0});
    std::vector<Vector3D> parallel(Nsample, {0, 0, 0});
    std::vector<std::size_t> feature_index(Nsample, 0);
    std::vector<std::size_t> plc_index(Nsample, 0);

    for(std::size_t i=0; i<Nsample; ++i){
        double const sample = rand_gen()*total_len; // sample a point uniformly [0, total_length)

        // find the index of the edge it lies
        auto const iter_lower_bound = std::lower_bound(start_len.begin(), start_len.end(), sample);
        std::size_t edge_index = std::distance(start_len.begin(), iter_lower_bound) - 1;

        // if sample is on a Vertex resample
        if(std::abs(start_len[edge_index] - sample) < 1e-14){
            //! TODO: make eps a user given parameter. 
            i--;
            continue;
        }

        // find point on edge
        Edge const& edge = edges[edge_index];
        Vector3D const& edge_vec = (edge->vertex2->vertex - edge->vertex1->vertex);
        double const factor = (sample - start_len[edge_index]) / abs(edge_vec);
        Vector3D const& point = edge->vertex1->vertex + factor*edge_vec;

        points[i] = point;
        parallel[i] = (edge_vec / abs(edge_vec));
        feature_index[i] = edge->crease_index;
        plc_index[i] = edge->index;
    }

    return std::tuple(points, parallel, feature_index, plc_index);
}

std::tuple<std::vector<Vector3D>, std::vector<Vector3D>, std::vector<std::size_t>, std::vector<std::size_t>> Trees::superSampleFaces(std::vector<Face> const& faces, std::size_t const Nsample){
    // generate a random number generator
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937, boost::uniform_01<>> rand_gen(rng, zeroone);

    std::vector<double> start_area(faces.size(), 0.0);

    // calculate the total area and the start area of each face
    // i.e. the total area up to the face
    double total_area = 0.0;
    for(std::size_t i=0; i<faces.size(); ++i){
        Face const& face = faces[i];
        start_area[i] = total_area;

        double const face_area = face->calcArea();
        total_area += face_area;
    }

    std::cout << "\nFace Samples: \n---------------------\n";
    std::vector<Vector3D> points(Nsample, {0, 0, 0});
    std::vector<Vector3D> normals(Nsample, {0, 0, 0});
    std::vector<std::size_t> feature_index(Nsample, 0);
    std::vector<std::size_t> plc_index(Nsample, 0);

    for(std::size_t i = 0; i<Nsample; ++i){
        double const sample_area = rand_gen()*total_area; // sample a number uniformly in [0, total_area)

        // find the face it lies on
        auto const iter_lower_bound = std::lower_bound(start_area.begin(), start_area.end(), sample_area);
        auto face_index = std::distance(start_area.begin(), iter_lower_bound) - 1;

        Face const& face = faces[face_index];
        

        //! TODO: add option for sampling for a face with more than 3 edges.
        if(face->vertices.size() != 3){
            std::cout << "ERROR: Algorith supports only triangular meshes for now" << std::endl;
            exit(1);
        }

        // sample a point in a triangle
        double const sqrt_r1 = std::sqrt(rand_gen());
        double const r2 = rand_gen();

        if(sqrt_r1 < 1e-14 || r2 < 1e-14){
            //! TODO: make eps a user given parameter. 
            i--;
            continue;
        }

        Vector3D const& A = face->vertices[0]->vertex;
        Vector3D const& B = face->vertices[1]->vertex;
        Vector3D const& C = face->vertices[2]->vertex;

        // sample point formula is taken from https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle-in-3d
        Vector3D const& point = (1.0-sqrt_r1)*A + (sqrt_r1*(1.0-r2))*B + (r2*sqrt_r1)*C;

        points[i] = point;
        normals[i] = face->calcNormal();
        feature_index[i] = face->patch_index;
        plc_index[i] = face->index;
    }

    return std::tuple(points, normals, feature_index, plc_index);
}