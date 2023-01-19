
#include "trees.hpp"
#include <boost/random.hpp>
#include <algorithm>

Trees::Trees(): kd_vertices(nullptr), 
                kd_edges(nullptr), 
                kd_faces(nullptr),
                vertices_points(),
                edges_points(),
                faces_points() {}

void Trees::loadPLC(PL_Complex const& plc, std::size_t const Nsample_edges, std::size_t const Nsample_faces){
    
    std::size_t const Npoints = plc.vertices.size();
    vertices_points = pointsFromVertices(plc.vertices);

    kd_vertices = std::make_shared<ANNkd_tree>(vertices_points, Npoints, 1, ANN_KD_SUGGEST);

    edges_points = superSampleEdges(plc.edges, Nsample_edges);

    kd_edges = std::make_shared<ANNkd_tree>(edges_points, Nsample_edges, 1, ANN_KD_SUGGEST);

    faces_points = superSampleFaces(plc.faces, Nsample_faces);

    kd_faces = std::make_shared<ANNkd_tree>(faces_points, Nsample_faces, 1, ANN_KD_SUGGEST);    
}

ANNpointArray Trees::pointsFromVertices(std::vector<Vertex> const& vertices){
    std::size_t const Npoints = vertices.size();
    ANNpointArray points = annAllocPts(Npoints, 3);

    for(std::size_t i = 0; i<Npoints; ++i){
        Vector3D vertex = vertices[i]->vertex;
        
        points[i][0] = vertex.x;
        points[i][1] = vertex.y;
        points[i][2] = vertex.z;
    }
    
    return points;
}

ANNpointArray Trees::superSampleEdges(std::vector<Edge> const& edges, std::size_t const Nsample){
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937, boost::uniform_01<>> rand_gen(rng, zeroone);

    std::vector<double> start_len(edges.size(), 0.0);

    double total_len = 0;
    for(unsigned int i=0; i<edges.size(); ++i){
        Edge const& edge = edges[i];
        start_len[i] = total_len;

        double const edge_len = abs(edge->vertex2->vertex - edge->vertex1->vertex);
        total_len += edge_len;
    }

    std::cout << "\nEdge Samples : \n---------------------\n";
    ANNpointArray points = annAllocPts(Nsample, 3);

    for(std::size_t i=0; i<Nsample; ++i){
        double const sample = rand_gen()*total_len;

        auto const iter_lower_bound = std::lower_bound(start_len.begin(), start_len.end(), sample);
        std::size_t edge_index = std::distance(start_len.begin(), iter_lower_bound) - 1;

        if(std::abs(start_len[edge_index] - sample) < 1e-14){
            //! TODO: make eps a user given parameter. 
            i--;
            continue;
        }

        Edge const& edge = edges[edge_index];
        Vector3D const& edge_vec = (edge->vertex2->vertex - edge->vertex1->vertex);

        double const factor = (sample - start_len[edge_index]) / abs(edge_vec);
        Vector3D const& point = edge->vertex1->vertex + factor*edge_vec;

        points[i][0] = point.x;
        points[i][1] = point.y;
        points[i][2] = point.z;        
    }

    return points;
}

ANNpointArray Trees::superSampleFaces(std::vector<Face> const& faces, std::size_t const Nsample){
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937, boost::uniform_01<>> rand_gen(rng, zeroone);

    std::vector<double> start_area(faces.size(), 0.0);

    double total_area;
    for(unsigned int i=0; i<faces.size(); ++i){
        Face const& face = faces[i];
        start_area[i] = total_area;

        double const face_area = face->calcArea();
        total_area += face_area;
    }

    std::cout << "\nFace Samples: \n---------------------\n";
    ANNpointArray points = annAllocPts(Nsample, 3);

    for(std::size_t i = 0; i<Nsample; ++i){
        double const sample_area = rand_gen()*total_area;

        auto const iter_lower_bound = std::lower_bound(start_area.begin(), start_area.end(), sample_area);
        std::size_t face_index = std::distance(start_area.begin(), iter_lower_bound) - 1;

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

        Vector3D point = (1.0-sqrt_r1)*A + (sqrt_r1*(1.0-r2))*B + (r2*sqrt_r1)*C;

        points[i][0] = point.x;
        points[i][1] = point.y;
        points[i][2] = point.z;

    }

    return points;
}



