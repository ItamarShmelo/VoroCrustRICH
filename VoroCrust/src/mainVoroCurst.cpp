#include "PLC/PL_Complex.hpp"
#include "VoroCrustAlgorithm.hpp"
#include "vorocrust_vtk.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <array>

std::vector<Vector3D> read_vertices(std::string filename);
std::vector<std::vector<unsigned int>> read_faces(std::string filename);

int main(int argc, char *argv[]){
    // std::vector<Vector3D> vertices{ Vector3D(0, 0, 0), 
    //                                 Vector3D(1, 0, 0), 
    //                                 Vector3D(1, 1, 0)};

    // PL_Complex plc_triangle = PL_Complex(vertices);
    // plc_triangle.addFace(std::vector<unsigned int>{0,1,2});
    // std::cout << plc_triangle.repr() << std::endl;

    // VoroCrustAlgorithm alg_triangle(plc_triangle, M_PI*0.1, M_PI*0.1, 0.1, 0.3, 0.13);
    
    // alg_triangle.run();

    //vorocrust_vtk::write_vtu_PL_Complex("triangle.vtu", alg_triangle.plc);
    //vorocrust_vtk::write_vtu_trees("triangle_trees.vtu", alg_triangle.trees);
    //vorocrust_vtk::write_ballTree("triangle_sharp_corners_sampling.vtp", alg_triangle.trees.ball_kd_vertices);
    //vorocrust_vtk::write_ballTree("triangle_sharp_edges_sampling.vtp", alg_triangle.trees.ball_kd_edges);

    std::cout << "\n\nFINISH PART ONE\n" << std::endl;
    
    // getchar();
    
    std::cout << "\nRead From File\n------------------------\n\n" << std::endl;
    auto vertices_from_file = read_vertices("data/fox/vertices.txt");

    PL_Complex plc_from_file(vertices_from_file);

    for(auto& vertex : plc_from_file.vertices){
        std::cout << "vertex " << vertex->index << ": " << vertex->repr() << "\n";
    }

    auto faces_from_file = read_faces("data/fox/faces.txt");

    int i = 0;
    for(auto& face_indices : faces_from_file){
        i++;
        std::cout << "face " << i << ": v1 = " << face_indices[0] << ", v2 = " << face_indices[1] << ", v3 = " << face_indices[2] << "\n";
        plc_from_file.addFace(face_indices);
    }

    std::cout << plc_from_file.repr() << std::endl;

    std::cout << std::endl;

    VoroCrustAlgorithm alg_fox(plc_from_file, M_PI*0.1, M_PI*0.1, 10., 0.3, 0.13);
    
    std::cout << alg_fox.repr() << std::endl;

    alg_fox.run();
    

    
    Vector3D query(100, 250, 50);

    int search_NN = alg_fox.trees.VC_kd_sharp_corners.nearestNeighbor(query);
    int search_kNN = alg_fox.trees.VC_kd_sharp_corners.kNearestNeighbors(query, 1)[0];

    if(search_kNN != search_NN){
        std::cout << "\nNN != kNN when k=1\n" << std::endl;

        exit(1);
    }


    std::vector<Vector3D> centeroids(alg_fox.plc.faces.size(), {0, 0, 0});
    std::vector<Vector3D> normals(alg_fox.plc.faces.size(), {0, 0, 0});

    for(std::size_t i=0; i<alg_fox.plc.faces.size(); ++i){
        centeroids[i] = alg_fox.plc.faces[i]->calculateCenteroid();
        normals[i] = -1*alg_fox.plc.faces[i]->calcNormal();
    }

    
    std::vector<Vector3D> vertex1(alg_fox.plc.sharp_edges.size(), {0, 0, 0});
    std::vector<Vector3D> edge_vectors(alg_fox.plc.sharp_edges.size(), {0, 0, 0});

    for(std::size_t i=0; i<alg_fox.plc.sharp_edges.size(); ++i){
        vertex1[i] = alg_fox.plc.sharp_edges[i]->vertex1->vertex;
        edge_vectors[i] = alg_fox.plc.sharp_edges[i]->vertex2->vertex - alg_fox.plc.sharp_edges[i]->vertex1->vertex;
    }
    // vorocrust_vtk::write_vtu_PL_Complex("fox.vtu", alg_fox.plc);
    
    // vorocrust_vtk::write_arbitrary_oriented_vectors("fox_face_noramls.vtp", centeroids,  normals, "normals", 10.0);
    
    // vorocrust_vtk::write_arbitrary_oriented_vectors("fox_face_creases.vtp", vertex1,  edge_vectors, "creases", 1.0);

    // vorocrust_vtk::write_vtu_trees("fox_trees.vtu", alg_fox.trees);
    
    // vorocrust_vtk::write_nearestNeighbor("fox_nearest_vertices.vtu", alg_fox.trees.VC_kd_sharp_corners, query);
    // vorocrust_vtk::write_kNearestNeighbors("fox_k_nearest_vertices.vtu", alg_fox.trees.VC_kd_sharp_corners, query, 25);
    // vorocrust_vtk::write_radiusSearch("fox_radius_10.vtu", alg_fox.trees.VC_kd_faces, query, 10.0);
    // vorocrust_vtk::write_radiusSearch("fox_radius_20.vtu", alg_fox.trees.VC_kd_faces, query, 20.0);

    // vorocrust_vtk::write_nearestNeighborToSegment("fox_nearest_to_segment.vtu", alg_fox.trees.VC_kd_sharp_corners, {Vector3D(-200, 250, 0), Vector3D(200, 250, 0)}, 1e3);

    // vorocrust_vtk::write_nearestNeighbor("fox_nearest_edges.vtu", alg_fox.trees.VC_kd_sharp_edges, query);
    // vorocrust_vtk::write_nearestNeighbor("fox_nearest_faces.vtu", alg_fox.trees.VC_kd_faces, query);

    // vorocrust_vtk::write_ballTree("fox_sharp_corners_sampling.vtp", alg_fox.trees.ball_kd_vertices);
    // vorocrust_vtk::write_ballTree("fox_sharp_edges_sampling.vtp", alg_fox.trees.ball_kd_edges);

    std::cout << "\n\nFINISH PART TWO\n\n CLICK TO END " << std::endl;
    getchar();
    return 0;

}

std::vector<Vector3D> read_vertices(std::string filename){
    std::fstream myfile(filename, std::ios_base::in);
    
    auto number_of_lines = std::count(std::istreambuf_iterator<char>(myfile), std::istreambuf_iterator<char>(), '\n');
    std::cout << "Reading vertices data from file: " << filename << ", number of lines (vertices): " << number_of_lines << std::endl;
    myfile.clear();
    myfile.seekg(0);
    std::vector<Vector3D> vertices_from_file = std::vector<Vector3D>(number_of_lines, Vector3D());

    double x, y, z;
    unsigned int i = 0;
    while (myfile >> x)
    {   
        myfile >> y;
        myfile >> z;

        vertices_from_file[i].Set(x, y, z);
        i++;
    }

    return vertices_from_file;
}

std::vector<std::vector<unsigned int>> read_faces(std::string filename){
    std::fstream myfile(filename, std::ios_base::in);

    auto number_of_lines = std::count(std::istreambuf_iterator<char>(myfile), std::istreambuf_iterator<char>(), '\n');
    std::cout << "Reading vertices data from file: " << filename << ", number of lines (vertices): " << number_of_lines << "\n";

    myfile.clear();
    myfile.seekg(0);

    std::vector<std::vector<unsigned int>> faces_from_file = std::vector<std::vector<unsigned int>>(number_of_lines, std::vector<unsigned int>({0, 0, 0}));
    
    double v1, v2, v3;
    
    unsigned int i=0;
    while(myfile >> v1){

        myfile >> v2;
        myfile >> v3;

        faces_from_file[i][0]=static_cast<unsigned int>(v1);
        faces_from_file[i][1]=static_cast<unsigned int>(v2);
        faces_from_file[i][2]=static_cast<unsigned int>(v3);
        i++;
    }

    return faces_from_file;
}



