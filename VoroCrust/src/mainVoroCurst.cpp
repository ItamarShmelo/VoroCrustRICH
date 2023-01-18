#include "PL_Complex.hpp"
#include "VoroCrustAlgorithm.hpp"
#include "vorocrust_vtk.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

std::vector<Vector3D> read_vertices(std::string filename);
std::vector<std::vector<unsigned int>> read_faces(std::string filename);

int main(int argc, char *argv[]){
    
    std::vector<Vector3D> vertices{ Vector3D(0, 0, 0), 
                                    Vector3D(1, 0, 0), 
                                    Vector3D(1, 1, 0), 
                                    Vector3D(0, 1, 0), 
                                    Vector3D(0, 2, 0),
                                    Vector3D(1, 2, 0)};

    PL_Complex plc = PL_Complex(vertices);
    plc.addFace(std::vector<unsigned int>{0,1,2,3});
    plc.addFace(std::vector<unsigned int>{2,3,4,5});

    std::cout << plc.repr() << std::endl;

    std::cout << "\n\nchange one vertex" << std::endl;
    plc.vertices[2]->vertex.x = 2;
    std::cout << plc.repr() << std::endl;

    std::cout << "\n\nCheck VoroCrustAlgorithm " << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    plc.vertices[2]->vertex.x = 1;
    VoroCrustAlgorithm alg(plc, M_PI/10.0, M_PI_4, 1., 0.8);

    std::cout << alg.repr() << std::endl;

    std::cout << "\n\nWrite VTK File for PLC\n-------------------------" << std::endl;
    
    vorocrust_vtk::write_vtu_PL_Complex("plc.vtu", plc);
    alg.run();
    
    std::cout << "\n\nFINISH PART ONE\n" << std::endl;

    std::cout << "\nRead From File\n------------------------\n\n" << std::endl;
    auto vertices_from_file = read_vertices("data/vertices.txt");

    PL_Complex plc_from_file(vertices_from_file);

    for(auto& vertex : plc_from_file.vertices){
        std::cout << "vertex " << vertex->index << ": " << vertex->repr() << "\n";
    }

    auto faces_from_file = read_faces("data/faces.txt");

    int i = 0;
    for(auto& face_indices : faces_from_file){
        i++;
        std::cout << "face " << i << ": v1 = " << face_indices[0] << ", v2 = " << face_indices[1] << ", v3 = " << face_indices[2] << "\n";
        plc_from_file.addFace(face_indices);
    }

    std::cout << plc_from_file.repr() << std::endl;

    std::cout << std::endl;

    VoroCrustAlgorithm alg_cat(plc_from_file, M_PI*0.1, M_PI*0.1, 1., 0.8);
    
    std::cout << alg_cat.repr() << std::endl;

    alg_cat.run();
    
    vorocrust_vtk::write_vtu_PL_Complex("cat.vtu", alg_cat.plc);

    std::cout << "\n\nFINISH PART TWO\n" << std::endl;

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



