#include "PLC/PL_Complex.hpp"
#include "VoroCrustAlgorithm.hpp"
#include "vorocrust_vtk.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <array>
#include <fenv.h>

std::vector<Vector3D> read_vertices(std::string filename);
std::vector<std::vector<unsigned int>> read_faces(std::string filename);
void write_points(std::string filename, std::vector<Vector3D> points);

void triangle(){
    std::vector<Vector3D> vertices{ Vector3D(0, 0, 0), 
                                    Vector3D(3, 0, 0), 
                                    Vector3D(0, 3, 0)};

    PL_Complex plc_triangle = PL_Complex(vertices);
    plc_triangle.addFace(std::vector<unsigned int>{0,1,2});
    // std::cout << plc_triangle.repr() << std::endl;

    VoroCrustAlgorithm alg_triangle(plc_triangle, M_PI*0.1, M_PI*0.1, 0.1, 0.3, 0.13);
    
    alg_triangle.run();

    std::string dirname = "./triangle";
    std::filesystem::create_directories(dirname);
    vorocrust_vtk::write_vtu_PL_Complex(dirname+"/triangle.vtu", *alg_triangle.plc);
    vorocrust_vtk::write_vtu_trees(dirname+"/triangle_trees.vtu", alg_triangle.trees);
    vorocrust_vtk::write_ballTree(dirname+"/triangle_sharp_corners_sampling.vtp", alg_triangle.trees.ball_kd_vertices);
    vorocrust_vtk::write_ballTree(dirname+"/triangle_sharp_edges_sampling.vtp", alg_triangle.trees.ball_kd_edges);
    vorocrust_vtk::write_ballTree(dirname+"/triangle_faces_sampling.vtp", alg_triangle.trees.ball_kd_faces);

    std::vector<Vector3D> const& seeds = alg_triangle.getSeeds();

    vorocrust_vtk::write_points(dirname+"/triangle_seeds.vtu", seeds);

    write_points(dirname+"/seeds", seeds);

    // std::cout << "\n\nFINISH PART ONE\n" << std::endl;
}

void fox(){
    // std::cout << "\nRead From File\n------------------------\n\n" << std::endl;
    auto vertices_from_file = read_vertices("data/fox/vertices.txt");

    PL_Complex plc_from_file(vertices_from_file);

    for(auto& vertex : plc_from_file.vertices){
        // // std::cout << "vertex " << vertex->index << ": " << vertex->repr() << "\n";
    }

    auto faces_from_file = read_faces("data/fox/faces.txt");

    int i = 0;
    for(auto& face_indices : faces_from_file){
        i++;
        // std::cout << "face " << i << ": v1 = " << face_indices[0] << ", v2 = " << face_indices[1] << ", v3 = " << face_indices[2] << "\n";
        plc_from_file.addFace(face_indices);
    }


    // std::cout << std::endl;

    VoroCrustAlgorithm alg_fox(plc_from_file, M_PI*0.1, M_PI*0.1, 100., 0.3, 0.13);
    
    // // std::cout << alg_fox.repr() << std::endl;

    alg_fox.run();
    

    
    Vector3D query(100, 250, 50);

    int search_NN = alg_fox.trees.VC_kd_sharp_corners.nearestNeighbor(query);
    int search_kNN = alg_fox.trees.VC_kd_sharp_corners.kNearestNeighbors(query, 1)[0];

    if(search_kNN != search_NN){
        // std::cout << "\nNN != kNN when k=1\n" << std::endl;

        exit(1);
    }


    std::vector<Vector3D> centeroids(alg_fox.plc->faces.size(), {0, 0, 0});
    std::vector<Vector3D> normals(alg_fox.plc->faces.size(), {0, 0, 0});

    for(std::size_t i=0; i<alg_fox.plc->faces.size(); ++i){
        centeroids[i] = alg_fox.plc->faces[i]->calculateCenteroid();
        normals[i] = -1*alg_fox.plc->faces[i]->calcNormal();
    }

    
    std::vector<Vector3D> vertex1(alg_fox.plc->sharp_edges.size(), {0, 0, 0});
    std::vector<Vector3D> edge_vectors(alg_fox.plc->sharp_edges.size(), {0, 0, 0});

    for(std::size_t i=0; i<alg_fox.plc->sharp_edges.size(); ++i){
        vertex1[i] = alg_fox.plc->sharp_edges[i]->vertex1->vertex;
        edge_vectors[i] = alg_fox.plc->sharp_edges[i]->vertex2->vertex - alg_fox.plc->sharp_edges[i]->vertex1->vertex;
    }

    std::string dirname = "./fox";
    std::filesystem::create_directories(dirname);
    vorocrust_vtk::write_vtu_PL_Complex(dirname+"/fox.vtu", *alg_fox.plc);
    
    vorocrust_vtk::write_arbitrary_oriented_vectors(dirname+"/fox_face_noramls.vtp", centeroids,  normals, "normals", 10.0);
    
    vorocrust_vtk::write_arbitrary_oriented_vectors(dirname+"/fox_face_creases.vtp", vertex1,  edge_vectors, "creases", 1.0);

    vorocrust_vtk::write_vtu_trees(dirname+"/fox_trees.vtu", alg_fox.trees);
    
    vorocrust_vtk::write_nearestNeighbor(dirname+"/fox_nearest_vertices.vtu", alg_fox.trees.VC_kd_sharp_corners, query);
    vorocrust_vtk::write_kNearestNeighbors(dirname+"/fox_k_nearest_vertices.vtu", alg_fox.trees.VC_kd_sharp_corners, query, 25);
    vorocrust_vtk::write_radiusSearch(dirname+"/fox_radius_10.vtu", alg_fox.trees.VC_kd_faces, query, 10.0);
    vorocrust_vtk::write_radiusSearch("fox_radius_20.vtu", alg_fox.trees.VC_kd_faces, query, 20.0);

    vorocrust_vtk::write_nearestNeighborToSegment(dirname+"/fox_nearest_to_segment.vtu", alg_fox.trees.VC_kd_sharp_corners, {Vector3D(-200, 250, 0), Vector3D(200, 250, 0)}, 1e3);

    vorocrust_vtk::write_nearestNeighbor(dirname+"/fox_nearest_edges.vtu", alg_fox.trees.VC_kd_sharp_edges, query);
    vorocrust_vtk::write_nearestNeighbor(dirname+"/fox_nearest_faces.vtu", alg_fox.trees.VC_kd_faces, query);

    vorocrust_vtk::write_ballTree(dirname+"/fox_sharp_corners_sampling.vtp", alg_fox.trees.ball_kd_vertices);
    vorocrust_vtk::write_ballTree(dirname+"/fox_sharp_edges_sampling.vtp", alg_fox.trees.ball_kd_edges);
    vorocrust_vtk::write_ballTree(dirname+"/fox_sharp_faces_sampling.vtp", alg_fox.trees.ball_kd_faces);

    std::vector<std::vector<Vector3D>> leftover_eligble(alg_fox.facesDriver.eligble_faces.size(), std::vector<Vector3D>());

    for(std::size_t i=0 ;i<alg_fox.facesDriver.eligble_faces.size(); ++i){
        leftover_eligble[i] = alg_fox.facesDriver.eligble_faces[i].face;
    }

    vorocrust_vtk::write_vtu_faces(dirname+"/fox_leftover_eligble_faces.vtu", leftover_eligble);

    std::vector<Vector3D> seeds = alg_fox.getSeeds();
    vorocrust_vtk::write_points(dirname+"/fox_seeds.vtu", seeds);
    
    write_points(dirname+"/fox_seeds", seeds);

    // std::cout << "\n\nFINISH PART TWO\n\nCLICK TO END " << std::endl;
}

void astroid(){
    // std::cout << "\nRead From File\n------------------------\n\n" << std::endl;
    auto vertices_from_file = read_vertices("data/astroid/vertices.txt");

    PL_Complex plc_from_file(vertices_from_file);

    for(auto& vertex : plc_from_file.vertices){
        // // std::cout << "vertex " << vertex->index << ": " << vertex->repr() << "\n";
    }

    auto faces_from_file = read_faces("data/astroid/faces.txt");

    int i = 0;
    for(auto& face_indices : faces_from_file){
        i++;
        // // std::cout << "face " << i << ": v1 = " << face_indices[0] << ", v2 = " << face_indices[1] << ", v3 = " << face_indices[2] << "\n";
        plc_from_file.addFace(face_indices);
    }


    // std::cout << std::endl;

    VoroCrustAlgorithm alg_astroid(plc_from_file, M_PI*30./180., M_PI*30./180., 100., 0.3, 0.13);
    
    // // std::cout << alg_astroid.repr() << std::endl;

    alg_astroid.run();
    

    
    Vector3D query(100, 250, 50);

    int search_NN = alg_astroid.trees.VC_kd_sharp_corners.nearestNeighbor(query);
    int search_kNN = alg_astroid.trees.VC_kd_sharp_corners.kNearestNeighbors(query, 1)[0];

    if(search_kNN != search_NN){
        // std::cout << "\nNN != kNN when k=1\n" << std::endl;

        exit(1);
    }


    std::vector<Vector3D> centeroids(alg_astroid.plc->faces.size(), {0, 0, 0});
    std::vector<Vector3D> normals(alg_astroid.plc->faces.size(), {0, 0, 0});

    for(std::size_t i=0; i<alg_astroid.plc->faces.size(); ++i){
        centeroids[i] = alg_astroid.plc->faces[i]->calculateCenteroid();
        normals[i] = -1*alg_astroid.plc->faces[i]->calcNormal();
    }

    
    std::vector<Vector3D> vertex1(alg_astroid.plc->sharp_edges.size(), {0, 0, 0});
    std::vector<Vector3D> edge_vectors(alg_astroid.plc->sharp_edges.size(), {0, 0, 0});

    for(std::size_t i=0; i<alg_astroid.plc->sharp_edges.size(); ++i){
        vertex1[i] = alg_astroid.plc->sharp_edges[i]->vertex1->vertex;
        edge_vectors[i] = alg_astroid.plc->sharp_edges[i]->vertex2->vertex - alg_astroid.plc->sharp_edges[i]->vertex1->vertex;
    }

    std::string dirname = "./astroid";
    std::filesystem::create_directories(dirname);
    vorocrust_vtk::write_vtu_PL_Complex(dirname+"/astroid.vtu", *alg_astroid.plc);
    
    vorocrust_vtk::write_arbitrary_oriented_vectors(dirname+"/astroid_face_noramls.vtp", centeroids,  normals, "normals", 10.0);
    
    vorocrust_vtk::write_arbitrary_oriented_vectors(dirname+"/astroid_face_creases.vtp", vertex1,  edge_vectors, "creases", 1.0);

    vorocrust_vtk::write_vtu_trees(dirname+"/astroid_trees.vtu", alg_astroid.trees);
    
    vorocrust_vtk::write_nearestNeighbor(dirname+"/astroid_nearest_vertices.vtu", alg_astroid.trees.VC_kd_sharp_corners, query);
    vorocrust_vtk::write_kNearestNeighbors(dirname+"/astroid_k_nearest_vertices.vtu", alg_astroid.trees.VC_kd_sharp_corners, query, 25);
    vorocrust_vtk::write_radiusSearch(dirname+"/astroid_radius_10.vtu", alg_astroid.trees.VC_kd_faces, query, 10.0);
    vorocrust_vtk::write_radiusSearch("astroid_radius_20.vtu", alg_astroid.trees.VC_kd_faces, query, 20.0);

    vorocrust_vtk::write_nearestNeighborToSegment(dirname+"/astroid_nearest_to_segment.vtu", alg_astroid.trees.VC_kd_sharp_corners, {Vector3D(-200, 250, 0), Vector3D(200, 250, 0)}, 1e3);

    vorocrust_vtk::write_nearestNeighbor(dirname+"/astroid_nearest_edges.vtu", alg_astroid.trees.VC_kd_sharp_edges, query);
    vorocrust_vtk::write_nearestNeighbor(dirname+"/astroid_nearest_faces.vtu", alg_astroid.trees.VC_kd_faces, query);

    vorocrust_vtk::write_ballTree(dirname+"/astroid_sharp_corners_sampling.vtp", alg_astroid.trees.ball_kd_vertices);
    vorocrust_vtk::write_ballTree(dirname+"/astroid_sharp_edges_sampling.vtp", alg_astroid.trees.ball_kd_edges);
    vorocrust_vtk::write_ballTree(dirname+"/astroid_sharp_faces_sampling.vtp", alg_astroid.trees.ball_kd_faces);

    std::vector<std::vector<Vector3D>> leftover_eligble(alg_astroid.facesDriver.eligble_faces.size(), std::vector<Vector3D>());

    for(std::size_t i=0 ;i<alg_astroid.facesDriver.eligble_faces.size(); ++i){
        leftover_eligble[i] = alg_astroid.facesDriver.eligble_faces[i].face;
    }

    vorocrust_vtk::write_vtu_faces(dirname+"/astroid_leftover_eligble_faces.vtu", leftover_eligble);

    std::vector<Vector3D> seeds = alg_astroid.getSeeds();
    vorocrust_vtk::write_points(dirname+"/astroid_seeds.vtu", seeds);
    
    write_points(dirname+"/astroid_seeds", seeds);

    std::cout << "\n\nFINISH PART TWO\n\nCLICK TO END " << std::endl;
}


std::vector<Vector3D> load_seeds(std::string filename){
    std::vector<Vector3D> seeds;
    
    std::fstream myfile_x(filename + "_x.txt", std::ios_base::in);
    std::fstream myfile_y(filename + "_y.txt", std::ios_base::in);
    std::fstream myfile_z(filename + "_z.txt", std::ios_base::in);
    
    double x;
    double y;
    double z;
    while(myfile_x >> x){
        myfile_y >> y;
        myfile_z >> z;

        seeds.push_back(Vector3D(x, y, z));
        // // std::cout << "seed x: " << x << ", y: " << y << ", z: " << z << std::endl;
    }

    return seeds;
}

void box(){
    std::vector<Vector3D> vertices{
        Vector3D(0, 0, 0), // 0
        Vector3D(3, 0, 0), // 1
        Vector3D(0, 3, 0), // 2
        Vector3D(0, 0, 3), // 3
        Vector3D(3, 3, 0), // 4
        Vector3D(0, 3, 3), // 5
        Vector3D(3, 0, 3), // 6
        Vector3D(3, 3, 3), // 7
    };

    PL_Complex plc_box = PL_Complex(vertices);

    // square 1
    plc_box.addFace({0, 1, 4});
    plc_box.addFace({4, 2, 0});

    // square 2
    plc_box.addFace({0, 1, 6});
    plc_box.addFace({6, 3, 0});

    // square 3
    plc_box.addFace({0, 3, 5});
    plc_box.addFace({5, 2, 0});
    
    // square 4
    plc_box.addFace({4, 1, 6});
    plc_box.addFace({4, 6, 7});

    // square 5
    plc_box.addFace({4, 2, 5});
    plc_box.addFace({4, 5, 7});

    // square 6
    plc_box.addFace({3, 6, 7});
    plc_box.addFace({3, 5, 7});
    
    plc_box.detectFeatures(M_PI*0.1, M_PI*0.1);

    std::vector<Vector3D> seeds{{1.5, 1.3, 1.}};

    std::string dirname = "./box";
    std::filesystem::create_directories(dirname);
    
    vorocrust_vtk::write_vtu_PL_Complex(dirname+"/box.vtu", plc_box);
    vorocrust_vtk::write_points(dirname+"/seeds.vtu", seeds);

    VoroCrustAlgorithm alg_box(plc_box, M_PI*0.1, M_PI*0.1, 0.3, 0.3, 0.13);
    auto const& [in_seeds, out_seeds] =  alg_box.determineIfSeedsAreInsideOrOutside(seeds);
    vorocrust_vtk::write_points(dirname+"/box_seeds_in.vtu", in_seeds);
    vorocrust_vtk::write_points(dirname+"/box_seeds_out.vtu", out_seeds);
}

void loadfox_and_split_seeds(){
    // std::cout << "\nRead From File\n------------------------\n\n" << std::endl;
    auto vertices_from_file = read_vertices("data/fox/vertices.txt");

    PL_Complex plc_from_file(vertices_from_file);

    auto faces_from_file = read_faces("data/fox/faces.txt");

    int i = 0;
    for(auto& face_indices : faces_from_file){
        i++;
        plc_from_file.addFace(face_indices);
    }


    // std::cout << std::endl;

    VoroCrustAlgorithm alg_fox(plc_from_file, M_PI*0.1, M_PI*0.1, 100., 0.3, 0.13);

    auto const& seeds = load_seeds("fox/fox_seeds");

    auto const& [in_seeds, out_seeds] =  alg_fox.determineIfSeedsAreInsideOrOutside(seeds);
    write_points("fox/fox_in_seeds", in_seeds);
    write_points("fox/fox_out_seeds", out_seeds);

    vorocrust_vtk::write_points("fox/fox_seeds_in.vtu", in_seeds);
    vorocrust_vtk::write_points("fox/fox_seeds_out.vtu", out_seeds);


}

void loadastroid_and_split_seeds(){
    // std::cout << "\nRead From File\n------------------------\n\n" << std::endl;
    auto vertices_from_file = read_vertices("data/astroid/vertices.txt");

    PL_Complex plc_from_file(vertices_from_file);

    auto faces_from_file = read_faces("data/astroid/faces.txt");

    int i = 0;
    for(auto& face_indices : faces_from_file){
        i++;
        plc_from_file.addFace(face_indices);
    }


    // std::cout << std::endl;

    VoroCrustAlgorithm alg_astroid(plc_from_file, M_PI*30./180., M_PI*30./180., 100., 0.3, 0.13);

    auto const& seeds = load_seeds("astroid/astroid_seeds");

    auto const& [in_seeds, out_seeds] =  alg_astroid.determineIfSeedsAreInsideOrOutside(seeds);
    write_points("astroid/astroid_in_seeds", in_seeds);
    write_points("astroid/astroid_out_seeds", out_seeds);

    vorocrust_vtk::write_points("astroid/astroid_seeds_in.vtu", in_seeds);
    vorocrust_vtk::write_points("astroid/astroid_seeds_out.vtu", out_seeds);


}

int main(int argc, char *argv[]){
	// feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    // triangle();
    // getchar();
    // fox();
    // getchar();
    // loadfox_and_split_seeds();

    astroid();
    loadastroid_and_split_seeds();

    // box();

    
    return 0;

}

std::vector<Vector3D> read_vertices(std::string filename){
    std::fstream myfile(filename, std::ios_base::in);
    
    auto number_of_lines = std::count(std::istreambuf_iterator<char>(myfile), std::istreambuf_iterator<char>(), '\n');
    // std::cout << "Reading vertices data from file: " << filename << ", number of lines (vertices): " << number_of_lines << std::endl;
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
    // std::cout << "Reading vertices data from file: " << filename << ", number of lines (vertices): " << number_of_lines << "\n";

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

void write_points(std::string filename, std::vector<Vector3D> points){
    std::ofstream myfile_x(filename + "_x.txt");
    std::ofstream myfile_y(filename + "_y.txt");
    std::ofstream myfile_z(filename + "_z.txt");

    for(Vector3D const& p : points){
        myfile_x << std::setprecision(8) << p.x << "\n";
        myfile_y << std::setprecision(8) << p.y << "\n";
        myfile_z << std::setprecision(8) << p.z << "\n";
    }

    myfile_x.close();
    myfile_y.close();
    myfile_z.close();
}



