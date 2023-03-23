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

std::vector<Seed> unsorted_unique(std::vector<Seed> const& vec, double const tol);
std::vector<Vector3D> read_vertices(std::string filename);
std::vector<std::vector<std::size_t>> read_faces(std::string filename);
void write_seeds(std::string filename, std::vector<Seed> const& points);
void write_points(std::string filename, std::vector<Vector3D> const& points);
std::vector<Seed> load_seeds(std::string filename);

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

    VoroCrustAlgorithm alg_astroid(plc_from_file, M_PI*75./180., M_PI*75./180., 0.01, 0.3, 0.13);
    
    // // std::cout << alg_astroid.repr() << std::endl;

    alg_astroid.run();
    

    

    std::vector<Vector3D> centeroids(alg_astroid.plc->faces.size(), {0, 0, 0});
    std::vector<Vector3D> normals(alg_astroid.plc->faces.size(), {0, 0, 0});

    std::string dirname = "./astroid";
    std::filesystem::create_directories(dirname);
    vorocrust_vtk::write_vtu_PL_Complex(dirname+"/astroid.vtu", *alg_astroid.plc);

    vorocrust_vtk::write_vtu_trees(dirname+"/astroid_trees.vtu", alg_astroid.trees);

    vorocrust_vtk::write_ballTree(dirname+"/astroid_sharp_corners_sampling.vtp", alg_astroid.trees.ball_kd_vertices);
    vorocrust_vtk::write_ballTree(dirname+"/astroid_sharp_edges_sampling.vtp", alg_astroid.trees.ball_kd_edges);
    vorocrust_vtk::write_ballTree(dirname+"/astroid_sharp_faces_sampling.vtp", alg_astroid.trees.ball_kd_faces);

    auto const& seeds = alg_astroid.getSeeds();

    std::vector<Vector3D> seeds_points;
    seeds_points.reserve(seeds.size() + 100);

    for(auto const& seed : seeds)
        seeds_points.push_back(seed.p);

    vorocrust_vtk::write_points(dirname+"/astroid_seeds.vtu", seeds_points);
    
    write_seeds(dirname+"/astroid_seeds", seeds);

    std::cout << "\n\nFINISH PART TWO\n\nCLICK TO END " << std::endl;
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

    VoroCrustAlgorithm alg_astroid(plc_from_file, M_PI*75./180., M_PI*75./180., 0.01, 0.3, 0.13);

    auto const& seeds = load_seeds("astroid/astroid_seeds");

    // auto const& [in_seeds, out_seeds] =  alg_astroid.calcVolumeSeedsUniform(seeds, 50, 50, 50);
    auto const& [in_seeds, out_seeds] =  alg_astroid.calcVolumeSeedsNonUniform(seeds, 0.00125);

    write_points("astroid/astroid_in_seeds", in_seeds);
    write_points("astroid/astroid_out_seeds", out_seeds);

    vorocrust_vtk::write_points("astroid/astroid_seeds_in.vtu", in_seeds);
    vorocrust_vtk::write_points("astroid/astroid_seeds_out.vtu", out_seeds);


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
        // // std::cout << "face " << i << ": v1 = " << face_indices[0] << ", v2 = " << face_indices[1] << ", v3 = " << face_indices[2] << "\n";
        plc_from_file.addFace(face_indices);
    }


    // std::cout << std::endl;

    VoroCrustAlgorithm alg_fox(plc_from_file, M_PI*10./180., M_PI*10./180., 100.0, 0.3, 0.13);
    
    // // std::cout << alg_fox.repr() << std::endl;

    alg_fox.run();
    

    

    std::vector<Vector3D> centeroids(alg_fox.plc->faces.size(), {0, 0, 0});
    std::vector<Vector3D> normals(alg_fox.plc->faces.size(), {0, 0, 0});

    std::string dirname = "./fox";
    std::filesystem::create_directories(dirname);
    vorocrust_vtk::write_vtu_PL_Complex(dirname+"/fox.vtu", *alg_fox.plc);

    vorocrust_vtk::write_vtu_trees(dirname+"/fox_trees.vtu", alg_fox.trees);

    vorocrust_vtk::write_ballTree(dirname+"/fox_sharp_corners_sampling.vtp", alg_fox.trees.ball_kd_vertices);
    vorocrust_vtk::write_ballTree(dirname+"/fox_sharp_edges_sampling.vtp", alg_fox.trees.ball_kd_edges);
    vorocrust_vtk::write_ballTree(dirname+"/fox_sharp_faces_sampling.vtp", alg_fox.trees.ball_kd_faces);

    auto const& seeds = alg_fox.getSeeds();

    std::vector<Vector3D> seeds_points;
    seeds_points.reserve(seeds.size() + 100);

    for(auto const& seed : seeds)
        seeds_points.push_back(seed.p);

    vorocrust_vtk::write_points(dirname+"/fox_seeds.vtu", seeds_points);
    
    write_seeds(dirname+"/fox_seeds", seeds);

    auto turnToVector3D = [](auto&& seed_vec){
        std::vector<Vector3D> vec;
        vec.reserve(seed_vec.size()*2);
        for(auto const& seed : seed_vec) vec.push_back(std::move(seed.p));
        return vec;
    };

    // auto const& [in_seeds, out_seeds] = alg_fox.determineIfSeedsAreInsideOrOutside(seeds);
    // write_points("fox/fox_in_seeds", turnToVector3D(in_seeds));
    // write_points("fox/fox_out_seeds", turnToVector3D(out_seeds));
    
    // auto const& [in_seeds, out_seeds] =  alg_fox.calcVolumeSeedsUniform(seeds, 100, 100, 100);
    auto const& [in_seeds, out_seeds] =  alg_fox.calcVolumeSeedsNonUniform(seeds, 2);

    write_points("fox/fox_in_seeds", in_seeds);
    write_points("fox/fox_out_seeds", out_seeds);

    vorocrust_vtk::write_points("fox/fox_seeds_in.vtu", in_seeds);
    vorocrust_vtk::write_points("fox/fox_seeds_out.vtu", out_seeds);
    std::cout << "\n\nFINISH PART TWO\n\nCLICK TO END " << std::endl;
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

    VoroCrustAlgorithm alg_fox(plc_from_file, M_PI*30./180., M_PI*30./180., 1, 0.3, 0.13);

    auto const& seeds = load_seeds("fox/fox_seeds");


    auto turnToVector3D = [](auto&& seed_vec){
        std::vector<Vector3D> vec;
        vec.reserve(seed_vec.size()*2);
        for(auto const& seed : seed_vec) vec.push_back(std::move(seed.p));
        return vec;
    };

    // auto const& [in_seeds, out_seeds] = alg_fox.determineIfSeedsAreInsideOrOutside(seeds);
    // write_points("fox/fox_in_seeds", turnToVector3D(in_seeds));
    // write_points("fox/fox_out_seeds", turnToVector3D(out_seeds));
    
    // auto const& [in_seeds, out_seeds] =  alg_fox.calcVolumeSeedsUniform(seeds, 100, 100, 100);
    auto const& [in_seeds, out_seeds] =  alg_fox.calcVolumeSeedsNonUniform(seeds, 1);

    write_points("fox/fox_in_seeds", in_seeds);
    write_points("fox/fox_out_seeds", out_seeds);

    vorocrust_vtk::write_points("fox/fox_seeds_in.vtu", in_seeds);
    vorocrust_vtk::write_points("fox/fox_seeds_out.vtu", out_seeds);

}

std::vector<Seed> load_seeds(std::string filename){
    std::vector<Seed> seeds;
    
    std::fstream myfile_x(filename + "_x.txt", std::ios_base::in);
    std::fstream myfile_y(filename + "_y.txt", std::ios_base::in);
    std::fstream myfile_z(filename + "_z.txt", std::ios_base::in);
    std::fstream myfile_r(filename + "_r.txt", std::ios_base::in);

    
    double x;
    double y;
    double z;
    double r;
    while(myfile_x >> x){
        myfile_y >> y;
        myfile_z >> z;
        myfile_r >> r;


        seeds.push_back(Seed(Vector3D(x, y, z), r));
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
    // vorocrust_vtk::write_seeds(dirname+"/seeds.vtu", seeds);

    VoroCrustAlgorithm alg_box(plc_box, M_PI*0.1, M_PI*0.1, 0.3, 0.3, 0.13);
    // auto const& [in_seeds, out_seeds] =  alg_box.determineIfSeedsAreInsideOrOutside(seeds);
    // vorocrust_vtk::write_points(dirname+"/box_seeds_in.vtu", in_seeds);
    // vorocrust_vtk::write_points(dirname+"/box_seeds_out.vtu", out_seeds);
}

int main(int argc, char *argv[]){
    std::cout << "RUNNING MAIN!!" << std::endl;
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    // triangle();
    // getchar();
    fox();
    // loadfox_and_split_seeds();

    // astroid();
    // loadastroid_and_split_seeds();

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
    std::size_t i = 0;
    while (myfile >> x)
    {   
        myfile >> y;
        myfile >> z;

        vertices_from_file[i].Set(x, y, z);
        i++;
    }

    return vertices_from_file;
}

std::vector<std::vector<std::size_t>> read_faces(std::string filename){
    std::fstream myfile(filename, std::ios_base::in);

    auto number_of_lines = std::count(std::istreambuf_iterator<char>(myfile), std::istreambuf_iterator<char>(), '\n');
    // std::cout << "Reading vertices data from file: " << filename << ", number of lines (vertices): " << number_of_lines << "\n";

    myfile.clear();
    myfile.seekg(0);

    std::vector<std::vector<std::size_t>> faces_from_file = std::vector<std::vector<std::size_t>>(number_of_lines, std::vector<std::size_t>({0, 0, 0}));
    
    double v1, v2, v3;
    
    std::size_t i=0;
    while(myfile >> v1){

        myfile >> v2;
        myfile >> v3;

        faces_from_file[i][0]=static_cast<std::size_t>(v1);
        faces_from_file[i][1]=static_cast<std::size_t>(v2);
        faces_from_file[i][2]=static_cast<std::size_t>(v3);
        i++;
    }

    return faces_from_file;
}

void write_points(std::string filename, std::vector<Vector3D> const& points){
    std::ofstream myfile_x(filename + "_x.txt");
    std::ofstream myfile_y(filename + "_y.txt");
    std::ofstream myfile_z(filename + "_z.txt");

    for(auto const& p : points){
        myfile_x << std::setprecision(15) << p.x << "\n";
        myfile_y << std::setprecision(15) << p.y << "\n";
        myfile_z << std::setprecision(15) << p.z << "\n";
    }

    myfile_x.close();
    myfile_y.close();
    myfile_z.close();
}


void write_seeds(std::string filename, std::vector<Seed> const& seeds){
    std::ofstream myfile_x(filename + "_x.txt");
    std::ofstream myfile_y(filename + "_y.txt");
    std::ofstream myfile_z(filename + "_z.txt");
    std::ofstream myfile_r(filename + "_r.txt");

    for(auto const& seed : seeds){
        myfile_x << std::setprecision(15) << seed.p.x << "\n";
        myfile_y << std::setprecision(15) << seed.p.y << "\n";
        myfile_z << std::setprecision(15) << seed.p.z << "\n";
        myfile_r << std::setprecision(15) << seed.radius << "\n";
    }

    myfile_x.close();
    myfile_y.close();
    myfile_z.close();
    myfile_r.close();
}



