#include "FacesRMPS.hpp"

FacesRMPS::FacesRMPS(double const maxRadius_, 
                     double const L_Lipschitz_, 
                     double const alpha_, 
                     double const sharpTheta_, 
                     std::shared_ptr<PL_Complex> const& plc_) : 
                                    maxRadius(maxRadius_), 
                                    L_Lipschitz(L_Lipschitz_), 
                                    alpha(alpha_), 
                                    sharpTheta(sharpTheta_), 
                                    uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()), 
                                    plc(plc_), 
                                    eligble_faces() {}

void FacesRMPS::loadFaces(std::vector<Face> const& faces) {
    if(not eligble_faces.empty()){
        std::cout << "eligble faces is not empty when loading faces" << std::endl;
        exit(1);
    }

    for(Face const& face : faces){
        std::vector<Vector3D> vertices;
        for(Vertex const& vertex : face->vertices) {
            vertices.push_back(vertex->vertex);
        }

        eligble_faces.push_back(EligbleFace(vertices, face->patch_index, face->index, face->calcArea()));
    }
}

std::pair<double const, std::vector<double> const> FacesRMPS::calculateTotalAreaAndStartAreaOfEligbleFaces() const {
    std::size_t const num_of_eligble_faces = eligble_faces.size();

    std::vector<double> start_area(num_of_eligble_faces, 0.0);
    double total_area = 0.0;

    for(std::size_t i=0; i<num_of_eligble_faces; ++i){
        EligbleFace const& face = eligble_faces[i];
        start_area[i] = total_area;
        total_area += face.area;
    }

    return std::pair<double const, std::vector<double> const>(total_area, start_area);
}