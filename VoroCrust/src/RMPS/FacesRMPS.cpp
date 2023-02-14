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

void FacesRMPS::divideEligbleFaces() {
    //! ASSUMPTION: This fundtion assumes that every face is a triangle!
    std::size_t const eligble_faces_size = eligble_faces.size();

    std::vector<EligbleFace> new_eligble_faces(4*eligble_faces_size, EligbleFace(std::vector<Vector3D>(3, Vector3D(0, 0, 0)), 0.0, 0.0, 0.0));

    for(std::size_t i=0; i < eligble_faces_size; ++i){
        EligbleFace const& face = eligble_faces[i];
        Vector3D const& midpoint_01 = 0.5*(face[0] + face[1]);
        Vector3D const& midpoint_12 = 0.5*(face[1] + face[2]);
        Vector3D const& midpoint_20 = 0.5*(face[2] + face[0]);

        double const sub_triangle_area = 0.25*face.area;

        new_eligble_faces[4*i][0] = face[0];
        new_eligble_faces[4*i][1] = midpoint_01;
        new_eligble_faces[4*i][2] = midpoint_20;

        new_eligble_faces[4*i+1][0] = face[1];
        new_eligble_faces[4*i+1][1] = midpoint_12;
        new_eligble_faces[4*i+1][2] = midpoint_01;

        new_eligble_faces[4*i+2][0] = face[2];
        new_eligble_faces[4*i+2][1] = midpoint_20;
        new_eligble_faces[4*i+2][2] = midpoint_12;

        new_eligble_faces[4*i+3][0] = midpoint_01;
        new_eligble_faces[4*i+3][1] = midpoint_12;
        new_eligble_faces[4*i+3][2] = midpoint_20;
        
        //! MAYBEINSERTTOLOOP:
        new_eligble_faces[4*i].patch_index = face.patch_index;
        new_eligble_faces[4*i+1].patch_index = face.patch_index;
        new_eligble_faces[4*i+2].patch_index = face.patch_index;
        new_eligble_faces[4*i+3].patch_index = face.patch_index;

        new_eligble_faces[4*i]  .plc_index = face.plc_index;
        new_eligble_faces[4*i+1].plc_index = face.plc_index;
        new_eligble_faces[4*i+2].plc_index = face.plc_index;
        new_eligble_faces[4*i+3].plc_index = face.plc_index;

        new_eligble_faces[4*i].area = sub_triangle_area;
        new_eligble_faces[4*i+1].area = sub_triangle_area;
        new_eligble_faces[4*i+2].area = sub_triangle_area;
        new_eligble_faces[4*i+3].area = sub_triangle_area;
    }

    eligble_faces = new_eligble_faces;
}