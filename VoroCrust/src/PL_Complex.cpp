#include "PL_Complex.hpp"
#include <iostream>
#include <cmath>
#include "../../source/misc/utils.hpp"

PL_Complex::PL_Complex(std::vector<Vector3D> const& vertices_) : vertices(), faces(){
    std::size_t index = 0;
    for(auto& vertex : vertices_) {
        vertices.push_back(std::make_shared<VoroCrustVertex>(vertex, index));
        index++;
    }
}

std::shared_ptr<VoroCrustEdge> PL_Complex::addEdge(std::shared_ptr<VoroCrustVertex> const& v1, std::shared_ptr<VoroCrustVertex> const& v2){
    for (auto& edge : edges){
        if (edge->checkIfEqual(v1, v2)){
            return edge;
        }
    }

    auto new_edge_ptr = std::make_shared<VoroCrustEdge>(v1, v2, edges.size());
    
    new_edge_ptr->vertex1->addEdge(new_edge_ptr);
    new_edge_ptr->vertex2->addEdge(new_edge_ptr);

    edges.push_back(new_edge_ptr);
    return new_edge_ptr;
}

void PL_Complex::addFace(std::vector<unsigned int> const& indices){

    if (indices.size() != unique(indices).size()){
        std::cout << "ERROR: Repeated indices in PL_Complex::addFace!" << std::endl;
        exit(1);
    }

    if(indices.size() < 3){
        std::cout << "ERROR: Face has to have at least three vertices!" << std::endl;
        exit(1);
    }

    std::vector<std::shared_ptr<VoroCrustVertex>> face_vertices = std::vector<std::shared_ptr<VoroCrustVertex>>();

    for(auto& index : indices)
        face_vertices.push_back(vertices[index]);

    std::shared_ptr<VoroCrustFace> new_face_ptr = std::make_shared<VoroCrustFace>(face_vertices, faces.size());

    for (std::shared_ptr<VoroCrustVertex> vertex_ptr : face_vertices){
        vertex_ptr->addFace(new_face_ptr);
    }

    for (std::size_t i = 0; i < new_face_ptr->vertices.size(); ++i){
        auto new_edge_ptr = this->addEdge(new_face_ptr->vertices[i], new_face_ptr->vertices[ (i+1) % (new_face_ptr->vertices.size())]);
        new_edge_ptr->addFace(new_face_ptr);
        new_face_ptr->addEdge(new_edge_ptr);
        
        if(new_edge_ptr->faces.size() == 2){
            new_face_ptr->neighbors.push_back(new_edge_ptr->faces[0]);
            new_edge_ptr->faces[0]->neighbors.push_back(new_edge_ptr->faces[1]);

            std::cout << "Neigbors :: faces " << new_face_ptr->index << " and face " << new_edge_ptr->faces[0]->index << std::endl;
        }
    }

    faces.push_back(new_face_ptr);
}

bool PL_Complex::checkAllVerticesAreOnFace(){
    for (auto& vertex : vertices)
        if (vertex->faces.size() == 0){
            std::cout << "ERROR: vertex " << vertex->index << " is not part of a face!!";
            return false;
        }

    return true;
}

bool PL_Complex::checkIfALLFacesAreFlat(){
    for (auto& face : faces){
        // if face is a triangle hence flat by definition cycle
        if(face->vertices.size() == 3) continue;

        // span{v1, v2} define the plane the face is assumed to be on. v2 is orthogonal to v1.
        Vector3D const& v1 = face->vertices[1]->vertex - face->vertices[0]->vertex;
        Vector3D const& v_temp = face->vertices[2]->vertex - face->vertices[1]->vertex;
        Vector3D const& v2 = v_temp - (ScalarProd(v1, v_temp)/ScalarProd(v1, v1))*v1;

        std::cout << "ScalarProduct(v1, v2) = " << ScalarProd(v2, v1) << std::endl;
        
        for(std::size_t i=2; i<face->vertices.size(); ++i){
            Vector3D const& v = face->vertices[(i+1)%face->vertices.size()]->vertex - face->vertices[i]->vertex;
            Vector3D const& v_perp = v - ((ScalarProd(v,v1) / ScalarProd(v1,v1)) * v1 + (ScalarProd(v,v2) / ScalarProd(v2,v2)) * v2);

            std::cout << "Edge: " << i+1 << " size of perp : " << ScalarProd(v_perp, v_perp) << std::endl;

            if(ScalarProd(v_perp, v_perp) > 1e-14){
                std::cout << "ERROR: face is not flat! Make sure all the vertices of face " << face->index << " are on the same plane!!" << std::endl;
                std::cout << "Edge: " << i+1 << " size of perp : " << ScalarProd(v_perp, v_perp) << std::endl;
                std::cout << "v = " << v.x << ", " << v.y << ", " << v.z << ", " << std::endl;
                std::cout << "v1 = " << v1.x << ", " << v1.y << ", " << v1.z << ", " << std::endl;
                std::cout << "v2 = " << v2.x << ", " << v2.y << ", " << v2.z << ", " << std::endl;
                std::cout << ScalarProd(v,v1) << ", " << ScalarProd(v,v2) << ", " << ScalarProd(v1,v1) << ", " << ScalarProd(v2,v2) << std::endl;
                return false;
            }
        }
    }
    return true;
}

void PL_Complex::detectFeatures(double const sharpTheta, double const flatTheta){
    /* Detect Sharp Edges */
    for(auto& edge : edges){
        if(edge->faces.size() == 1){
            sharp_edges.push_back(edge);
            edge->isSharp = true;
            continue;
        }

        double const dihedralAngle = edge->calcDihedralAngle();

        std::cout << "Edge " << edge->index << ", dihedral angle = " << dihedralAngle / M_PI << "*pi" << std::endl;

        if(dihedralAngle < (M_PI - sharpTheta)){
            sharp_edges.push_back(edge);
            edge->isSharp = true;
            continue;
        }

        if(dihedralAngle < (M_PI - flatTheta)){
            std::cout << "ERROR: dihedral angle of edge" << edge->index << " is not sharp nor flat!!!!!" << std::endl;
            exit(1);
        }

        edge->isSharp = false;
    }
    std::cout << "Sharp Edges:\n"; 

    for(auto& edge : sharp_edges){
        std::cout << "Edge: " << edge->index <<"\n"; 
    }
    std::cout << std::endl;
    
    for(auto& vertex : vertices){
        std::vector<std::shared_ptr<VoroCrustEdge>> vertex_sharp_edges;

        for(auto& edge : vertex->edges)
            if(edge->isSharp)
                vertex_sharp_edges.push_back(edge);

        if(vertex_sharp_edges.size() > 2){
            sharp_corners.push_back(vertex);
            vertex->isSharp = true;
            continue;
        }

        if(vertex_sharp_edges.size() < 2){
            std::cout << "ERROR: vertex has less than 2 sharp edges bordering it, vertex:"<< vertex->index << std::endl;
            exit(1);
        }

        auto& edge1 = vertex_sharp_edges[0];
        auto& edge2 = vertex_sharp_edges[1];

        Vector3D v1, v2;
        
        v1 = edge1->vertex2->vertex - edge1->vertex1->vertex;
        v2 = edge2->vertex2->vertex - edge2->vertex1->vertex;
        
        if(vertex->index == edge1->vertex2->index){
            v1 = -1*v1;
        }

        if(vertex->index == edge2->vertex2->index){
            v2 = -1*v2;
        }

        double const angle = std::acos(ScalarProd(v1, v2)/(abs(v1)*abs(v2)));

        std::cout << "vertex " << vertex->index << ", angle between edge " << edge1->index << " and edge " << edge2->index << " = " << angle / M_PI << "*pi" << std::endl;

        if(angle < (M_PI-sharpTheta)){
            sharp_corners.push_back(vertex);
            vertex->isSharp = true;
            continue;
        }

        vertex->isSharp = false;
    }

    std::cout << "\n\nSummary Sharp Features";
    std::cout << "\nSharp Edges:";
    for(auto& edge : sharp_edges) std::cout << "\n edge " << edge->index;
    std::cout << "\n\nSharp Corners:";
    for(auto& corner : sharp_corners) std::cout << "\n vertex " << corner->index;
    std::cout << std::endl;
}

std::string PL_Complex::repr() const{
    std::ostringstream s;


    for(auto& face : faces){
        s << "Face " << face->index << ": \n" << face->repr() << "\n";
    }

    return s.str();
}