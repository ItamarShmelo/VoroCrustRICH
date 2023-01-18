#include "PL_Complex.hpp"
#include <iostream>
#include <cmath>
#include <queue>
#include "../../../source/misc/utils.hpp"

PL_Complex::PL_Complex(std::vector<Vector3D> const &vertices_) : vertices(),
                                                                 edges(),
                                                                 faces(),
                                                                 sharp_edges(),
                                                                 sharp_corners(),
                                                                 creases(),
                                                                 patches()
{
    std::size_t index = 0;
    for (auto &vertex : vertices_)
    {
        vertices.push_back(std::make_shared<VoroCrustVertex>(vertex, index));
        index++;
    }
}

Edge PL_Complex::addEdge(Vertex const &v1, Vertex const &v2)
{
    // first check if edge was already created by a different face.
    for (auto &edge : edges)
    {
        if (edge->checkIfEqual(v1, v2))
        {
            return edge;
        }
    }

    // create new edge starting at v1 and ending at v1
    auto new_edge_ptr = std::make_shared<VoroCrustEdge>(v1, v2, edges.size());

    // add the edge to the vertices constructing it
    new_edge_ptr->vertex1->addEdge(new_edge_ptr);
    new_edge_ptr->vertex2->addEdge(new_edge_ptr);

    // update edges vector
    edges.push_back(new_edge_ptr);
    return new_edge_ptr;
}

void PL_Complex::addFace(std::vector<unsigned int> const &indices)
{

    if (indices.size() != unique(indices).size())
    {
        std::cout << "ERROR: Repeated indices in PL_Complex::addFace!" << std::endl;
        exit(1);
    }

    if (indices.size() < 3)
    {
        std::cout << "ERROR: Face has to have at least three vertices!" << std::endl;
        exit(1);
    }

    std::vector<Vertex> face_vertices = std::vector<Vertex>();

    for (auto &index : indices)
        face_vertices.push_back(vertices[index]);

    Face new_face_ptr = std::make_shared<VoroCrustFace>(face_vertices, faces.size());

    // add new face to the vertices constucting it
    for (Vertex vertex_ptr : face_vertices)
    {
        vertex_ptr->addFace(new_face_ptr);
    }

    // create new edges and add the to edges of new face
    for (std::size_t i = 0; i < new_face_ptr->vertices.size(); ++i)
    {
        auto new_edge_ptr = this->addEdge(new_face_ptr->vertices[i], new_face_ptr->vertices[(i + 1) % (new_face_ptr->vertices.size())]);
        new_edge_ptr->addFace(new_face_ptr);
        new_face_ptr->addEdge(new_edge_ptr);

        // if edge was already created update the neighbors array of the faces on the sides of the edge
        if (new_edge_ptr->faces.size() == 2)
        {
            new_face_ptr->neighbors.push_back(new_edge_ptr->faces[0]);
            new_edge_ptr->faces[0]->neighbors.push_back(new_edge_ptr->faces[1]);

            std::cout << "Neigbors :: faces " << new_face_ptr->index << " and face " << new_edge_ptr->faces[0]->index << std::endl;
        }
    }

    faces.push_back(new_face_ptr);
}

bool PL_Complex::checkAllVerticesAreOnFace()
{
    for (auto &vertex : vertices)
        if (vertex->faces.size() == 0)
        {
            std::cout << "ERROR: vertex " << vertex->index << " is not part of a face!!";
            return false;
        }

    return true;
}

bool PL_Complex::checkIfALLFacesAreFlat()
{
    for (auto &face : faces)
    {
        // if Face is a triangle hence flat by definition cycle
        if (face->vertices.size() == 3)
            continue;

        // span{v1, v2} define the plane the face is assumed to be on. v2 is orthogonal to v1.
        Vector3D const &v1 = face->vertices[1]->vertex - face->vertices[0]->vertex;
        Vector3D const &v_temp = face->vertices[2]->vertex - face->vertices[1]->vertex;
        Vector3D const &v2 = v_temp - (ScalarProd(v1, v_temp) / ScalarProd(v1, v1)) * v1; // make v2 orthogonal using Grahm-Shmidt

        std::cout << "ScalarProduct(v1, v2) = " << ScalarProd(v2, v1) << std::endl;

        // Calculate the magnitute of ||v_perp||=||v - projection(v)|| and check that it is less then some eps.
        // projection(v) is the projection on span{v1, v2}
        for (std::size_t i = 2; i < face->vertices.size(); ++i)
        {
            Vector3D const &v = face->vertices[(i + 1) % face->vertices.size()]->vertex - face->vertices[i]->vertex;
            Vector3D const &v_perp = v - ((ScalarProd(v, v1) / ScalarProd(v1, v1)) * v1 + (ScalarProd(v, v2) / ScalarProd(v2, v2)) * v2);

            std::cout << "Edge: " << i + 1 << " size of perp : " << abs(v_perp) << std::endl;

            if (abs(v_perp) > 1e-14)
            {
                std::cout << "ERROR: face is not flat! Make sure all the vertices of face " << face->index << " are on the same plane!!" << std::endl;
                std::cout << "Edge: " << i + 1 << " size of perp : " << abs(v_perp) << std::endl;
                return false;
            }
        }
    }
    return true;
}

void PL_Complex::detectFeatures(double const sharpTheta, double const flatTheta)
{
    /* Detect Sharp Edges */
    for (auto &edge : edges)
    {
        // if Edge is incident to only one face it is sharp.
        if (edge->faces.size() == 1)
        {
            sharp_edges.push_back(edge);
            edge->isSharp = true;
            continue;
        }

        // if the dihedral angle of an Edge is less than `PI-sharpTheta` it is sharp.
        double const dihedralAngle = edge->calcDihedralAngle();

        std::cout << "Edge " << edge->index << ", dihedral angle = " << dihedralAngle / M_PI << "*pi" << std::endl;

        if (dihedralAngle < (M_PI - sharpTheta))
        {
            sharp_edges.push_back(edge);
            edge->isSharp = true;
            continue;
        }

        if (dihedralAngle < (M_PI - flatTheta))
        {
            std::cout << "ERROR: dihedral angle of edge" << edge->index << " is not sharp nor flat!!!!!" << std::endl;
            exit(1);
        }

        edge->isSharp = false;
    }

    std::cout << "Sharp Edges:\n";

    for (auto &edge : sharp_edges)
    {
        std::cout << "Edge: " << edge->index << "\n";
    }
    std::cout << std::endl;

    /* Detect Sharp Corners */
    for (auto &vertex : vertices)
    {
        std::vector<Edge> vertex_sharp_edges;

        // find sharp Edges incident to Vertex
        for (auto &edge : vertex->edges)
            if (edge->isSharp)
                vertex_sharp_edges.push_back(edge);

        // if Vertex is shared by 0 sharp Edges it is not a sharp corner.
        if(vertex_sharp_edges.empty()){
            vertex->isSharp = false;
            continue;
        }

        // if a Vertex is shared by more than 2 sharp Edges then it is a sharp corner.
        if (vertex_sharp_edges.size() > 2)
        {
            sharp_corners.push_back(vertex);
            vertex->isSharp = true;
            continue;
        }

        // if Vertex is shared by 2 sharp edges if the angle between them is less than `PI-sharpTheta`
        // then the vertex is sharp.
        if (vertex_sharp_edges.size() == 2)
        {
            
            auto &edge1 = vertex_sharp_edges[0];
            auto &edge2 = vertex_sharp_edges[1];

            Vector3D v1, v2;

            // the vectors representing the edges
            v1 = edge1->vertex2->vertex - edge1->vertex1->vertex;
            v2 = edge2->vertex2->vertex - edge2->vertex1->vertex;

            // gurentee that both v1, v2 start at current Vertex.
            if (vertex->index == edge1->vertex2->index)
            {
                v1 = -1 * v1;
            }

            if (vertex->index == edge2->vertex2->index)
            {
                v2 = -1 * v2;
            }

            double const angle = std::acos(ScalarProd(v1, v2) / (abs(v1) * abs(v2)));

            std::cout << "vertex " << vertex->index << ", angle between edge " << edge1->index << " and edge " << edge2->index << " = " << angle / M_PI << "*pi" << std::endl;

            if (angle < (M_PI - sharpTheta))
            {
                sharp_corners.push_back(vertex);
                vertex->isSharp = true;
                continue;
            }

            vertex->isSharp = false;
            continue;
        } 

        // if Vertex is shared by only one sharp edge then it is sharp.
        vertex->isSharp = true;
    }

    std::cout << "\n\nSummary Sharp Features";
    std::cout << "\nSharp Edges:";
    for (auto &edge : sharp_edges)
        std::cout << "\n edge " << edge->index;
    std::cout << "\n\nSharp Corners:";
    for (auto &corner : sharp_corners)
        std::cout << "\n vertex " << corner->index;
    std::cout << std::endl;

    buildCreases();

    buildSurfacePatches();
}

void PL_Complex::buildCreases()
{
    
    for (auto &edge : sharp_edges)
    {
        if (edge->isCreased)
            continue;

        Crease const &new_crease = createCrease(edge);
        creases.push_back(new_crease);

        std::cout << "\nCrease number : " << creases.size() << "\n";
        std::cout << "Edges: ";

        for (Edge const &crease_edge : new_crease)
        {
            std::cout << crease_edge->index << " -> ";
        }

        std::cout << new_crease[0]->index << std::endl;
    }
}

Crease PL_Complex::createCrease(Edge const &edge)
{
    /* Creates the Creases using the flood fill algorithm across vertices which are not sharp corners */

    Crease crease;
    std::queue<Edge> queue;

    queue.push(edge);

    while (not queue.empty())
    {
        Edge const &curr_edge = queue.front();
        queue.pop();

        if (not curr_edge->isCreased)
        {
            crease.push_back(curr_edge);
            curr_edge->isCreased = true;

            if (not curr_edge->vertex1->isSharp)
            {

                for (Edge const &vertex_edge : curr_edge->vertex1->edges)
                {
                    if (vertex_edge->isSharp && not vertex_edge->isCreased)
                    {
                        queue.push(vertex_edge);
                        break;
                    }
                }
            }

            if (not curr_edge->vertex2->isSharp)
            {
                for (Edge const &vertex_edge : curr_edge->vertex2->edges)
                {
                    if (vertex_edge->isSharp && not vertex_edge->isCreased)
                    {
                        queue.push(vertex_edge);
                        break;
                    }
                }
            }
        }
    }

    return crease;
}

void PL_Complex::buildSurfacePatches()
{   
    for (Face &face : faces)
    {
        if (face->isPatched) 
            continue;
        
        SurfacePatch const& new_patch = createSurfacePatch(face);
        patches.push_back(new_patch);

        std::cout << "\n SurfacePatch number : " << patches.size() << "\n"; 
        std::cout << "Faces : ";

        for(Face const& patch_face : new_patch){
            std::cout << patch_face->index << ", ";
        } 
        std::cout << std::endl;        
    }
}

SurfacePatch PL_Complex::createSurfacePatch(Face const &face)
{
    /* Create Surface Patch using the flood fill algorithm across Faces which share a flat Edge */
    SurfacePatch patch;
    std::queue<Face> queue;

    queue.push(face);

    while(not queue.empty()){
        Face const& curr_face = queue.front();
        queue.pop();

        if(not curr_face->isPatched){
            patch.push_back(curr_face);
            curr_face->isPatched = true;

            for(Edge const& edge : curr_face->edges){
                if(edge->isSharp) continue;

                for(Face const& edge_face : edge->faces){
                    if(not edge_face->isPatched){
                        queue.push(edge_face);
                    }
                }
            }
        }
    }

    return patch;
}

std::string PL_Complex::repr() const
{
    std::ostringstream s;

    for (auto &face : faces)
    {
        s << "Face " << face->index << ": \n"
          << face->repr() << "\n";
    }

    return s.str();
}