#include "PL_Complex.hpp"
#include <iostream>
#include <cmath>
#include <numeric>
#include <functional>
#include <queue>
#include "../unsorted_unique.hpp"

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
    for (auto const& edge : edges)
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
    //! WARNING: unique needs vector to be sorted so this does nothing actually... 
    if (indices.size() != unsorted_unique<unsigned int>(indices).size())
    {
        std::cout << "ERROR: Repeated indices in PL_Complex::addFace!" << std::endl;
        exit(1);
    }

    if (indices.size() < 3)
    {
        std::cout << "ERROR: Face has to have at least three vertices!" << std::endl;
        exit(1);
    }

    std::vector<Vertex> face_vertices;

    for (unsigned int const index : indices)
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

            // std::cout << "Neigbors :: faces " << new_face_ptr->index << " and face " << new_edge_ptr->faces[0]->index << std::endl;
        }
    }

    faces.push_back(new_face_ptr);
}

bool PL_Complex::checkAllVerticesAreUnique() {
    std::vector<Vector3D> vec_vertices;
    vec_vertices.reserve(vertices.size());

    for(auto const& vertex : vertices)
        vec_vertices.push_back(vertex->vertex);
    
    auto const& new_vec_vertices = unsorted_unique<Vector3D>(vec_vertices);

    if(new_vec_vertices.size() != vec_vertices.size()){
        std::cout << "ERROR: some vertices are not unique" << std::endl;
        std::cout << "number : " << (vec_vertices.size()-new_vec_vertices.size()) << std::endl;
        return false;
    }

    return true;
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

        // if a Vertex is shared by more than 2 sharp Edges or by one sharp edge then it is a sharp corner.
        if (vertex_sharp_edges.size() > 2 || vertex_sharp_edges.size() == 1)
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
            crease_edge->crease_index = creases.size() - 1;
            std::cout << crease_edge->index << " -> ";
        }

        std::cout << new_crease[0]->index << std::endl;
    }
}

void PL_Complex::orderCrease(Crease &crease){
    Edge start;
    bool foundStart;
    // look for the start of the crease (if crease is a circle it will take the last element as the start)
    for(std::size_t i=0; i<crease.size(); ++i){
        start = crease[i];
        foundStart = true;
        
        for(std::size_t j=0; j<crease.size(); ++j){
            if(start->vertex1->index == crease[j]->vertex2->index){
                foundStart = false;
            }
        }

        if(foundStart) break;
    }

    Crease orderedCrease;

    orderedCrease.push_back(start);

    // build an ordered crease i.e v2[i-1] = v1[i]
    for(std::size_t i=0; i<crease.size()-1; ++i){
        bool found = false;
        for(std::size_t j=0; j<crease.size(); ++j){
            if(orderedCrease[i]->vertex2->index == crease[j]->vertex1->index){
                orderedCrease.push_back(crease[j]);
                found = true;
                break;
            }
        }
        if(not found){
            std::cout << "ERROR IN ORDERING THE CREASES" << std::endl;
            exit(1);
        }
    }

    crease = orderedCrease;
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
                        vertex_edge->orientWithRespectTo(curr_edge);
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
                        vertex_edge->orientWithRespectTo(curr_edge);
                        queue.push(vertex_edge);
                        break;
                    }
                }
            }
        }
    }

    orderCrease(crease);

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
            patch_face->patch_index = patches.size()-1;
            std::cout << patch_face->index << ", ";
        } 
        std::cout << std::endl;        
    }
}

void SurfacePatch::findCreasesAndCorners(){
    for(Face const& face : this->patch){
        for(Edge const& edge : face->edges){
            if(not edge->isSharp) continue;

            if(std::find(patch_creases.begin(), patch_creases.end(), edge->crease_index) == patch_creases.end()){
                patch_creases.push_back(edge->crease_index);
            }
        }

        for(Vertex const& vertex : face->vertices){
            if(not vertex->isSharp) continue;

            if(std::find(patch_corners.begin(), patch_corners.end(), vertex->index) == patch_corners.end()){
                patch_corners.push_back(vertex->index);
            }
        }
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
                        edge_face->orientWithRespectTo(curr_face);
                        queue.push(edge_face);
                    }
                }
            }
        }
    }
    
    patch.findCreasesAndCorners();
    return patch;
}

std::array<double, 6> PL_Complex::getBoundingBox() const {
    double ll_x, ll_y, ll_z; // lower left
    double ur_x, ur_y, ur_z; // upper right
    ll_x = ll_y = ll_z = std::numeric_limits<double>::max();
    ur_x = ur_y = ur_z = -std::numeric_limits<double>::max();

    for(Vertex const& vertex : vertices){
        Vector3D const& p = vertex->vertex;

        ll_x = std::min(p.x, ll_x);
        ll_y = std::min(p.y, ll_y);
        ll_z = std::min(p.z, ll_z);

        ur_x = std::max(p.x, ur_x);
        ur_y = std::max(p.y, ur_y);
        ur_z = std::max(p.z, ur_z);    
    }

    return std::array<double, 6>({ll_x, ll_y, ll_z, ur_x, ur_y, ur_z});
}

PL_Complex::Location PL_Complex::determineLocation(Vector3D const& p) const {
    // determine if a point is in or out using the ray casting algorithm
    int count = 0;
    for(Face const& face : faces) {
        
        if(face->isPointCompletelyOffFace(p)) continue;

        auto const& [success, p_inter] = face->pointXYaxisRayIntersectsAt(p);

        if(not success) continue;

        if(p_inter.z > p.z && face->pointIsInsideFace(p_inter)){
            count++;
        }
    }

    if(count % 2 == 0){
        return Location::OUT;
    } else {
        return Location::IN;
    }
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