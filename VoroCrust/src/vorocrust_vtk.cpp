#include "vorocrust_vtk.hpp"

#include "../../source/misc/utils.hpp"
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkPolyhedron.h>
#include <vtkDataArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkNew.h>

namespace vorocrust_vtk {
    void write_vtu_PL_Complex(std::filesystem::path const& filename, PL_Complex const& plc){
        auto const& vertices = plc.vertices;
        auto const& faces = plc.faces;

        std::size_t const num_vertices = vertices.size();
        std::size_t const num_faces = faces.size();

        vtkNew<vtkUnstructuredGrid> ugrid;

        vtkNew<vtkPoints> points;
        points->SetNumberOfPoints(num_vertices);

        for (std::size_t p = 0; p < num_vertices; ++p){
            points->SetPoint(vertices[p]->index, vertices[p]->vertex.x, vertices[p]->vertex.y, vertices[p]->vertex.z);
        }
        ugrid->SetPoints(points);

        
        ugrid->Allocate(num_faces);

        std::vector<vtkIdType> point_array_in_cell;
        for(std::size_t face_index=0; face_index < num_faces; ++face_index){
            vtkNew<vtkIdList> faces_vtk;
            point_array_in_cell.clear();

            auto& face = faces[face_index];
            size_t const Nvertices_in_face = face->vertices.size();

            faces_vtk->InsertNextId(Nvertices_in_face);
            for (std::size_t i = 0; i < Nvertices_in_face; ++i) {
                faces_vtk->InsertNextId(face->vertices[i]->index);
                point_array_in_cell.push_back(face->vertices[i]->index);
            }

            std::sort(point_array_in_cell.begin(), point_array_in_cell.end());
            point_array_in_cell = unique(point_array_in_cell);
            vtkIdType* ptIds = &point_array_in_cell[0];
            ugrid->InsertNextCell(VTK_POLYHEDRON, point_array_in_cell.size(), ptIds, 1, faces_vtk->GetPointer(0));        
        }

        std::vector<SurfacePatch> patches = plc.patches;

        vtkNew<vtkIntArray> data;
        data->SetName("Patch Index");
        data->SetNumberOfComponents(1);
        data->SetNumberOfValues(num_faces);
        
        for(std::size_t face_index=0; face_index < num_faces; ++face_index){
            Face const& face = faces[face_index];
            data->SetValue(face_index, -2);
        }
        
        std::cout << "Yo Yo Patches size " << patches.size() << std::endl;
        
        for(unsigned int p_index=0; p_index < patches.size(); ++p_index){
            SurfacePatch const& patch = patches[p_index];
            for(Face const& face : patch){
                data->SetValue(face->index, p_index);
            }
        }

        ugrid->GetCellData()->AddArray(data);


        vtkNew<vtkXMLUnstructuredGridWriter> writer;
        writer->SetCompressionLevel(1);
        writer->SetFileName(filename.c_str());
        writer->SetInputData(ugrid);
        writer->Write();
    }
}