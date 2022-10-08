#ifndef WRITE_VTK_HPP_
#define WRITE_VTK_HPP_

#include <vector>
#include <string>

namespace write_vtk{

/**
 * @brief writes an ASCII vtk snapshot filefile
 * 
 * @param file_name name of the file (without the ".vtk" suffix)
 * @param x_vertices 
 * @param y_vertices 
 * @param cells_num_vertices 
 * @param cell_variable_names 
 * @param cell_variables 
 * @param time 
 * @param cycle 
 */
void write_vtk(std::string const& file_name,
			   std::vector<double> const& x_vertices,
			   std::vector<double> const& y_vertices,
			   std::vector<std::size_t> const& cells_num_vertices, //number of vertices each cell has
			   std::vector<std::string> const& cell_variable_names,
			   std::vector<std::vector<double>> const& cell_variables,
			   double const time,
			   std::size_t cycle);

} //namespace

#endif /* WRITE_VTK_HPP_ */