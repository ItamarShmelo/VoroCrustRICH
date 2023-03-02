#include "hdf_write.hpp"
#include <filesystem>
#include "../../misc/hdf5_utils.hpp"
#include "../../misc/int2str.hpp"
#include "write_vtu_3d.hpp"

using namespace H5;

namespace
{

  template<class T> vector<T> read_vector_from_hdf5
  (const Group& file,
   const string& caption,
   const DataType& datatype)
  {
    DataSet dataset = file.openDataSet(caption);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out, nullptr);
    const size_t NX = static_cast<size_t>(dims_out[0]);
    vector<T> result(NX);
    dataset.read(&result[0], datatype);
    return result;
  }

  vector<double> read_double_vector_from_hdf5
  (const Group& file, string const& caption)
  {
    return read_vector_from_hdf5<double>
      (file,
       caption,
       PredType::NATIVE_DOUBLE);
  }

  vector<int> read_int_vector_from_hdf5
  (const Group& file,
   const string& caption)
  {
    return read_vector_from_hdf5<int>
      (file,
       caption,
       PredType::NATIVE_INT);
  }

  vector<size_t> read_sizet_vector_from_hdf5
  (const Group& file,
   const string& caption)
  {
    return read_vector_from_hdf5<size_t>
      (file,
       caption,
       PredType::NATIVE_ULLONG);
  }
}

DiagnosticAppendix3D::~DiagnosticAppendix3D() {}

Snapshot3D::Snapshot3D(void) :
  mesh_points(),
#ifdef RICH_MPI
  proc_points(),
#endif
  volumes(),
  cells(),
  time(),
  cycle(),
  tracerstickernames(),ll(Vector3D()),ur(Vector3D()){}

Snapshot3D::Snapshot3D(const Snapshot3D& source) :
  mesh_points(source.mesh_points),
#ifdef RICH_MPI
  proc_points(source.proc_points),
#endif
  volumes(source.volumes),
  cells(source.cells),
  time(source.time),
  cycle(source.cycle),
  tracerstickernames(source.tracerstickernames),ll(source.ll),ur(source.ur)
{}

void WriteVoronoi(Voronoi3D const& tri, std::string const& filename)
{
  H5File file(H5std_string(filename), H5F_ACC_TRUNC);
  vector<double> x, y, z, vx, vy, vz;
  vector<size_t> Nfaces;
  vector<size_t> Nvert;
  vector<size_t> FacesInCell;
  vector<size_t> VerticesInFace;
  size_t Npoints = tri.GetPointNo();
  for (size_t i = 0; i < Npoints; ++i)
    {
      x.push_back(tri.GetMeshPoint(i).x);
      y.push_back(tri.GetMeshPoint(i).y);
      z.push_back(tri.GetMeshPoint(i).z);
    }
  write_std_vector_to_hdf5(file, x, "mesh_point_x");
  write_std_vector_to_hdf5(file, y, "mesh_point_y");
  write_std_vector_to_hdf5(file, z, "mesh_point_z");
  x.clear();
  y.clear();
  z.clear();
  for (size_t i = 0; i < tri.GetTotalPointNumber(); ++i)
    {
      x.push_back(tri.GetMeshPoint(i).x);
      y.push_back(tri.GetMeshPoint(i).y);
      z.push_back(tri.GetMeshPoint(i).z);
    }
  write_std_vector_to_hdf5(file, x, "all_mesh_point_x");
  write_std_vector_to_hdf5(file, y, "all_mesh_point_y");
  write_std_vector_to_hdf5(file, z, "all_mesh_point_z");
  for (size_t i = 0; i < Npoints; ++i)
    {
      Nfaces.push_back(tri.GetCellFaces(i).size());
      for (size_t j = 0; j < Nfaces.back(); ++j)
	FacesInCell.push_back(tri.GetCellFaces(i)[j]);
    }
  IntType datatype(PredType::NATIVE_ULLONG);
  datatype.setOrder(H5T_ORDER_LE);
  write_std_vector_to_hdf5(file, Nfaces, "Number_of_faces_in_cell", datatype);
  write_std_vector_to_hdf5(file, FacesInCell, "Faces_in_cell", datatype);
  Npoints = tri.GetFacePoints().size();
  for (size_t i = 0; i < Npoints; ++i)
    {
      vx.push_back(tri.GetFacePoints()[i].x);
      vy.push_back(tri.GetFacePoints()[i].y);
      vz.push_back(tri.GetFacePoints()[i].z);
    }
  write_std_vector_to_hdf5(file, vx, "vertice_x");
  write_std_vector_to_hdf5(file, vy, "vertice_y");
  write_std_vector_to_hdf5(file, vz, "vertice_z");
  Npoints = tri.GetTotalFacesNumber();
  for (size_t i = 0; i < Npoints; ++i)
    {
      Nvert.push_back(tri.GetPointsInFace(i).size());
      for (size_t j = 0; j < Nvert.back(); ++j)
	VerticesInFace.push_back(tri.GetPointsInFace(i)[j]);
    }
  write_std_vector_to_hdf5(file, Nvert, "Number_of_vertices_in_face", datatype);
  write_std_vector_to_hdf5(file, VerticesInFace, "Vertices_in_face", datatype);
}

void WriteSnapshot3D(HDSim3D const& sim, std::string const& filename,
		     const vector<DiagnosticAppendix3D*>& appendices
#ifdef RICH_MPI
		     ,bool mpi_write
#endif // RICH_MPI
, bool const write_vtu
		     )
{
#ifdef RICH_MPI
  int rank = 0;
  int ws = 0;
  if (mpi_write)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &ws);
    }
#endif
  H5File file;
#ifdef RICH_MPI
  if (rank == 0)
    {
#endif // RICH_MPI
      H5File file2(H5std_string(filename), H5F_ACC_TRUNC);
      file2.close();
      file.openFile(H5std_string(filename), H5F_ACC_RDWR);
#ifdef RICH_MPI
    }
#endif // RICH_MPI
  std::vector<std::vector<double> > vtu_cell_variables;
  std::vector<std::string> vtu_cell_variable_names;
  
  std::vector<std::string> vtu_cell_vectors_names;
	std::vector<std::vector<Vector3D> > vtu_cell_vectors;
  Tessellation3D const& tess = sim.getTesselation();
  Group writegroup;
#ifdef RICH_MPI
  if (mpi_write)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      int dummy = 0;
      if (rank > 0)
	{
	  MPI_Recv(&dummy, 1, MPI_INT, rank - 1, 343, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  file.openFile(H5std_string(filename), H5F_ACC_RDWR);
	}
      file.createGroup("/rank" + int2str(rank));
      writegroup = file.openGroup("/rank" + int2str(rank));
    }
  else
    writegroup = file.openGroup("/");
#else
  writegroup = file.openGroup("/");
#endif
	
  vector<ComputationalCell3D> const& cells = sim.getCells();

  size_t Ncells = tess.GetPointNo();

  std::vector<double> box(6);
  box[0] = tess.GetBoxCoordinates().first.x;
  box[1] = tess.GetBoxCoordinates().first.y;
  box[2] = tess.GetBoxCoordinates().first.z;
  box[3] = tess.GetBoxCoordinates().second.x;
  box[4] = tess.GetBoxCoordinates().second.y;
  box[5] = tess.GetBoxCoordinates().second.z;
#ifdef RICH_MPI
  if(rank==0)
#endif // RICH_MPI
    write_std_vector_to_hdf5(file, box, "Box");

  vector<double> temp(Ncells);
  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetMeshPoint(i).x;
  write_std_vector_to_hdf5(writegroup, temp, "X");

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetMeshPoint(i).y;
  write_std_vector_to_hdf5(writegroup, temp, "Y");

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetMeshPoint(i).z;
  write_std_vector_to_hdf5(writegroup, temp, "Z");
  if(write_vtu)
  {
    vtu_cell_vectors_names.push_back("Coordinates");
    std::vector<Vector3D> vel(Ncells);
    for (size_t i = 0; i < Ncells; ++i)
      vel[i] = tess.GetMeshPoint(i);
    vtu_cell_vectors.push_back(vel);
  }
  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetCellCM(i).x;
  write_std_vector_to_hdf5(writegroup, temp, "CMx");

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetCellCM(i).y;
  write_std_vector_to_hdf5(writegroup, temp, "CMy");

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetCellCM(i).z;
  write_std_vector_to_hdf5(writegroup, temp, "CMz");

#ifdef RICH_MPI
  // MPI
  if (rank == 0)
    {
      Group mpi = writegroup.createGroup("/mpi");
      Tessellation3D const& tproc = sim.getProcTesselation();
      Ncells = tproc.GetPointNo();
      temp.resize(Ncells);
      for (size_t i = 0; i < Ncells; ++i)
	temp[i] = tproc.GetMeshPoint(i).x;
      write_std_vector_to_hdf5(mpi, temp, "proc_X");

      for (size_t i = 0; i < Ncells; ++i)
	temp[i] = tproc.GetMeshPoint(i).y;
      write_std_vector_to_hdf5(mpi, temp, "proc_Y");

      for (size_t i = 0; i < Ncells; ++i)
	temp[i] = tproc.GetMeshPoint(i).z;
      write_std_vector_to_hdf5(mpi, temp, "proc_Z");
      mpi.close();
    }
#endif

  Ncells = tess.GetPointNo();
  temp.resize(Ncells);
  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].density;
  write_std_vector_to_hdf5(writegroup, temp, "Density");
  if(write_vtu)
  {
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("Density");
  }

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].pressure;
  write_std_vector_to_hdf5(writegroup, temp, "Pressure");
  if(write_vtu)
  {
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("Pressure");
  }

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].internal_energy;
  write_std_vector_to_hdf5(writegroup, temp, "InternalEnergy");
  if(write_vtu)
  {
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("InternalEnergy");
  }

  std::vector<size_t> ids(Ncells);
  for (size_t i = 0; i < Ncells; ++i)
    ids[i] = cells[i].ID;
  write_std_vector_to_hdf5(writegroup, ids, "ID");
  if(write_vtu)
  {
    for (size_t i = 0; i < Ncells; ++i)
      temp[i] = ids[i];
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("ID");
  }
  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].velocity.x;
  write_std_vector_to_hdf5(writegroup, temp, "Vx");

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].velocity.y;
  write_std_vector_to_hdf5(writegroup, temp, "Vy");

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].velocity.z;
  write_std_vector_to_hdf5(writegroup, temp, "Vz");

  if(write_vtu)
  {
    vtu_cell_vectors_names.push_back("Velocity");
    std::vector<Vector3D> vel(Ncells);
    for (size_t i = 0; i < Ncells; ++i)
      vel[i] = cells[i].velocity;
    vtu_cell_vectors.push_back(vel);
  }

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].temperature;
  write_std_vector_to_hdf5(writegroup, temp, "Temperature");
  if(write_vtu)
  {
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("Temperature");
  }

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = cells[i].Erad;
  write_std_vector_to_hdf5(writegroup, temp, "Erad");
  if(write_vtu)
  {
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("Erad");
  }

  Group tracers, stickers;
#ifdef RICH_MPI
  if (mpi_write)
    {
      tracers = writegroup.createGroup("/rank" + int2str(rank) + "/tracers");
      stickers = writegroup.createGroup("/rank" + int2str(rank) + "/stickers");
    }
  else
    {
      tracers = writegroup.createGroup("/tracers");
      stickers = writegroup.createGroup("/stickers");
    }
#else
  tracers = writegroup.createGroup("/tracers");
  stickers = writegroup.createGroup("/stickers");
#endif
  for (size_t j = 0; j < ComputationalCell3D::tracerNames.size(); ++j)
    {
      for (size_t i = 0; i < Ncells; ++i)
	      temp[i] = cells[i].tracers[j];
      write_std_vector_to_hdf5(tracers, temp, ComputationalCell3D::tracerNames[j]);
      if(write_vtu)
      {
        vtu_cell_variables.push_back(temp);
        vtu_cell_variable_names.push_back(ComputationalCell3D::tracerNames[j]);
      }
    }

  for (size_t j = 0; j <ComputationalCell3D::stickerNames.size(); ++j)
    {
      for (size_t i = 0; i < Ncells; ++i)
	      temp[i] = cells[i].stickers[j];
      write_std_vector_to_hdf5(stickers, temp, ComputationalCell3D::stickerNames[j]);
      if(write_vtu)
      {
        vtu_cell_variables.push_back(temp);
        vtu_cell_variable_names.push_back(ComputationalCell3D::stickerNames[j]);
      }
    }

  for (size_t i = 0; i < Ncells; ++i)
    temp[i] = tess.GetVolume(i);
  write_std_vector_to_hdf5(writegroup, temp, "Volume");
  if(write_vtu)
  {
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("Volume");
  }
#ifdef RICH_MPI
  if (rank == 0)
    {
#endif // RICH_MPI
      vector<double> time(1, sim.getTime());
      write_std_vector_to_hdf5(file, time, "Time");

      vector<int> cycle(1, static_cast<int>(sim.getCycle()));
      write_std_vector_to_hdf5(file, cycle, "Cycle");
#ifdef RICH_MPI
    }
#endif // RICH_MPI
  // Appendices
  for (size_t i = 0; i < appendices.size(); ++i)
    write_std_vector_to_hdf5(writegroup,(*(appendices.at(i)))(sim),appendices.at(i)->getName());
#ifdef RICH_MPI
  if (mpi_write)
    {
      if (rank < (ws - 1))
	{
	  int dummy = 0;
	  tracers.close();
	  stickers.close();
	  writegroup.close();
	  file.close();
	  MPI_Send(&dummy, 1, MPI_INT, rank + 1, 343, MPI_COMM_WORLD);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  else
    file.close();
#else
  file.close();
#endif
  if(write_vtu)
  {
    std::filesystem::path vtu_name(filename);
    vtu_name.replace_extension("vtu");
    write_vtu3d::write_vtu_3d(vtu_name, vtu_cell_variable_names, vtu_cell_variables, vtu_cell_vectors_names, vtu_cell_vectors, sim.getTime(), sim.getCycle(), tess);
  }
}

std::vector<Vector3D> ReadVoronoiPoints(std::string const& filename)
{
  std::vector<Vector3D> res;
  H5File file(filename, H5F_ACC_RDONLY);
  Group read_location = file.openGroup("/");
  const vector<double> x = read_double_vector_from_hdf5(read_location, "mesh_point_x");
  const vector<double> y = read_double_vector_from_hdf5(read_location, "mesh_point_y");
  const vector<double> z = read_double_vector_from_hdf5(read_location, "mesh_point_z");
  size_t const N = x.size();
  for(size_t i = 0; i < N; ++i)
    res.push_back(Vector3D(x[i], y[i], z[i]));
  return res;
}

Snapshot3D ReadSnapshot3D(const string& fname
#ifdef RICH_MPI
			  , bool mpi_write,int fake_rank
#endif
			  )
{
  Snapshot3D res;
  H5File file(fname, H5F_ACC_RDONLY);
  Group read_location = file.openGroup("/");
#ifdef RICH_MPI
  if (mpi_write)
    {
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (fake_rank >= 0)
	rank = fake_rank;
      read_location = file.openGroup("/rank" + int2str(rank));
    }
#endif
  // Mesh points
  {
    const vector<double> x = read_double_vector_from_hdf5(read_location, "X");
    const vector<double> y = read_double_vector_from_hdf5(read_location, "Y");
    const vector<double> z = read_double_vector_from_hdf5(read_location, "Z");
    res.mesh_points.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i)
      res.mesh_points.at(i) = Vector3D(x.at(i), y.at(i), z.at(i));
  }

#ifdef RICH_MPI
  // MPI
  Group mpi = file.openGroup("mpi");
  const vector<double> px = read_double_vector_from_hdf5(mpi, "proc_X");
  const vector<double> py = read_double_vector_from_hdf5(mpi, "proc_Y");
  const vector<double> pz = read_double_vector_from_hdf5(mpi, "proc_Z");
  res.proc_points.resize(px.size());
  for (size_t i = 0; i < px.size(); ++i)
    res.proc_points.at(i) = Vector3D(px.at(i), py.at(i), pz.at(i));
#endif
  const vector<double> box = read_double_vector_from_hdf5(file, "Box");
  res.ll.Set(box[0], box[1], box[2]);
  res.ur.Set(box[3], box[4], box[5]);
  // Hydrodynamic
  {
    const vector<double> density = read_double_vector_from_hdf5(read_location, "Density");
    const vector<double> Erad = read_double_vector_from_hdf5(read_location, "Erad");
    const vector<double> temperature = read_double_vector_from_hdf5(read_location, "Temperature");
    const vector<double> pressure = read_double_vector_from_hdf5(read_location, "Pressure");
    const vector<double> energy = read_double_vector_from_hdf5(read_location, "InternalEnergy");
    vector<size_t> IDs(density.size(), 0);
    hsize_t objcount = read_location.getNumObjs();
    for (hsize_t i = 0; i < objcount; ++i)
      {
	std::string name = read_location.getObjnameByIdx(i);
	if (name.compare(std::string("ID"))==0)
	  IDs = read_sizet_vector_from_hdf5(read_location, "ID");
      }
    const vector<double> x_velocity = read_double_vector_from_hdf5(read_location, "Vx");
    const vector<double> y_velocity = read_double_vector_from_hdf5(read_location, "Vy");
    const vector<double> z_velocity = read_double_vector_from_hdf5(read_location, "Vz");

    Group g_tracers = read_location.openGroup("tracers");
    Group g_stickers = read_location.openGroup("stickers");
    vector<vector<double> > tracers(g_tracers.getNumObjs());
    vector<string> tracernames(tracers.size());
    for (hsize_t n = 0; n < g_tracers.getNumObjs(); ++n)
      {
	const H5std_string name = g_tracers.getObjnameByIdx(n);
	tracernames[n] = name;
	tracers[n] = read_double_vector_from_hdf5(g_tracers, name);
      }

    vector<vector<int> > stickers(g_stickers.getNumObjs());
    vector<string> stickernames(stickers.size());
    for (hsize_t n = 0; n < g_stickers.getNumObjs(); ++n)
      {
	const H5std_string name = g_stickers.getObjnameByIdx(n);
	stickernames[n] = name;
	stickers[n] = read_int_vector_from_hdf5(g_stickers, name);
      }
    res.tracerstickernames.first = tracernames;
    res.tracerstickernames.second = stickernames;
    res.cells.resize(density.size());
    for (size_t i = 0; i < res.cells.size(); ++i)
      {
	res.cells.at(i).density = density.at(i);
  res.cells.at(i).Erad = Erad.at(i);
  res.cells.at(i).temperature = temperature.at(i);
	res.cells.at(i).pressure = pressure.at(i);
	res.cells.at(i).internal_energy = energy.at(i);
	res.cells.at(i).ID = IDs[i];
	res.cells.at(i).velocity.x = x_velocity.at(i);
	res.cells.at(i).velocity.y = y_velocity.at(i);
	res.cells.at(i).velocity.z = z_velocity.at(i);
	for (size_t j = 0; j < tracernames.size(); ++j)
	  res.cells.at(i).tracers.at(j) = tracers.at(j).at(i);
	for (size_t j = 0; j < stickernames.size(); ++j)
	  res.cells.at(i).stickers.at(j) = stickers.at(j).at(i) == 1;
      }
  }

  // Misc
  {
    const vector<double> time = read_double_vector_from_hdf5(file, "Time");
    res.time = time.at(0);

    const vector<double> volume = read_double_vector_from_hdf5(read_location, "Volume");
    res.volumes = volume;
    const vector<int> cycle = read_int_vector_from_hdf5(file, "Cycle");
    res.cycle = cycle.at(0);
  }

  return res;
}

#ifdef RICH_MPI
Snapshot3D ReDistributeData3D(string const& filename, Tessellation3D const& proctess, size_t snapshot_number,bool mpi_write)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Read the data
  Snapshot3D snap/*,temp*/;
  for (int i = 0; i < static_cast<int>(snapshot_number); ++i)
    {
      /*
      if(mpi_write)
	temp = ReadSnapshot3D(filename, true,i);
      else
	temp = ReadSnapshot3D(filename + int2str(i) + ".h5", false);
      */
      Snapshot3D temp(mpi_write ?
		      ReadSnapshot3D(filename, true, i) :
		      ReadSnapshot3D(filename + int2str(i) + ".h5", false));
      if (i == 0)
	{
	  snap.time = temp.time;
	  snap.cycle = temp.cycle;
	  snap.tracerstickernames = temp.tracerstickernames;
	  snap.ll = temp.ll;
	  snap.ur = temp.ur;
	  snap.proc_points = temp.proc_points;
	}
      size_t N = temp.cells.size();
      for (size_t j = 0; j < N; ++j)
	{
	  if (PointInPoly(proctess, temp.mesh_points[j], static_cast<size_t>(rank)))
	    {
	      snap.cells.push_back(temp.cells[j]);
	      snap.mesh_points.push_back(temp.mesh_points[j]);
	    }
	}
    }
  return snap;
}
#endif
