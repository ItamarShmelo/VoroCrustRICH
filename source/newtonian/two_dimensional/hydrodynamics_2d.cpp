#include "hydrodynamics_2d.hpp"
#include "../../tessellation/calc_face_vertex_velocity.hpp"

using std::max;

int get_other_index(const Edge& edge, const int index)
{
  if(edge.neighbors.first==index && edge.neighbors.second!=index)
    return edge.neighbors.second;
  else if(edge.neighbors.second==index && edge.neighbors.first!=index)
    return edge.neighbors.first;
  else
    throw UniversalError("Something went wrong in Hydrodynamics2D::get_other_index");
}

namespace {
  Vector2D calc_representing_point(Tessellation const& tess,
				   int index,
				   bool cm_flag)
  {
    if(cm_flag)
      return tess.GetCellCM(index);
    else
      return tess.GetMeshPoint(index);
  }

  Primitive initialize_single_cell(Tessellation const& tess,
				   int index,
				   bool cm_flag,
				   SpatialDistribution const& density,
				   SpatialDistribution const& pressure,
				   SpatialDistribution const& xvelocity,
				   SpatialDistribution const& yvelocity,
				   EquationOfState const& eos)
  {
    const Vector2D r = calc_representing_point(tess,index,cm_flag);
    return CalcPrimitive(density(r),
			 pressure(r),
			 Vector2D(xvelocity(r),
				  yvelocity(r)),
			 eos);
  }

  class CellInitializer: public LazyList<Primitive>
  {
  public:

    CellInitializer(Tessellation const& tess,
		    bool cm_flag,
		    SpatialDistribution const& density,
		    SpatialDistribution const& pressure,
		    SpatialDistribution const& xvelocity,
		    SpatialDistribution const& yvelocity,
		    EquationOfState const& eos):
      tess_(tess), cm_flag_(cm_flag),
      density_(density),
      pressure_(pressure),
      xvelocity_(xvelocity),
      yvelocity_(yvelocity),
      eos_(eos) {}

    Primitive operator[](size_t n) const
    {
      return initialize_single_cell(tess_,
				    static_cast<int>(n),
				    cm_flag_,
				    density_,
				    pressure_,
				    xvelocity_,
				    yvelocity_,
				    eos_);
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

    ~CellInitializer(void) {}

  private:
    Tessellation const& tess_;
    const bool cm_flag_;
    SpatialDistribution const& density_;
    SpatialDistribution const& pressure_;
    SpatialDistribution const& xvelocity_;
    SpatialDistribution const& yvelocity_;
    EquationOfState const& eos_;
  };
}

vector<Primitive> InitialiseCells
(SpatialDistribution const& density,
 SpatialDistribution const& pressure,
 SpatialDistribution const& xvelocity,
 SpatialDistribution const& yvelocity,
 EquationOfState const& eos,
 Tessellation const& tess,
 bool cm_value)
{
  return serial_generate
    (CellInitializer
     (tess, cm_value, density,
      pressure, xvelocity, yvelocity, eos));
}

namespace {
  class IntensiveInitializer: public LazyList<Conserved>
  {
  public:

    explicit IntensiveInitializer(vector<Primitive> const& cells):
      cells_(cells) {}

    Conserved operator[](size_t n) const
    {
      return Primitive2Conserved(cells_[n]);
    }

    size_t size(void) const
    {
      return cells_.size();
    }

  private:
    vector<Primitive> const& cells_;
  };
}

vector<Conserved> CalcConservedIntensive
(vector<Primitive> const& cells)
{
  return serial_generate(IntensiveInitializer(cells));
}

namespace {

  class CellEdgesGetter: public LazyList<Edge>
  {
  public:

    CellEdgesGetter(const Tessellation& tess, int n):
      tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

    size_t size(void) const
    {
      return edge_indices_.size();
    }

    Edge operator[](size_t i) const
    {
      return tess_.GetEdge(edge_indices_[i]);
    }

  private:
    const Tessellation& tess_;
    const vector<int> edge_indices_;
  };

  class ExtensiveInitializer: public LazyList<Conserved>
  {
  public:

    ExtensiveInitializer(const vector<Conserved>& intensive,
			 const Tessellation& tess,
			 const PhysicalGeometry& pg):
      intensive_(intensive), tess_(tess), pg_(pg)  {}

    Conserved operator[](size_t n) const
    {
      return pg_.calcVolume(serial_generate(CellEdgesGetter(tess_,static_cast<int>(n))))*
	intensive_[n];
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

  private:
    const vector<Conserved>& intensive_;
    const Tessellation& tess_;
    const PhysicalGeometry& pg_;
  };
}

vector<Conserved> CalcConservedExtensive
(const vector<Conserved>& cons_int,
 const Tessellation& tess,
 const PhysicalGeometry& pg)
{
  return serial_generate(ExtensiveInitializer(cons_int, tess, pg));
}

double TimeStepForCell(Primitive const& cell,double width,
		       vector<Vector2D> const& face_velocites)
{
  double max_fv=0;
  for(size_t i=0;i<face_velocites.size();++i)
    max_fv = max(max_fv,abs(face_velocites[i]-cell.Velocity));
  return width/(cell.SoundSpeed+max_fv);
}

namespace {
  class NewPointPosition: public LazyList<Vector2D>
  {
  public:

    NewPointPosition(Tessellation const& tess,
		     vector<Vector2D> const& point_velocity,
		     double dt):
      tess_(tess),
      point_velocity_(point_velocity),
      dt_(dt) {}

    Vector2D operator[](size_t n) const
    {
      return tess_.GetMeshPoint(static_cast<int>(n))+dt_*point_velocity_[n];
    }

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

  private:
    Tessellation const& tess_;
    vector<Vector2D> const& point_velocity_;
    const double dt_;
  };
}

void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation& tessellation,vector<Vector2D> oldpoints)
{
  if(oldpoints.empty())
    {
      tessellation.Update(serial_generate
			  (NewPointPosition(tessellation,
					    pointvelocity,
					    dt)));
    }
  else
    {
      for(size_t i=0;i<oldpoints.size();++i)
	oldpoints[i]+=pointvelocity[i]*dt;
      tessellation.Update(oldpoints);
    }
}

#ifdef RICH_MPI
void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
		    double dt, Tessellation& tessellation,

		    Tessellation const& vproc,
		    vector<Vector2D> oldpoints)
{
  if(oldpoints.empty())
    tessellation.Update(serial_generate
			(NewPointPosition(tessellation,
					  pointvelocity,
					  dt)),vproc);
  else
    {
      for(size_t i=0;i<oldpoints.size();++i)
	oldpoints[i]+=pointvelocity[i]*dt;
      tessellation.Update(oldpoints,vproc);
    }
}
#endif // RICH_MPI

namespace {
  class IntensiveCalculator: public LazyList<Conserved>
  {
  public:

    IntensiveCalculator(const Tessellation& tess,
			const vector<Conserved>& extensive,
			const PhysicalGeometry& pg):
      tess_(tess), extensive_(extensive), pg_(pg) {}

    size_t size(void) const
    {
      return extensive_.size();
    }

    Conserved operator[](size_t i) const
    {
      return extensive_[i]/
	pg_.calcVolume(serial_generate(CellEdgesGetter(tess_,static_cast<int>(i))));
    }

  private:
    const Tessellation& tess_;
    const vector<Conserved>& extensive_;
    const PhysicalGeometry& pg_;
  };
}

vector<Conserved> calc_conserved_intensive
(const Tessellation& tess,
 const vector<Conserved>& extensive,
 const PhysicalGeometry& pg)
{
  return serial_generate(IntensiveCalculator(tess,extensive,pg));
}

void UpdateConservedIntensive(Tessellation const& tessellation,
			      vector<Conserved> const& conservedextensive,
			      vector<Conserved>& conservedintensive)
{
  conservedintensive.resize(conservedextensive.size());
  for(int i = 0; i<tessellation.GetPointNo(); i++){
    conservedintensive[static_cast<size_t>(i)] = conservedextensive[static_cast<size_t>(i)]/
      tessellation.GetVolume(i);
  }
}

namespace {

  std::pair<Conserved,bool> calc_safe_conserved(Conserved const& raw,
						bool density_floor,
						double min_density,
						double min_pressure,
						Primitive const& old,
						EquationOfState const& eos)
  {
    std::pair<Conserved, bool> res;
    res.first = raw;
    if(density_floor)
      {
	if (res.first.Mass < min_density)
	  {
	    res.first.Mass = min_density;
	    res.first.Momentum = old.Velocity*min_density;
	    const double kinetic_energy = 0.5*pow(abs(res.first.Momentum / res.first.Mass), 2);
	    res.first.Energy = res.first.Mass*kinetic_energy
	      + res.first.Mass*eos.dp2e(res.first.Mass, min_pressure);
	    res.second = true;
	  }
	const double kinetic_energy = 0.5*pow(abs(res.first.Momentum / res.first.Mass), 2);
	const double thermal_energy = res.first.Energy / res.first.Mass - kinetic_energy;
	const double pressure = eos.de2p(res.first.Mass, thermal_energy);
	if (pressure < min_pressure || res.second)
	  {
	    res.first.Energy = res.first.Mass*kinetic_energy 
	      + res.first.Mass*eos.dp2e(res.first.Mass, min_pressure);
	    res.second = true;
	  }
      }
    return res;
  }

#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
  __attribute__((noreturn))
#endif
  void update_primitives_rethrow(int cell_index,
				 UniversalError& eo)
  {
    eo.AddEntry("UpdatePrimitive data starts here",0);
    eo.AddEntry("cell index",static_cast<double>(cell_index));
    throw eo;
  }

  std::pair<Primitive,bool> regular_cell_evolve(Conserved const& intensive,
						bool density_floor,
						double min_density,
						double min_pressure,
						Primitive const& old,
						EquationOfState const& eos)
  {
    const std::pair<Conserved, bool> temp = calc_safe_conserved
      (intensive,density_floor, min_density,
       min_pressure, old, eos);
    return std::pair<Primitive,bool> (Conserved2Primitive(temp.first, eos),temp.second);
  }
}

vector<bool> UpdatePrimitives
(vector<Conserved> const& conservedintensive,
 EquationOfState const& eos,vector<Primitive>& cells,
 vector<CustomEvolution*> const& CellsEvolve,vector<Primitive> &old_cells,
 bool densityfloor,double densitymin,double pressuremin,Tessellation const&
 tess,double time,vector<vector<double> > const& tracers)
{
  cells.resize(conservedintensive.size());
  vector<bool> bres(cells.size(), false);
  for(int i=0;i < tess.GetPointNo(); i++)
    {
      try
	{
	  if (CellsEvolve[static_cast<size_t>(i)] == 0 || CellsEvolve[static_cast<size_t>(i)]->DensityFloorRelevant())
	    {
	      Primitive old_cell = densityfloor ? old_cells[static_cast<size_t>(i)] : cells[static_cast<size_t>(i)];
	      std::pair<Primitive, bool > res = regular_cell_evolve
		(conservedintensive[static_cast<size_t>(i)], densityfloor,
		 densitymin, pressuremin, old_cell, eos);
	      if (CellsEvolve[static_cast<size_t>(i)] != 0 && !res.second)
		cells[static_cast<size_t>(i)] = CellsEvolve[static_cast<size_t>(i)]->UpdatePrimitive
		  (conservedintensive,eos, old_cells, i, tess, time, tracers);
	      else
		cells[static_cast<size_t>(i)] = res.first;
	      bres[static_cast<size_t>(i)] = res.second;
	    }
	  else
	    cells[static_cast<size_t>(i)]=CellsEvolve[static_cast<size_t>(i)]->UpdatePrimitive
	      (conservedintensive,
	       eos,old_cells,i,tess,time,tracers);
	}
      catch(UniversalError& eo)
	{
	  eo.AddEntry("x momentum per unit volume",conservedintensive[static_cast<size_t>(i)].Momentum.x);
	  eo.AddEntry("y momentum per unit volume",conservedintensive[static_cast<size_t>(i)].Momentum.y);
	  eo.AddEntry("thermal energy per unit mass",conservedintensive[static_cast<size_t>(i)].Energy);
	  eo.AddEntry("Cell volume",tess.GetVolume(i));
	  eo.AddEntry("Cell x location",tess.GetMeshPoint(i).x);
	  eo.AddEntry("Cell y location",tess.GetMeshPoint(i).y);
#ifdef RICH_MPI
	  eo.AddEntry("Error in CPU",static_cast<double>(get_mpi_rank()));
#endif
	  update_primitives_rethrow(i,eo);
	}
    }
  return bres;
}

Primitive RotatePrimitive(Vector2D const& normaldir,
			  Vector2D const& paraldir,
			  Primitive const& p)
{
  Primitive res = p;
  res.Velocity.Set(Projection(p.Velocity,normaldir),
		   Projection(p.Velocity,paraldir));
  return res;
}

Conserved RotateFluxBack(Conserved const& c,
			 Vector2D const& normaldir,
			 Vector2D const& paraldir)
{
  Conserved res = c;
  res.Momentum = c.Momentum.x*normaldir/abs(normaldir)+
    c.Momentum.y*paraldir/abs(paraldir);
  return res;
}

Conserved FluxInBulk(Vector2D const& normaldir,
		     Vector2D const& paraldir,
		     Primitive const& left,
		     Primitive const& right,
		     Vector2D const& edge_velocity,
		     RiemannSolver const& rs)
{
  const Primitive rotated_left = RotatePrimitive(normaldir, paraldir, left);
  const Primitive rotated_right = RotatePrimitive(normaldir, paraldir, right);
  const double normal_speed = Projection(edge_velocity,normaldir);
  const Conserved res = rs(rotated_left, rotated_right, normal_speed);
  return RotateFluxBack(res, normaldir, paraldir);
}

void ExternalForceContribution
(const Tessellation& tess,
 const PhysicalGeometry& pg,
 const CacheData& cd,
 const vector<ComputationalCell>& cells,
 const vector<Extensive>& fluxes,
 const vector<Vector2D>& point_velocities,
 const SourceTerm& force,
 double t,
 double dt,
 vector<Extensive>& extensives)
{
  const vector<Extensive> diff = force
    (tess,pg,cd,cells,fluxes,point_velocities,t);
  for(size_t i=0;i<static_cast<size_t>(tess.GetPointNo());++i){
    extensives[i].mass += dt*diff[i].mass;
    extensives[i].momentum += dt*diff[i].momentum;
    extensives[i].energy += dt*diff[i].energy;
    for(boost::container::flat_map<std::string,double>::const_iterator it =
	  diff[i].tracers.begin();
	it!=diff[i].tracers.end(); ++it)
      extensives[i].tracers[it->first] += dt*it->second;
  }
}

vector<Vector2D> get_all_mesh_points
(Tessellation const& tess)
{
  vector<Vector2D> res(static_cast<size_t>(tess.GetPointNo()));
  for(int i=0;i<static_cast<int>(tess.GetPointNo());++i)
    res[static_cast<size_t>(i)] = tess.GetMeshPoint(i);
  return res;
}

vector<CustomEvolution*> convert_indices_to_custom_evolution
(const CustomEvolutionManager& cem,const vector<size_t>& indices)
{
  vector<CustomEvolution*> res(indices.size());
  for(size_t i=0;i<res.size();++i)
    res[i] = cem.getFunction(indices[i]);
  return res;
}

vector<Primitive> make_eos_consistent
(vector<Primitive> const& vp,
 EquationOfState const& eos)
{
  vector<Primitive> res = vp;
  for(int i=0;i<static_cast<int>(vp.size());++i)
    res[static_cast<size_t>(i)] = make_eos_consistent(vp[static_cast<size_t>(i)],eos);
  return res;
}

vector<double> GetMaxKineticEnergy(Tessellation const& tess,vector<Primitive> const&
				   cells,vector<CustomEvolution*> const& /*customevolve*/)
{
  const int n=tess.GetPointNo();
  vector<double> res;
  res.resize(static_cast<size_t>(n));
  for(int j=0;j<n;++j)
    {
      vector<int> neightemp=tess.GetNeighbors(j);
      vector<int> neigh;
      for(size_t i=0;i<neightemp.size();++i)
	if(neightemp[i]>=0)
	  neigh.push_back(neightemp[i]);
      double e=pow(abs(cells[static_cast<size_t>(j)].Velocity-cells[static_cast<size_t>(neigh[0])].Velocity),2);
      for(int i=1;i<static_cast<int>(neigh.size());++i)
	{// This could be made much faster by writing the expression implicitly
	  e=max(e,pow(abs(cells[static_cast<size_t>(j)].Velocity-cells[static_cast<size_t>(neigh[static_cast<size_t>(i)])].Velocity),2));
	}
      res[static_cast<size_t>(j)]=0.5*e;
    }
  return res;
}

vector<double> GetForceEnergy(Tessellation const& tess,
			      vector<double> const& g)
{
  vector<double> res;
  int n=int(g.size());
  res.resize(static_cast<size_t>(n));
  for(int i=0;i<n;++i)
    res[static_cast<size_t>(i)]=g[static_cast<size_t>(i)]*tess.GetWidth(i);
  return res;
}

void FixPressure(vector<Conserved> &intensive,vector<vector<double> > const& entropy,
		 EquationOfState const& eos,vector<double> const& Ek,
		 vector<double> const& Ef,double as,double bs,vector<CustomEvolution*>
		 const& customevolve,Tessellation const& tess,//vector<Conserved> &extensive,
		 vector<char> const& shockedcells,bool densityfloor)
{
  int n=tess.GetPointNo();
  double Et,Ek2;
  double temp;
  for(int i=0;i<n;++i)
    {
      if(customevolve[static_cast<size_t>(i)]==0||customevolve[static_cast<size_t>(i)]->TimeStepRelevant())
	{
	  if (intensive[static_cast<size_t>(i)].Mass < 0)
	    continue;
	  //Make intensive
	  temp=entropy[static_cast<size_t>(i)][0]/(tess.GetVolume(i)*intensive[static_cast<size_t>(i)].Mass);
	  Ek2=0.5*pow(abs(intensive[static_cast<size_t>(i)].Momentum)/intensive[static_cast<size_t>(i)].Mass,2);
	  Et=intensive[static_cast<size_t>(i)].Energy/intensive[static_cast<size_t>(i)].Mass-Ek2;
	  if((Et<as*Ek[static_cast<size_t>(i)])||(Et<bs*Ef[static_cast<size_t>(i)]))
	    {
	      if((shockedcells[static_cast<size_t>(i)]==0)||Et<0)
		{
		  Et=eos.dp2e(intensive[static_cast<size_t>(i)].Mass,
			      eos.sd2p(temp,intensive[static_cast<size_t>(i)].Mass));
		  if(Et<0&&!densityfloor)
		    {
		      UniversalError eo("Negative thermal enegry");
		      eo.AddEntry("Cell index",i);
		      eo.AddEntry("Thermal energy",Et);
		      eo.AddEntry("ShockedStatus",shockedcells[static_cast<size_t>(i)]);
		      eo.AddEntry("Extensive entropy",entropy[static_cast<size_t>(i)][0]);
		      eo.AddEntry("The density",intensive[static_cast<size_t>(i)].Mass);
		      throw eo;
		    }
		  intensive[static_cast<size_t>(i)].Energy=intensive[static_cast<size_t>(i)].Mass*(Et+Ek2);
		  //extensive[i].Energy=tess.GetVolume(i)*intensive[i].Energy;
		}
	    }
	}
    }
}

bool NearBoundary(int index,Tessellation const& tess,
		  vector<CustomEvolution*> const& /*customevolve*/)
{
  vector<int> neigh=tess.GetNeighbors(index);
  int n=int(neigh.size());
  for(int i=0;i<n;++i)
    {
      if(neigh[static_cast<size_t>(i)]<0)
	return true;
      /*if(customevolve[neigh[i]]!=0)
	return true;*/
    }
  return false;
}

namespace {
  vector<double> scalar_mult(const vector<double>& v,
			     double s)
  {
    if(v.empty())
      return vector<double>();
    vector<double> res(v.size());
    for(size_t i=0;i<v.size();++i)
      res[i] = s*v[i];
    return res;
  }
}

namespace {
  class ExtensiveTracerCalculator: public LazyList<vector<double> >
  {
  public:

    ExtensiveTracerCalculator(const vector<vector<double> >& tracers,
			      const Tessellation& tess,
			      const vector<Primitive>& cells,
			      const PhysicalGeometry& pg):
      tracers_(tracers), tess_(tess), cells_(cells), pg_(pg) {}

    size_t size(void) const
    {
      if(tracers_.empty())
	return 0;
      else
	return static_cast<size_t>(tess_.GetPointNo());
    }

    vector<double> operator[](size_t i) const
    {
      return scalar_mult
	(tracers_[i],
	 pg_.calcVolume(serial_generate(CellEdgesGetter(tess_,static_cast<int>(i))))*
	 cells_[i].Density);
    }

  private:
    const vector<vector<double> >& tracers_;
    const Tessellation& tess_;
    const vector<Primitive>& cells_;
    const PhysicalGeometry& pg_;
  };
}

vector<vector<double> > calc_extensive_tracer
(const vector<vector<double> >& intensive_tracer,
 const Tessellation& tess,
 const vector<Primitive>& cells,
 const PhysicalGeometry& pg)
{
  return serial_generate(ExtensiveTracerCalculator(intensive_tracer,
						   tess,
						   cells,
						   pg));
}

void MakeTracerExtensive(vector<vector<double> > const &tracer,
			 Tessellation const& tess,
			 vector<Primitive> const& cells,
			 vector<vector<double> > &result)
{
  const size_t n=static_cast<size_t>(tess.GetPointNo());
  result.resize(n);
  for(size_t i=0;i<n;++i)
    result[i] = scalar_mult(tracer[i],
			    tess.GetVolume(static_cast<int>(i))*cells[i].Density);
}

namespace {
  vector<double> scalar_div(const vector<double>& v,
			    const double s)
  {
    vector<double> res(v.size());
    for(size_t i=0;i<res.size();++i)
      res[i] = v[i]/s;
    return res;
  }

  class IntensiveTracerCalculator: public LazyList<vector<double> >
  {
  public:

    IntensiveTracerCalculator(const vector<vector<double> >& extensive,const Tessellation& tess,
			      const vector<Primitive>& cells,const PhysicalGeometry& pg,
			      const vector<vector<double> > &old_trace_intensive,const vector<bool> &min_density,
			      const vector<CustomEvolution*> &cevolve) :
      extensive_(extensive),tess_(tess), cells_(cells), pg_(pg),
      old_trace_intensive_(old_trace_intensive), min_density_(min_density), cevolve_(cevolve){}

    size_t size(void) const
    {
      if(extensive_.empty())
	return 0;
      else
	return static_cast<size_t>(tess_.GetPointNo());
    }

    vector<double> operator[](size_t i) const
    {
      if (min_density_[i] && !cevolve_[i])
	return old_trace_intensive_[i];
      else
	{
	  const double mass = cells_[i].Density*
	    pg_.calcVolume(serial_generate(CellEdgesGetter(tess_, static_cast<int>(i))));
	  return scalar_div(extensive_[i], mass);
	}
    }

  private:
    const vector<vector<double> >& extensive_;
    const Tessellation& tess_;
    const vector<Primitive>& cells_;
    const PhysicalGeometry& pg_;
    const vector<vector<double> > &old_trace_intensive_;
    const vector<bool> &min_density_;
    const vector<CustomEvolution*> &cevolve_;
  };
}

void MakeTracerIntensive(vector<vector<double> > &tracer,
			 const vector<vector<double> >& extensive,
			 const Tessellation& tess,
			 const vector<Primitive>& cells,
			 const PhysicalGeometry& pg,vector<bool> const& min_density_on,vector<vector<double> > const& old_trace,
			 vector<CustomEvolution*> const& cevolve)
{
  tracer = serial_generate(IntensiveTracerCalculator(extensive,tess,cells,pg,old_trace,min_density_on,cevolve));
}

void UpdateTracerExtensive(vector<vector<double> > &tracerextensive,
			   vector<vector<double> > const& tracerchange,vector<CustomEvolution*> const&
			   CellsEvolve,vector<Primitive> const& cells,Tessellation const& tess,
			   double time)
{
  for(size_t i=0;i<tracerextensive.size();++i)
    if(CellsEvolve[i])
      tracerextensive[i]=CellsEvolve[i]->UpdateTracer
	(static_cast<int>(i),tracerextensive,tracerchange,cells,tess,time);
    else
      for(size_t j=0;j<tracerextensive[i].size();++j)
	tracerextensive[i][j]+=tracerchange[i][j];
}

void TracerResetCalc
(double alpha,SpatialDistribution const& originalD,
 SpatialDistribution const& originalP,SpatialDistribution const& originalVx,
 SpatialDistribution const& originalVy,vector<SpatialDistribution const*> const& originalTracers,vector<Primitive> &cells,
 Tessellation const& tess,vector<vector<double> > &tracer,
 int tracerindex,EquationOfState const& eos,vector<CustomEvolution*>
 const& cevolve,bool coldflows)
{
  const int n = tess.GetPointNo();
  if(n<1)
    return;
  Vector2D velocity;
  if(tracer.empty())
    return;
  if(tracerindex>=static_cast<int>(tracer[0].size())||tracerindex<0)
    throw UniversalError("Error in tracerReset, wrong dimension for tracer");
  for(int i=0;i<n;++i)
    {
      bool customforce=false;
      if(cevolve[static_cast<size_t>(i)])
	customforce=cevolve[static_cast<size_t>(i)]->ShouldForceTracerReset();
      if((tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]<alpha)||customforce)
	{
	  velocity.Set(originalVx(tess.GetCellCM(i)),
		       originalVy(tess.GetCellCM(i)));
	  cells[static_cast<size_t>(i)]=CalcPrimitive(originalD(tess.GetCellCM(i)),
						      originalP(tess.GetCellCM(i)),velocity,eos);
	  if(tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]<0)
	    tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]=0;
	  for (size_t j = 0;j<tracer[static_cast<size_t>(i)].size();++j)
	    if (coldflows&&j == 0)
	      tracer[static_cast<size_t>(i)][j] = eos.dp2s(cells[static_cast<size_t>(i)].Density, cells[static_cast<size_t>(i)].Pressure);
	    else
	      if (static_cast<int>(j) != tracerindex)
		tracer[static_cast<size_t>(i)][j] = originalTracers[j]->operator()(tess.GetCellCM(i));
	}
    }
  return;
}

void GetPointToRemove(Tessellation const& tess,Vector2D const& point,
		      double R,vector<int> & PointToRemove,int Inner)
{
  int n=tess.GetPointNo();
  PointToRemove.clear();
  for(int i=Inner;i<n;++i)
    {
      // Check if point is completly engulfed
      bool test=true;
      vector<int> neigh=tess.GetNeighbors(i);
      for(int j=0;j<static_cast<int>(neigh.size());++j)
	if(neigh[static_cast<size_t>(j)]>=Inner)
	  test=false;
      // Is point inside a radius?
      if(abs(point-tess.GetMeshPoint(i))<R||test)
	PointToRemove.push_back(i);
    }
  return;
}

void FixAdvection(vector<Conserved>& extensive,
		  vector<Conserved> const& intensive,Tessellation const& tessold,
		  Tessellation const& tessnew,vector<Vector2D> const& facevelocity,
		  double dt,vector<Vector2D> const& /*pointvelocity*/)
{
  int n=tessold.GetTotalSidesNumber();
  int npoints=tessold.GetPointNo();
  vector<double> Rold(static_cast<size_t>(npoints)),Rnew(static_cast<size_t>(npoints));
  vector<vector<Vector2D> > pold(static_cast<size_t>(npoints)),pnew(static_cast<size_t>(npoints));
  for(int i=0;i<npoints;++i)
    {
      Rold[static_cast<size_t>(i)]=tessold.GetWidth(i);
      Rnew[static_cast<size_t>(i)]=tessnew.GetWidth(i);
      ConvexHull(pold[static_cast<size_t>(i)],&tessold,i);
      ConvexHull(pnew[static_cast<size_t>(i)],&tessnew,i);
    }

  PolygonOverlap polyoverlap;
  double eps=1e-7;
  for(int i=0;i<n;++i)
    {
      Edge const& edge=tessold.GetEdge(i);
      int n0=edge.neighbors.first;
      int n1=edge.neighbors.second;
      if(n0<0||n1<0)
	continue;
      Vector2D norm(tessold.GetMeshPoint(n1)-tessold.GetMeshPoint(n0));
      norm=norm/abs(norm);
      norm=norm*edge.GetLength();
      double dv_dt=ScalarProd(facevelocity[static_cast<size_t>(i)],norm)*dt;
      /*		vector<Vector2D> poly0,poly1;
			ConvexHull(poly0,&tessold,tessold.GetOriginalIndex(n0));
			ConvexHull(poly1,&tessnew,tessold.GetOriginalIndex(n1));
			if(n0>=npoints)
			{
			const Vector2D diff(tessold.GetMeshPoint(tessold.GetOriginalIndex(n0))
			-tessold.GetMeshPoint(n0));
			int N=static_cast<int>(poly0.size());
			for(int j=0;j<N;++j)
			poly0[j]-=diff;
			}
			if(n1>=npoints)
			{
			const Vector2D diff(tessnew.GetMeshPoint(tessnew.GetOriginalIndex(n1))
			-tessnew.GetMeshPoint(n1));
			int N=static_cast<int>(poly1.size());
			for(int j=0;j<N;++j)
			poly1[j]-=diff;
			}*/
      double real_dv1=polyoverlap.polygon_overlap_area
	    (pold[static_cast<size_t>(tessold.GetOriginalIndex(n0))],
	     pnew[static_cast<size_t>(tessold.GetOriginalIndex(n1))],
	     Rold[static_cast<size_t>(tessold.GetOriginalIndex(n0))]*eps,
	     Rnew[static_cast<size_t>(tessold.GetOriginalIndex(n1))]*eps);
      double real_dv0=polyoverlap.polygon_overlap_area
	    (pnew[static_cast<size_t>(tessold.GetOriginalIndex(n0))],
	     pold[static_cast<size_t>(tessold.GetOriginalIndex(n1))],
	     Rnew[static_cast<size_t>(tessold.GetOriginalIndex(n0))]*eps,
	     Rold[static_cast<size_t>(tessold.GetOriginalIndex(n1))]*eps);

      if(dv_dt>0)
	{
	  if(n0<npoints)
	    {
	      extensive[static_cast<size_t>(n0)]+=(real_dv0-dv_dt)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	      extensive[static_cast<size_t>(n0)]-=real_dv1*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	    }
	  if(n1<npoints)
	    {
	      extensive[static_cast<size_t>(n1)]+=(dv_dt-real_dv0)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	      extensive[static_cast<size_t>(n1)]+=real_dv1*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	    }
	}
      else
	{
	  if(n0<npoints)
	    {
	      extensive[static_cast<size_t>(n0)]-=(real_dv1+dv_dt)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	      extensive[static_cast<size_t>(n0)]+=real_dv0*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	    }
	  if(n1<npoints)
	    {
	      extensive[static_cast<size_t>(n1)]-=(-dv_dt-real_dv1)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
	      extensive[static_cast<size_t>(n1)]-=real_dv0*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
	    }
	}
    }
}

double determine_time_step(double hydro_time_step,
			   double external_dt,
			   double current_time,
			   double end_time)
{
  double dt = hydro_time_step;
  if(external_dt>0)
    dt = std::min(external_dt,dt);
  if(end_time>0)
    dt = std::min(end_time-current_time,dt);

#ifdef RICH_MPI
  double dt_temp = dt;
  MPI_Reduce(&dt_temp,&dt,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

  return dt;
}
