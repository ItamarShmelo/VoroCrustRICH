/*! \file AMR3D.hpp
\brief Abstract class for amr in 3D
\author Elad Steinberg
*/

#ifndef AMR3D_HPP
#define AMR3D_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"
#include "../common/equation_of_state.hpp"
#include "3D/tesselation/Tessellation3D.hpp"
#include <boost/scoped_ptr.hpp>
#include "LinearGauss3D.hpp"
#include "hdsim_3d.hpp"

//! \brief Abstract class for cell update scheme in amr
class AMRCellUpdater3D
{
public:

	/*! \brief Calculates the computational cell
	\param extensive Extensive conserved variables
	\param eos Equation of state
	\param volume Cell volume
	\param old_cell Old computational cell
	\param tracerstickernames The names of the tracers and stickers
	\return Computational cell
	*/
	virtual ComputationalCell3D ConvertExtensiveToPrimitve3D(Conserved3D& extensive, const EquationOfState& eos,
		double volume, ComputationalCell3D const& old_cell) const = 0;

	//! \brief Class destructor
	virtual ~AMRCellUpdater3D(void);

  /*! \brief Null constructor
   */
  AMRCellUpdater3D(void);

  /*! \brief Copy constructor
   */
  AMRCellUpdater3D(const AMRCellUpdater3D&);

  /*! \brief Copy constructor
    \return Reference to new object
   */
  AMRCellUpdater3D& operator=(const AMRCellUpdater3D&);
};

//! \brief Abstract class for extensive update scheme in amr
class AMRExtensiveUpdater3D
{
public:

	/*! \brief Calculates the computational cell
	\param cell Computational cell
	\param eos Equation of state
	\param volume Cell volume
	\param tracerstickernames The names of the tracers and stickers
	\param slope Gradients
	\param CMold Old centre of mass
	\param CMnew New centre of mass
	\return Extensive
	*/
	virtual Conserved3D ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
		double volume, Slope3D const& slope,Vector3D const& CMold,Vector3D const& CMnew) const = 0;

	//! \brief Class destructor
	virtual ~AMRExtensiveUpdater3D(void);

  /*! \brief Null constructor
   */
  AMRExtensiveUpdater3D(void);

  /*! \brief Copy constructor
   */
  AMRExtensiveUpdater3D(const AMRExtensiveUpdater3D&);

  /*! \brief Copy assignment
    \return Reference to new object
   */
  AMRExtensiveUpdater3D& operator=(const AMRExtensiveUpdater3D&);
};

//! \brief Simple class for extensive update scheme in amr
class SimpleAMRExtensiveUpdater3D : public AMRExtensiveUpdater3D
{
public:
	SimpleAMRExtensiveUpdater3D(void);

	Conserved3D ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
		double volume, Slope3D const& slope, Vector3D const& CMold, Vector3D const& CMnew) const override;
};

//! \brief Simple class for cell update scheme in amr
class SimpleAMRCellUpdater3D : public AMRCellUpdater3D
{
private:
	const vector<string>  toskip_;
public:
	/*!
	\brief class constructor
	\param toskip A list of sticker names to skip their cell update
	*/
	SimpleAMRCellUpdater3D(const vector<string>& toskip = vector<string>());

	ComputationalCell3D ConvertExtensiveToPrimitve3D(Conserved3D& extensive, const EquationOfState& eos,
		double volume, ComputationalCell3D const& old_cell) const override;
};

//! \brief Simple class for extensive update scheme in amr for SR
class SimpleAMRExtensiveUpdaterSR3D : public AMRExtensiveUpdater3D
{
public:
  /*! 
    \param cell Computational cell
    \param eos Equation of state
    \param volume Volume
    \param tracerstickernames Tracers and stickers names
    \param slope Gradients
    \param CMold Old centre of mass
    \param CMnew New centre of mass
    \return Conserved variables
   */
  Conserved3D ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
		double volume, Slope3D const& slope, Vector3D const& CMold, Vector3D const& CMnew) const override;
};

//! \brief Simple class for cell update scheme in amr for SR
class SimpleAMRCellUpdaterSR3D : public AMRCellUpdater3D
{
private:
	const double G_;
	const vector<string>  toskip_;
public:
	/*!
	\brief class constructor
	\param G The adiabatic index
	\param toskip A list of sticker names to skip their cell update
	*/
  SimpleAMRCellUpdaterSR3D(double G, const vector<string>& toskip);

	ComputationalCell3D ConvertExtensiveToPrimitve3D(Conserved3D& extensive, const EquationOfState& eos,
							 double volume, const ComputationalCell3D& old_cell) const override;
};

//! \brief Chooses which cells should be remove
class CellsToRemove3D
{
public:
	/*!
	\brief Finds the cells to remove
	\param tess The tesselation
	\param cells The computational cells
	\param time The sim time
	\param tracerstickernames The names of the tracers and stickers
	\return The indeces of cells to remove with a corresponding merit which decides if there are neighboring cells which one to choose to remove
	*/
	virtual std::pair<vector<size_t>, vector<double> > ToRemove(Tessellation3D const& tess,
		vector<ComputationalCell3D> const& cells, double time)const = 0;

	//! \brief Virtual destructor
	virtual ~CellsToRemove3D(void);
};

//! \brief Chooses which cells should be refined
class CellsToRefine3D
{
public:
	/*!
	\brief Finds the cells to refine
	\param tess The tesselation
	\param cells The computational cells
	\param time The sim time
	\param tracerstickernames The names of the tracers and stickers
	\return The indeces of cells to remove and the direction to split (can be given empty)
	*/
	virtual std::pair<vector<size_t>,vector<Vector3D> > ToRefine(Tessellation3D const& tess,
		vector<ComputationalCell3D> const& cells, double time)const = 0;

	//! \brief Virtual destructor
	virtual ~CellsToRefine3D(void);
};

//! \brief Base class for amr
class AMR3D
{
private:
	EquationOfState const& eos_;
	CellsToRefine3D const& refine_;
	CellsToRemove3D  const& remove_;
	SimpleAMRCellUpdater3D scu_;
	SimpleAMRExtensiveUpdater3D seu_;
	SpatialReconstruction3D &interp_;
	AMRCellUpdater3D* cu_;
	AMRExtensiveUpdater3D* eu_;
	AMR3D(AMR3D const& amr);
	AMR3D& operator=(AMR3D const&);
	
public:
	/*!
	\brief Runs the AMR
	\param sim The sim object
	*/
	void operator() (HDSim3D &sim);

	/*! \brief Class constructor
	\param refine Refinement scheme
	\param remove Removal scheme
	\param cu Cell updater
	\param eu Extensive updater
	\param eos Equation of state
	\param interp Interpolation scheme
	*/
	AMR3D(EquationOfState const& eos, CellsToRefine3D const& refine, CellsToRemove3D const& remove,SpatialReconstruction3D &interp, AMRCellUpdater3D* cu = nullptr,
		AMRExtensiveUpdater3D* eu = nullptr);
	//! Class destructor
	~AMR3D();
};

#endif // AMR3D_HPP
