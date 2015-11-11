#ifndef CONDITION_ACTION_SEQUENCE2_HPP
#define CONDITION_ACTION_SEQUENCE2_HPP 1

#include "condition_action_sequence.hpp"
#include "spatial_reconstruction.hpp"

//! \brief Second order flux calculator based on a series of conditions and actions
class ConditionActionSequence2 : public FluxCalculator
{
public:

	//! \brief Action taken to calculate flux
	class Action2
	{
	public:

		/*! \brief Calculates flux
		\param edge Interface between cells
		\param tess Tessellation
		\param cells Computational cells
		\param eos Equation of state
		\param aux Auxiliary variable for assymetric problems (true means the relevant cell is on the left side, false mean right)
		\param edge_values The interpolated values at the edge 
		\return Flux
		*/
		virtual Extensive operator()
			(const Edge& edge,
				const Tessellation& tess,
				const Vector2D& edge_velocity,
				const vector<ComputationalCell>& cells,
				const EquationOfState& eos,
				const bool aux,
				const pair<ComputationalCell,ComputationalCell> & edge_values) const = 0;

		virtual ~Action2(void);
	};

	/*! \brief Class constructor
	\param sequence Series of condition and action action pairs
	*/
	ConditionActionSequence2
		(const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence,
			const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence2::Action2*> >& sequence2,
			SpatialReconstruction const& interp);

	~ConditionActionSequence2(void);

	vector<Extensive> operator()
		(const Tessellation& tess,
			const vector<Vector2D>& edge_velocities,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& extensives,
			const CacheData& cd,
			const EquationOfState& eos,
			const double time,
			const double dt) const;

private:
	const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> > sequence_;
	const vector<pair<const ConditionActionSequence::Condition*, const Action2*> > sequence2_;
	const SpatialReconstruction & interp_;
};

//! \brief Calculates flux between two regular bulk cells
class RegularFlux2 : public ConditionActionSequence2::Action2
{
public:

	/*! \brief Class constructor
	\param rs Riemann solver
	*/
	RegularFlux2(const RiemannSolver& rs);

	Extensive operator()
		(const Edge& edge,
			const Tessellation& tess,
			const Vector2D& edge_velocity,
			const vector<ComputationalCell>& cells,
			const EquationOfState& eos,
			const bool aux,
			const pair<ComputationalCell, ComputationalCell> & edge_values) const;

private:

	const RiemannSolver& rs_;
};


//! \brief Calculates flux assuming rigid wall boundary conditions
class RigidWallFlux2 : public ConditionActionSequence2::Action2
{
public:

	/*! \brief Class constructor
	\param rs Riemann solver
	*/
	RigidWallFlux2(const RiemannSolver& rs);

	Extensive operator()
		(const Edge& edge,
			const Tessellation& tess,
			const Vector2D& edge_velocity,
			const vector<ComputationalCell>& cells,
			const EquationOfState& eos,
			const bool aux,
			const pair<ComputationalCell, ComputationalCell> & edge_values) const;

private:
	const RiemannSolver& rs_;
};

class Ratchet : public ConditionActionSequence2::Action2
{
public:

	/*! \brief Class constructor
	\param in If the ratchet allows inflow or outflow
	*/
  Ratchet(const RiemannSolver&rs, const bool in);

	Extensive operator()
		(const Edge& edge,
			const Tessellation& tess,
			const Vector2D& edge_velocity,
			const vector<ComputationalCell>& cells,
			const EquationOfState& eos,
			const bool aux,
			const pair<ComputationalCell, ComputationalCell> & edge_values) const;

private:
	const bool in_;
	const RigidWallFlux2 wall_;
	const FreeFlowFlux free_;
};

#endif // CONDITION_ACTION_SEQUENCE2_HPP
