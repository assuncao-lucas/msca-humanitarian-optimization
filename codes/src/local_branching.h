#ifndef LOCAL_BRANCHING_H
#define LOCAL_BRANCHING_H

#include <ilcplex/ilocplex.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <utility>
#include <list>
#include "instance.h"
#include "heuristic_solution.h"
#include "user_cut.h"

class LocalBranching
{
public:
	LocalBranching();
	~LocalBranching();
	void Init(Instance& instance);
	void Reset();
	boost::dynamic_bitset<> * Run(double time_limit, HeuristicSolution * initial_sol);
	LBHeuristicSolution solution_;
private:
	IloEnv env_;					// Cplex environment
	IloCplex cplex_;				// Cplex solver
	IloModel model_;
	IloRange lb_constraint_;

	IloNumVarArray f;
	IloNumVarArray x;
	IloNumVarArray y;
    IloNumVar slack_;

	Instance * curr_instance_;

	//void BuildLeftBranchOnX(IloModel&, double *, double*,int);

	void BuildLeftLBConstraintOnX(boost::dynamic_bitset<>& bitset, int k);
	//void BuildRightLBConstraintOnX(IloModel& model,boost::dynamic_bitset<>& bitset, int k);

	//void BuildLeftLBConstraintOnY(IloModel& model,boost::dynamic_bitset<>& bitset, int k);
	//void BuildRightLBConstraintOnY(IloModel& model,boost::dynamic_bitset<>& bitset, int k);
	void AddCutsAsConstraints(std::list<UserCutGeneral*> * initial_cuts);

	void BuildHeuristicSolution();

	void RetrieveAndRoundArcValuesBasic(boost::dynamic_bitset<>& bitset);
};

#endif
