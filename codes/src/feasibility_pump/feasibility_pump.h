#pragma once

#include <ilcplex/ilocplex.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <utility>
#include <list>
#include "src/instance.h"
#include "src/heuristic_solution.h"
#include "src/benders_generic_callback.h"

class FeasibilityPump
{
public:
	FeasibilityPump();
	~FeasibilityPump();
	void Init(Instance &);
	void Reset();
	void Run();
	FPHeuristicSolution solution_;

private:
	IloEnv env_;	 // Cplex environment
	IloCplex cplex_; // Cplex solver
	IloModel model_; // Cplex model
	IloObjective objective_;
	IloExpr original_obj_expr_;
	double original_obj_norm_;
	double curr_alpha_; /* ratio by which the new obj is combined from distance and original obj */
	double previous_alpha_;
	double previous_integrality_gap_;
	double curr_integrality_gap_;

	IloNumVarArray f_;
	IloNumVar slack_;
	MasterVariables master_vars_;

	boost::dynamic_bitset<> curr_int_y_;
	boost::dynamic_bitset<> previous_int_y_;
	boost::dynamic_bitset<> curr_int_x_;
	boost::dynamic_bitset<> previous_int_x_;
	IloNumArray curr_relax_y_;
	IloNumArray curr_relax_x_;

	std::list<std::pair<int, double>> curr_integrality_gaps_x_;
	std::list<std::pair<int, double>> curr_integrality_gaps_y_;

	bool found_int_y_;
	bool found_int_x_;

	Instance *curr_instance_;

	void RetrieveAndRoundVertexValues();
	void RetrieveAndRoundArcValuesBasic();
	void RetrieveAndRoundArcValuesBFS();
	void BuildPathsBFS(std::vector<std::list<std::pair<int, double>>> &residual_graph, std::list<std::list<int>> &paths);
	bool BuildPathsBFSIter(std::vector<std::list<std::pair<int, double>>> &residual_graph, std::list<std::list<int>> &paths,
						   std::list<int> curr_path, std::vector<char> &colors, size_t vertex);
	void BuildModel(double *, double *);
	void BuildHeuristicSolution();

	void SetNewObjStage1();
	void SetNewObjStage2();
};
