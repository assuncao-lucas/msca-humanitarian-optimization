#ifndef FOMRULATIONS_H_
#define FORMULATIONS_H_

#include <vector>
#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/matrix.hpp"
#include "src/general.h"
#include "src/heuristic_solution.h"


int a_var_to_index(int vertex, int budget, int num_vertices);

// dual problem variables.
int u_0_var_to_index(int vertex, int budget, int num_vertices);
int u_1_var_to_index(int vertex, int budget, int num_vertices);
int u_3_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs);
int u_4_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs);

std::pair<int,int> index_to_a_var(int index, int num_vertices);

template <class T> void SetSolutionStatus(IloCplex & cplex, Solution<T> & solution, bool solve_relax)
{
  if(solve_relax)
  {
    if((cplex.getCplexStatus() == IloCplex::Infeasible)||(cplex.getCplexStatus() == IloCplex::InfOrUnbd)) solution.lp_ = -1;
    else solution.lp_ = cplex.getObjValue();
  }else{
    if((cplex.getCplexStatus() == IloCplex::Infeasible)||(cplex.getCplexStatus() == IloCplex::InfOrUnbd)) solution.is_feasible_ = false;
    else
    {
      double cost = cplex.getObjValue();
      if(cplex.getCplexStatus() == IloCplex::Optimal)
      {
        solution.lb_ = cost;
        solution.ub_ = cost;
        solution.is_optimal_ =  true;
      }else
      {
        solution.lb_ = cost;
        solution.ub_ = cplex.getBestObjValue();
        if(!(double_less(solution.lb_, solution.ub_))) solution.is_optimal_ = true;
      }
    }
  }
  solution.num_nodes_ = cplex.getNnodes();
}

void setPriorityOrder(IloCplex & cplex, IloEnv & env, IloNumVarArray & x, Instance & instance, double * R0);
Solution<int>* optimize(IloCplex & cplex, IloEnv& env, IloModel& model, IloNumVarArray & y, IloNumVarArray & x, Instance& instance, double total_time_limit, bool solve_relax, bool callback, bool find_root_cuts, double * R0, double * Rn, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol);

static void PopulateByRowCommon(IloEnv& env, IloModel& model, IloNumVar & slack, IloNumVarArray & y, IloNumVarArray & x, Instance& instance, bool force_use_all_vehicles);

static void PopulateByRowCompactBaselineContinuousSpace(IloEnv& env, IloModel& model, IloNumVarArray & x, IloNumVarArray & y, IloNumVarArray &a, std::optional<const int *> x_values, std::optional<const int *> y_values, Instance& instance);
static void PopulateByRowDualCompactBaselineContinuousSpace(IloEnv& env, IloModel& model, IloNumVarArray &u_0, IloNumVarArray &u_1, IloNumVarArray &u_2, IloNumVarArray &u_3, std::optional<const int *> x_values, std::optional<const int *> y_values, Instance& instance);

static void PopulateByRowCompactBaseline(IloEnv& env, IloModel& model, IloNumVar & slack, IloNumVarArray & y, IloNumVarArray & x, IloNumVarArray & a, Instance& instance, double* R0, double * Rn, bool force_use_all_vehicles);
Solution<int>* CompactBaseline(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles);
Solution<double>* PrimalSubproblemCompactBaseline(Instance& inst, const int * x, const int * y, double * R0, double * Rn, double time_limit);
Solution<double>* DualSubproblemCompactBaseline(Instance& inst, const int * x, const int * y, double * R0, double * Rn, double time_limit);

static void PopulateByRowCompactSingleCommodity(IloEnv& env, IloModel& model, IloNumVar & slack, IloNumVarArray & y, IloNumVarArray & x, IloNumVarArray & f, Instance& instance, double* R0, double * Rn, bool force_use_all_vehicles);
Solution<int>* CompactSingleCommodity(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles);

#endif