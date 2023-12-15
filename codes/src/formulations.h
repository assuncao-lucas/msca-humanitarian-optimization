#ifndef FOMRULATIONS_H_
#define FORMULATIONS_H_

#include <vector>
#include <unordered_map>
#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/matrix.hpp"
#include "src/general.h"
#include "src/heuristic_solution.h"

struct DualVariables
{
  IloNumVarArray u_0;
  IloNumVarArray u_1;
  IloNumVarArray u_2;
  IloNumVarArray u_3;
  std::unordered_map<IloInt,IloExpr> ext_ray_coef;
};

int a_var_to_index(int vertex, int budget, int num_vertices);
int f_var_to_index(int arc_pos, int budget, int num_arcs);

// dual problem variables.
int u_0_var_to_index(int vertex, int budget, int num_vertices);
int u_1_var_to_index(int vertex, int budget, int num_vertices);
int u_2_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs);
int u_3_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs);

std::pair<int,int> index_to_a_var(int index, int num_vertices);

template <class T> void SetSolutionStatus(IloCplex & cplex, Solution<T> & solution, bool solve_relax)
{
  if(solve_relax)
  {
    if((cplex.getCplexStatus() == IloCplex::Infeasible)||(cplex.getCplexStatus() == IloCplex::InfOrUnbd)) solution.lp_ = -1;
    else solution.lp_ = cplex.getObjValue();
    if((cplex.getCplexStatus() == IloCplex::Optimal)||(cplex.getCplexStatus() == IloCplex::OptimalTol))
    {
      solution.is_optimal_ =  true;
    }
  }else{
    if((cplex.getCplexStatus() == IloCplex::Infeasible)||(cplex.getCplexStatus() == IloCplex::InfOrUnbd)) solution.is_feasible_ = false;
    else
    {
      double cost = cplex.getObjValue();
      if((cplex.getCplexStatus() == IloCplex::Optimal)||(cplex.getCplexStatus() == IloCplex::OptimalTol))
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
void SetPriorityOrder(IloCplex & cplex, IloEnv & env, IloNumVarArray & x, Instance & instance, double * R0);
static void AddNamesToDualVariables(DualVariables& dual_vars, Instance& inst);
static void SetDualVariablesProperties(IloEnv& master_env, DualVariables& dual_vars, IloNumVarArray& x, IloNumVarArray& y, Instance& instance);

static void PopulateByRowCommon(IloEnv& env, IloModel& model, IloNumVar & slack, IloNumVarArray & y, IloNumVarArray & x, Instance& instance, bool force_use_all_vehicles);
static void PopulateByRowCompactBaselineContinuousSpace(IloEnv& env, IloModel& model, IloNumVarArray & x, IloNumVarArray & y, IloNumVarArray &a, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, Instance& instance);
static void PopulateByRowDualCompactBaselineContinuousSpace(IloEnv& env, IloModel& model, DualVariables &dual_vars, Instance& instance);
static void PopulateByRowCompactBaseline(IloCplex &cplex, IloEnv& env, IloModel& model, IloNumVar & slack, IloNumVarArray & y, IloNumVarArray & x, IloNumVarArray & a, Instance& instance, double* R0, double * Rn, bool force_use_all_vehicles, bool export_model);

static void PopulateByRowCompactSingleCommodityContinuousSpace(IloEnv& env, IloModel& model, IloNumVarArray & x, IloNumVarArray & y, IloNumVarArray &f, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, double * R0, double * Rn, Instance& instance);
static void PopulateByRowCompactSingleCommodity(IloCplex& cplex, IloEnv& env, IloModel& model, IloNumVar& slack, IloNumVarArray & y, IloNumVarArray & x, IloNumVarArray & f, Instance& instance, double * R0, double* Rn, bool force_use_all_vehicles, bool export_model);

static void FillObjectiveExpressionDualCompactBaselineContinuousSpace(IloExpr& obj, DualVariables &dual_vars, IloNumArray& x_values, IloNumArray& y_values, Instance& instance);
static void AllocateDualVariables(IloEnv& env, DualVariables& dual_vars, Instance & instance);

Solution<int> * BendersCompactBaseline(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles, bool export_model);
IloBool SeparateBendersCut(IloEnv& master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar& dual_bound, IloNumArray & x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex& worker_cplex, DualVariables& dual_vars, Instance& instance, IloObjective& worker_obj, IloExpr & cut_expr);

Solution<int>* CompactBaseline(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles, bool export_model);
Solution<double>* PrimalSubproblemCompactBaseline(Instance& inst, IloNumArray& x_values, IloNumArray& y_values, double * R0, double * Rn, double time_limit, bool export_model);
Solution<double>* DualSubproblemCompactBaseline(Instance& inst, IloNumArray& x_values, IloNumArray& y_values, double * R0, double * Rn, double time_limit, bool export_model);
Solution<int>* CompactSingleCommodity(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles, bool export_model);

Solution<int>* optimize(IloCplex & cplex, IloEnv& env, IloModel& model, IloNumVarArray & y, IloNumVarArray & x, std::optional<std::reference_wrapper<IloNumVar>> dual_bound_opt, Instance& instance, double total_time_limit, bool solve_relax, bool callback, bool find_root_cuts, double * R0, double * Rn, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool export_model);
Solution<double> * optimizeLP(IloCplex & cplex, IloEnv& env, IloModel& model, Instance& instance, double total_time_limit, double * R0, double * Rn);

#endif