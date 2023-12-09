#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <list>
#include <queue>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <optional>
#include "src/graph.h"
#include "src/formulations.h"
#include "src/timer.h"
#include "src/user_cut.h"
ILOSTLBEGIN

// Implementation class for the user-defined user cut callback.
// The function BendersUserCallback allows to add Benders' cuts as user cuts.
//
ILOUSERCUTCALLBACK7(BendersUserCutCallback, IloCplex&, worker_cplex, IloObjective&, worker_obj,
                           DualVariables&, dual_vars, IloNumVarArray&, y, IloNumVarArray&, x,
                           IloNumVar&, dual_bound, Instance&, instance)
{
  // Skip the separation if not at the end of the cut loop
  if ( !isAfterCutLoop() )
    return;

  IloEnv master_env = getEnv();
  // Get the current solution.
  IloNumArray x_values(master_env);
  IloNumArray y_values(master_env);
  IloNum dual_bound_value = getValue(dual_bound);
  getValues(y_values, y);
  getValues(x_values, x);

  // Benders' cut separation.
  IloExpr cut_expr(master_env);
  IloBool sep_status = SeparateBendersCut(master_env,x,y,dual_bound,x_values,y_values,dual_bound_value, worker_cplex, dual_vars, instance, worker_obj, cut_expr);
  if(sep_status)
  {
    //std::cout << "found new Cut" << std::endl;
    add(cut_expr >= 0).end();
  }

  // Free memory.
  cut_expr.end();
  x_values.end();
  y_values.end();

  return;

} // END BendersUserCutCallback

// Implementation class for the user-defined lazy constraint callback.
// The function BendersLazyCallback allows to add Benders' cuts as lazy constraints.
//
ILOLAZYCONSTRAINTCALLBACK7(BendersLazyCallback, IloCplex&, worker_cplex, IloObjective&, worker_obj,
                           DualVariables&, dual_vars, IloNumVarArray&, y, IloNumVarArray&, x,
                           IloNumVar&, dual_bound, Instance&, instance)
{
  IloEnv master_env = getEnv();
  // Get the current solution.
  IloNumArray x_values(master_env);
  IloNumArray y_values(master_env);
  IloNum dual_bound_value = getValue(dual_bound);
  getValues(y_values, y);
  getValues(x_values, x);

  // Benders' cut separation.
  IloExpr cut_expr(master_env);
  IloBool sep_status = SeparateBendersCut(master_env,x,y,dual_bound,x_values,y_values,dual_bound_value, worker_cplex, dual_vars, instance, worker_obj, cut_expr);
  if(sep_status)
  {
    //std::cout << "found new Cut" << std::endl;
    add(cut_expr >= 0).end();
  }

  // Free memory.
  cut_expr.end();
  x_values.end();
  y_values.end();

  return;

} // END BendersLazyCallback

// This routine separates Benders' cuts violated by the current x solution.
// Violated cuts are found by solving the worker LP
//
IloBool SeparateBendersCut(IloEnv& master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar& dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex& worker_cplex, DualVariables& dual_vars, Instance & instance, IloObjective& worker_obj, IloExpr& cut_expr)
{
  // std::cout << "try to separate benders cut" << std::endl; 
  IloBool violatedCutFound = IloFalse;
  IloEnv worker_env = worker_cplex.getEnv();
  IloModel worker_model = worker_cplex.getModel();

  // Update the objective function in the worker LP.
  worker_model.remove(worker_obj);
  IloExpr obj_expr = worker_obj.getExpr();
  obj_expr.clear();
  FillObjectiveExpressionDualCompactBaselineContinuousSpace(obj_expr,dual_vars,x_values,y_values,instance);
  worker_obj.setExpr(obj_expr);
  worker_model.add(worker_obj);
  obj_expr.end(); 
  AddNamesToDualVariables(dual_vars,instance);
  worker_cplex.extract(worker_model);
  //worker_cplex.exportModel("worker_model_Benders_compact_baseline.lp");

  // Solve the worker LP
  worker_cplex.solve();

  auto status = worker_cplex.getStatus();
  IloNum new_dual_bound = -IloInfinity;

  if (status == IloAlgorithm::InfeasibleOrUnbounded || status == IloAlgorithm::Infeasible )
    std::cout << "DEU RUIM!" << std::endl;

  // Get the violated cut as an unbounded ray of the worker LP
  if ( status == IloAlgorithm::Unbounded)
  {
    //std::cout << "found new feasibility cut" << std::endl;
    IloNumVarArray var(worker_env);
    IloNumArray val(worker_env);

    worker_cplex.getRay(val, var);

    for(IloInt i = 0; i < val.getSize(); ++i)
    {
      IloInt var_id = var[i].getId();
      auto var_pos = dual_vars.ext_ray_coef.find(var_id);
      if ( var_pos != dual_vars.ext_ray_coef.end())
        cut_expr += operator*(val[i],dual_vars.ext_ray_coef[var_id]);
    }

    var.end();
    val.end();

    violatedCutFound = IloTrue;
  }

  // if (status == IloAlgorithm::Optimal)
  // {
  //   new_dual_bound = worker_cplex.getObjValue();
  //   std::cout << status << " " << new_dual_bound << " x " << dual_bound_value << std::endl;
  // }

  // if the problem is solved to optimality and the current master solution is violated, a new optimality cut is found.
  if ( status == IloAlgorithm::Optimal && double_greater(dual_bound_value,new_dual_bound) )
  {
    //std::cout << "found new optimality cut" << std::endl;
    IloNumArray u_0_values(worker_env);
    worker_cplex.getValues(u_0_values, dual_vars.u_0);
    IloNumArray u_1_values(worker_env);
    worker_cplex.getValues(u_1_values, dual_vars.u_1);
    IloNumArray u_2_values(worker_env);
    worker_cplex.getValues(u_2_values, dual_vars.u_2);
    IloNumArray u_3_values(worker_env);
    worker_cplex.getValues(u_3_values, dual_vars.u_3);

    const Graph * graph = instance.graph();
    const int num_vertices = graph->num_vertices();
    const int num_arcs = graph->num_arcs();
    const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
    const auto * vertices_info = graph->vertices_info();
    const int budget = instance.uncertainty_budget();
    const double route_limit = instance.limit();
    const int num_mandatory = instance.num_mandatory();

    for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
    {
      for(auto j: graph->AdjVerticesOut(0))
        cut_expr -= operator*(u_0_values[u_0_var_to_index(j,budget_iter,num_vertices)]*(*graph)[0][j]->distance(),(y[j]));

      for (int i = 0; i < num_vertices; ++i)
      {
        const auto vertex_info = vertices_info[i];
        double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
        if(i <= num_mandatory)
          cut_expr += u_1_values[u_1_var_to_index(i,budget_iter,num_vertices)]*route_limit;
        else
          cut_expr += operator*(u_1_values[u_1_var_to_index(i,budget_iter,num_vertices)]*vertex_deadline,y[i]);
      
        if (i > 0) // only profitable and mandatory.
        {
          for(auto j: graph->AdjVerticesOut(i))
          {
            GArc * arc = (*graph)[i][j];
            const int arc_pos = graph->pos(i,j);
            IloExpr coef1(master_env);
            coef1 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + arc->distance(),1-x[arc_pos]) - vertex_info.nominal_service_time_ - arc->distance());
            cut_expr += operator*(u_2_values[u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef1);
            if(budget_iter > 0)
            {
              IloExpr coef2(master_env);
              coef2 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance(),1-x[arc_pos]) - vertex_info.nominal_service_time_  - vertex_info.dev_service_time_ - arc->distance());
              cut_expr += operator*(u_3_values[u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef2);
            }
          }
        }
      }
    }

    cut_expr -= dual_bound;

    u_0_values.end();
    u_1_values.end();
    u_2_values.end();
    u_3_values.end();
    violatedCutFound = IloTrue;
  }

  return violatedCutFound;

} // END separate

int a_var_to_index(int vertex, int budget, int num_vertices)
{
    return budget * num_vertices + vertex;
}

int u_0_var_to_index(int vertex, int budget, int num_vertices)
{
  assert(num_vertices > 1);
  return budget * (num_vertices-1) + vertex - 1;
} 

int u_1_var_to_index(int vertex, int budget, int num_vertices)
{
  return budget * num_vertices + vertex;
}

int u_2_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs)
{
  assert(num_arcs > num_arcs_from_zero);
  return budget * (num_arcs-num_arcs_from_zero) + arc_pos - num_arcs_from_zero;
}

int u_3_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs)
{
  assert(num_arcs > num_arcs_from_zero);
  assert(budget > 0);
  return (budget-1) * (num_arcs-num_arcs_from_zero) + arc_pos - num_arcs_from_zero;
}

std::pair<int,int> index_to_a_var(int index, int num_vertices)
{
    int vertex = index%num_vertices;
    int budget = index/num_vertices;
    return std::pair<int,int>(vertex,budget);
}

static void AddNamesToDualVariables(DualVariables& dual_vars, Instance& inst)
{
  const Graph* graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = inst.uncertainty_budget();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));

  // add name to variables.
  // u_0.
  for(int i = 1; i < num_vertices; ++i)
  {
    for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
    {
      char strnum[28];
      sprintf(strnum,"u_0(%d)(%d)",i,budget_iter);
      dual_vars.u_0[u_0_var_to_index(i,budget_iter,num_vertices)].setName(strnum);
    }
  }

  // u_1.
  for(int i = 0; i < num_vertices; ++i)
  {
    for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
    {
      char strnum[28];
      sprintf(strnum,"u_1(%d)(%d)",i,budget_iter);
      dual_vars.u_1[u_1_var_to_index(i,budget_iter,num_vertices)].setName(strnum);
    }
  }

  // u_2.
  for(int i = 1; i < num_vertices; ++i)
  {
    for (int j: graph->AdjVerticesOut(i))
    {
      for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
      {
        char strnum[31];
        const int arc_pos = graph->pos(i,j);
        sprintf(strnum,"u_2(%d)(%d)(%d)",i,j,budget_iter);
        dual_vars.u_2[u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].setName(strnum);
      }
    }
  }

  // u_3.
  for(int i = 1; i < num_vertices; ++i)
  {
    for (int j: graph->AdjVerticesOut(i))
    {
      for (int budget_iter = 1; budget_iter <= budget; ++budget_iter)
      {
        char strnum[31];
        const int arc_pos = graph->pos(i,j);
        sprintf(strnum,"u_3(%d)(%d)(%d)",i,j,budget_iter);
        dual_vars.u_3[u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].setName(strnum);
      }
    }
  }
}

void SetPriorityOrder(IloCplex & cplex, IloEnv & env, IloNumVarArray & x, Instance & instance, double * R0)
{
  const Graph * graph = instance.graph();
  int v1 = 0, v2 = 0;

  for(std::list<int>::iterator it = (graph->AdjVerticesOut(v1)).begin(); it != (graph->AdjVerticesOut(v1)).end(); ++it)
  {
    v2 = *it;
    cplex.setPriority(x[graph->pos(v1,v2)],1.0);
  }

  /*IloNumArray pri_order(env, graph->num_arcs());

  for(v1 = 0; v1 < num_vertices; v1++)
  {
  // highest priorities for arcs leaving source
  if(v1 == 0) curr_pri = 1.0;
  else curr_pri = 0.0;

  for(std::list<int>::iterator it = ((*adj_lists)[v1]).begin(); it != ((*adj_lists)[v1]).end(); it++)
  {
  v2 = *it;
  pri_order[graph->pos(v1,v2)] = curr_pri + 1.0/(R0[v1] + ((*graph)[v1][v2])->dist());
}
}

cplex.setPriorities(x, pri_order);*/
}

Solution<int> * optimize(IloCplex & cplex, IloEnv& env, IloModel& model, IloNumVarArray & y, IloNumVarArray & x, std::optional<std::reference_wrapper<IloNumVar>> dual_bound_opt, Instance& instance, double total_time_limit, bool solve_relax, bool callback, bool find_root_cuts, double * R0, double * Rn, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol)
{
  const Graph* graph = instance.graph();
  int num_vertices = graph->num_vertices();
  Timestamp * ti = NewTimestamp(), *tf = NewTimestamp();
  Timer * timer = GetTimer();
  Solution<int> * solution = new Solution<int>(num_vertices);
  std::optional<IloEnv> worker_env = std::nullopt;
  std::optional<IloCplex> worker_cplex = std::nullopt;
  std::optional<IloObjective> worker_obj = std::nullopt;
  std::optional<DualVariables> worker_vars = std::nullopt;

  cplex.extract(model);
  cplex.setParam(IloCplex::Param::WorkMem,15000);
  cplex.setParam(IloCplex::IloCplex::Param::MIP::Strategy::File,3);

  if(callback)
  {
    // create subproblem!
    // The subproblem will be always the same, except for the objective function, which relies
    // on the values of X and y variables. Then, we can create a unique model beforehand.
    worker_env = IloEnv();
    worker_cplex = IloCplex(*worker_env);
    worker_cplex->setOut((*worker_env).getNullStream());
    worker_vars = DualVariables{};

    // Turn off the presolve reductions and set the CPLEX optimizer
    // to solve the worker LP with primal simplex method.
    worker_cplex->setParam(IloCplex::Param::Preprocessing::Reduce, 0);
    worker_cplex->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

    IloModel worker_model(*worker_env);
    worker_cplex->extract(worker_model);

    // Set empty initial objective function.
    worker_obj = IloObjective(*worker_env);
    worker_obj->setSense(IloObjective::Minimize);
    worker_model.add(*worker_obj);

    worker_vars = DualVariables{};
    AllocateDualVariables(*worker_env,*worker_vars,instance);

    PopulateByRowDualCompactBaselineContinuousSpace(*worker_env,worker_model,*worker_vars,instance);
    SetDualVariablesProperties(env,*worker_vars,x,y,instance);

    cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse); 
    // cplex.setParam(IloCplex::Param::Threads, 1); 
    // Turn on traditional search for use with control callbacks.
    cplex.setParam(IloCplex::Param::MIP::Strategy::Search,
                          IloCplex::Traditional);
    cplex.use(BendersLazyCallback(env,*worker_cplex,*worker_obj,*worker_vars,y,x,dual_bound_opt->get(),instance));
    cplex.use(BendersUserCutCallback(env,*worker_cplex,*worker_obj,*worker_vars,y,x,dual_bound_opt->get(),instance));
  }

  if(!K_MULTI_THREAD) cplex.setParam(IloCplex::Param::Threads, 1);

  if(solve_relax) total_time_limit = -1;

  if(!double_equals(total_time_limit,-1))
  {
    //std::cout << total_time_limit - instance.time_spent_in_preprocessing_ << std::endl;
    cplex.setParam(IloCplex::Param::ClockType,2);
    cplex.setParam(IloCplex::Param::TimeLimit, total_time_limit - instance.time_spent_in_preprocessing());
    //std::cout << total_time_limit << " - " << instance.time_spent_in_preprocessing_ << std::endl;
  }

  timer->Clock(ti);
  double curr_bound = std::numeric_limits<double>::infinity();

  if (!cplex.solve())
  {
    timer->Clock(tf);
    if(callback)
    {
      // delete subproblem objects.
      (*worker_cplex).end();
      (*worker_obj).end();
      (*worker_env).end();
    }
    if(solve_relax) solution->root_time_ = timer->ElapsedTime(ti,tf);
    else solution->milp_time_ += timer->ElapsedTime(ti,tf);
    if((cplex.getCplexStatus() == IloCplex::Infeasible)||(cplex.getCplexStatus() == IloCplex::InfOrUnbd)) solution->is_feasible_ = false;
    return solution;
  }

  curr_bound = cplex.getObjValue();

  timer->Clock(tf);
  if(solve_relax) solution->root_time_ = timer->ElapsedTime(ti,tf);
  else solution->milp_time_ += (timer->ElapsedTime(ti,tf) + instance.time_spent_in_preprocessing());

  SetSolutionStatus(cplex,*solution,solve_relax);

  delete(ti);
  ti = nullptr;
  delete(tf);
  tf = nullptr;

  // delete subproblem objects.
  if(callback)
  {
    (*worker_cplex).end();
    (*worker_obj).end();
    (*worker_env).end();
  }

  return solution;
}

Solution<double> * optimizeLP(IloCplex & cplex, IloEnv& env, IloModel& model, Instance& instance, double total_time_limit, double * R0, double * Rn)
{
  const Graph* graph = instance.graph();
  int num_vertices = graph->num_vertices();
  Timestamp * ti = NewTimestamp(), *tf = NewTimestamp();
  Timer * timer = GetTimer();
  Solution<double> * solution = new Solution<double>(num_vertices);

  cplex.extract(model);
  cplex.setParam(IloCplex::Param::WorkMem,15000);
  cplex.setParam(IloCplex::IloCplex::Param::MIP::Strategy::File,3);

  if(!K_MULTI_THREAD) cplex.setParam(IloCplex::Param::Threads, 1);

  if(!double_equals(total_time_limit,-1))
  {
    //std::cout << total_time_limit - instance.time_spent_in_preprocessing_ << std::endl;
    cplex.setParam(IloCplex::Param::ClockType,2);
    cplex.setParam(IloCplex::Param::TimeLimit, total_time_limit - instance.time_spent_in_preprocessing());
    //std::cout << total_time_limit << " - " << instance.time_spent_in_preprocessing_ << std::endl;
  }

  timer->Clock(ti);
  double curr_bound = std::numeric_limits<double>::infinity();

  if (!cplex.solve())
  {
    timer->Clock(tf);
    timer->ElapsedTime(ti,tf);
    if((cplex.getCplexStatus() == IloCplex::Infeasible)||(cplex.getCplexStatus() == IloCplex::InfOrUnbd)) solution->is_feasible_ = false;
    return solution;
  }

  curr_bound = cplex.getObjValue();

  timer->Clock(tf);
  solution->root_time_ = timer->ElapsedTime(ti,tf);

  SetSolutionStatus(cplex,*solution,false);

  delete(ti);
  ti = nullptr;
  delete(tf);
  tf = nullptr;

  return solution;
}

Solution<int> * BendersCompactBaseline(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles, bool export_model)
{
  const Graph * graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_routes = inst.num_vehicles();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));

  IloEnv env;
  IloCplex cplex(env);
  //cplex.setOut(env.getNullStream());
  IloModel model(env);

  IloNumVar slack(env, 0, force_use_all_vehicles? 0: num_routes, ILOFLOAT);
  IloNumVar dual_bound(env, -IloInfinity, 0, ILOFLOAT);
  IloNumVarArray x(env);
  IloNumVarArray y(env);

  if(solve_relax)
  {
        x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOFLOAT);
        y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOFLOAT);
  }
  else
  {
        x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOINT);
        y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOINT);
  }

  PopulateByRowCommon(env,model,slack,y,x,inst,force_use_all_vehicles);
	
  const int num_mandatory = inst.num_mandatory();
  const auto vertices_info = graph->vertices_info();

  // add objective function.
  IloExpr obj(env);

  for(int i = num_mandatory + 1; i < num_vertices; ++i)
  {
    const auto& vertex_info = vertices_info[i];
    obj += operator*(vertex_info.profit_,y[i]);
  }

  obj += dual_bound;
  model.add(IloMaximize(env, obj));
  obj.end();

  if(export_model)
  {
    // add name to variables.
    for(int i = 0; i < num_vertices; ++i)
    {
      char strnum3[15];
      sprintf(strnum3,"y(%d)",i);
      y[i].setName(strnum3);

      for(const auto j: graph->AdjVerticesOut(i))
      {
        char strnum[26];
        sprintf(strnum,"x(%d)(%d)",i,j);
        x[graph->pos(i,j)].setName(strnum);
      }
    }

    slack.setName("slack");
    dual_bound.setName("dual_bound");

    cplex.extract(model);
    cplex.exportModel("model_Benders_compact_baseline.lp");
  }

  auto result = optimize(cplex,env,model,y,x,dual_bound,inst,time_limit,solve_relax,callback,find_root_cuts, R0, Rn, initial_cuts, initial_sol);

  // // print solution.
  // IloNumArray y_solution(env);
  // cplex.getValues(y_solution,y);
  // for (int i = 0; i < num_vertices; ++i)
  //   std::cout << "y[" << i << "]: " << y_solution[i] << std::endl;

  // IloNumArray x_solution(env);
  // cplex.getValues(x_solution,x);
  // for (int i = 0; i < num_vertices; ++i)
  //   for(const int &j: graph->AdjVerticesOut(i))
  //     std::cout << "x[" << i << "," << j << "]: " << x_solution[graph->pos(i,j)] << std::endl;

  // std::cout << "slack: " << cplex.getValue(slack) << std::endl;
  // std::cout << "dual bound: " << cplex.getValue(dual_bound) << std::endl;
  // x_solution.end();
  // y_solution.end();
  
  cplex.end();
  env.end();
  return result;
}

Solution<int> * CompactBaseline(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles, bool export_model)
{
  const Graph * graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_routes = inst.num_vehicles();
  const int budget = inst.uncertainty_budget();

  IloEnv env;
  IloCplex cplex(env);
  cplex.setOut(env.getNullStream());
  IloModel model(env);

  IloNumVar slack(env, 0, force_use_all_vehicles? 0: num_routes, ILOFLOAT);
  IloNumVarArray x(env);
  IloNumVarArray y(env);
  // budget + 1 to consider level 0 of budget! 0,..., budget
  IloNumVarArray a(env, (budget+1)*num_vertices, 0, IloInfinity, ILOFLOAT);

  if(solve_relax)
  {
        x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOFLOAT);
        y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOFLOAT);
    }
  else
  {
        x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOINT);
        y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOINT);
    }

	PopulateByRowCompactBaseline(cplex,env,model,slack,y,x,a,inst,R0,Rn,force_use_all_vehicles,export_model);

  auto result = optimize(cplex,env,model,y,x,std::nullopt,inst,time_limit,solve_relax,callback,find_root_cuts, R0, Rn, initial_cuts, initial_sol);

  // // print solution.
  // IloNumArray y_solution(env);
  // cplex.getValues(y_solution,y);
  // for (int i = 0; i < num_vertices; ++i)
  //   std::cout << "y[" << i << "]: " << y_solution[i] << std::endl;

  // IloNumArray x_solution(env);
  // cplex.getValues(x_solution,x);
  // for (int i = 0; i < num_vertices; ++i)
  //   for(const int &j: graph->AdjVerticesOut(i))
  //     std::cout << "x[" << i << "," << j << "]: " << x_solution[graph->pos(i,j)] << std::endl;

  // IloNumArray a_solution(env);
  // cplex.getValues(a_solution,a);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;

  // x_solution.end();
  // y_solution.end();
  // a_solution.end();

  cplex.end();
  env.end();
  return result;
}

Solution<double> * PrimalSubproblemCompactBaseline(Instance& inst, IloNumArray& x_values, IloNumArray& y_values, double * R0, double * Rn, double time_limit, bool export_model)
{
  const Graph * graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = inst.uncertainty_budget();

  IloEnv env;
  IloCplex cplex(env);
  cplex.setOut(env.getNullStream());
  IloModel model(env);

  // budget + 1 to consider level 0 of budget! 0,..., budget
  IloNumVarArray a(env, (budget+1)*num_vertices, 0, IloInfinity, ILOFLOAT);

  // empty.
  IloNumVarArray y(env);
  IloNumVarArray x(env);

  PopulateByRowCompactBaselineContinuousSpace(env,model,x,y,a,x_values,y_values,inst);

   // add objective function.
  IloExpr obj(env);

  const int num_mandatory = inst.num_mandatory();
  const auto* vertices_info = inst.graph()->vertices_info();

  for(int i = num_mandatory + 1; i < num_vertices; ++i)
  {
    const auto& vertex_info = vertices_info[i];
    obj -= operator*(vertex_info.decay_ratio_,a[a_var_to_index(i,budget,num_vertices)]);
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  if(export_model)
  {
    // add name to variables.
    for(int i = 0; i < num_vertices; ++i)
    {
      for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
      {
        char strnum[26];
        sprintf(strnum,"a(%d)(%d)",i,budget_iter);
        a[a_var_to_index(i,budget_iter,num_vertices)].setName(strnum);
      }
    }

    cplex.extract(model);
    cplex.exportModel("primal_continuous_model_compact_baseline.lp");
  }

  auto result = optimizeLP(cplex,env,model,inst,time_limit, R0, Rn);

  // // print solution.
  // IloNumArray a_solution(env);
  // cplex.getValues(a_solution,a);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;
  // a_solution.end();
  cplex.end();
  env.end();
  return result;
}

static void SetDualVariablesProperties(IloEnv& master_env, DualVariables& dual_vars, IloNumVarArray& x, IloNumVarArray& y, Instance& instance)
{
  const Graph * graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto * vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();
  const double route_limit = instance.limit();
  const int num_mandatory = instance.num_mandatory();

  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for(auto j: graph->AdjVerticesOut(0))
    {
      IloExpr exp1(master_env);
      exp1 -= operator*((*graph)[0][j]->distance(),y[j]);
      dual_vars.ext_ray_coef[dual_vars.u_0[u_0_var_to_index(j,budget_iter,num_vertices)].getId()] = exp1;
    }

    for (int i = 0; i < num_vertices; ++i)
    {
      const auto vertex_info = vertices_info[i];
      double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
      
      IloExpr exp2(master_env);
      if(i <= num_mandatory)
      {
        exp2 += route_limit;
        dual_vars.ext_ray_coef[dual_vars.u_1[u_1_var_to_index(i,budget_iter,num_vertices)].getId()] = exp2;
      }else
      {
        exp2 +=operator*(vertex_deadline,y[i]);
        dual_vars.ext_ray_coef[dual_vars.u_1[u_1_var_to_index(i,budget_iter,num_vertices)].getId()] = exp2;
      }
      if (i > 0) // only profitable and mandatory.
      {
        for(auto j: graph->AdjVerticesOut(i))
        {
          GArc * arc = (*graph)[i][j];
          const int arc_pos = graph->pos(i,j);
          IloExpr exp3(master_env);
          exp3 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + arc->distance(),1-x[arc_pos]) - vertex_info.nominal_service_time_ - arc->distance());
          dual_vars.ext_ray_coef[dual_vars.u_2[u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].getId()] = exp3;

          if(budget_iter > 0)
          {
            IloExpr exp4(master_env);
            exp4 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance(),1-x[arc_pos]) - vertex_info.nominal_service_time_  - vertex_info.dev_service_time_ - arc->distance());
            dual_vars.ext_ray_coef[dual_vars.u_3[u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].getId()] = exp4;
          }
        }
      }
    }
  }
}

static void FillObjectiveExpressionDualCompactBaselineContinuousSpace(IloExpr& obj, DualVariables & dual_vars, IloNumArray& x_values, IloNumArray& y_values, Instance& instance)
{
  const Graph * graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto * vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();
  const double route_limit = instance.limit();
  const int num_mandatory = instance.num_mandatory();

  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for(auto j: graph->AdjVerticesOut(0))
      obj -= operator*(dual_vars.u_0[u_0_var_to_index(j,budget_iter,num_vertices)], (*graph)[0][j]->distance())*(y_values[j]);

    for (int i = 0; i < num_vertices; ++i)
    {
      const auto vertex_info = vertices_info[i];
      double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
      
      if(i <= num_mandatory)
        obj += operator*(dual_vars.u_1[u_1_var_to_index(i,budget_iter,num_vertices)], route_limit);
      else
        obj += operator*(dual_vars.u_1[u_1_var_to_index(i,budget_iter,num_vertices)], vertex_deadline * y_values[i]);
    
      if (i > 0) // only profitable and mandatory.
      {
        for(auto j: graph->AdjVerticesOut(i))
        {
          GArc * arc = (*graph)[i][j];
          const int arc_pos = graph->pos(i,j);
          double coef1 = (vertex_deadline + vertex_info.nominal_service_time_ + arc->distance())*(1-x_values[arc_pos]) - vertex_info.nominal_service_time_ - arc->distance();
          obj += operator*(dual_vars.u_2[u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef1);

          if(budget_iter > 0)
          {
            double coef2 = (vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance())*(1-x_values[arc_pos]) - vertex_info.nominal_service_time_  - vertex_info.dev_service_time_ - arc->distance();
            obj += operator*(dual_vars.u_3[u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef2);
          }
        }
      }
    }
  }
}

static void AllocateDualVariables(IloEnv& env, DualVariables& dual_vars, Instance & instance)
{
  const Graph * graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = instance.uncertainty_budget();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));

  // budget + 1 to consider level 0 of budget! 0,..., budget
  // u_0 defined for all i \in S U P, 0...budget
  dual_vars.u_0 = IloNumVarArray(env, (budget+1)*(num_vertices-1), 0, IloInfinity, ILOFLOAT);
  // u_1 defined for all i \in N, 0...budget
  dual_vars.u_1 = IloNumVarArray(env, (budget+1)*num_vertices, 0, IloInfinity, ILOFLOAT);
  // u_2 defined for all arcs (i,j) i != 0, 0...budget
  dual_vars.u_2 = IloNumVarArray(env, (budget+1)*(num_arcs-num_arcs_from_origin), 0, IloInfinity, ILOFLOAT);
  // u_3 defined for all arcs (i,j) i != 0, 1...budget (considering 0...budget -1)
  dual_vars.u_3 = IloNumVarArray(env, budget*(num_arcs-num_arcs_from_origin), 0, IloInfinity, ILOFLOAT);
}

Solution<double> * DualSubproblemCompactBaseline(Instance& inst, IloNumArray& x_values, IloNumArray& y_values, double * R0, double * Rn, double time_limit, bool export_model)
{
  const Graph * graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();

  IloEnv env;
  IloCplex cplex(env);
  cplex.setOut(env.getNullStream());
  IloModel model(env);

  DualVariables dual_vars;
  AllocateDualVariables(env,dual_vars,inst);

  PopulateByRowDualCompactBaselineContinuousSpace(env,model,dual_vars,inst);

  // add objective function.
  IloExpr obj_expr(env);
  FillObjectiveExpressionDualCompactBaselineContinuousSpace(obj_expr,dual_vars,x_values,y_values,inst);
  model.add(IloMinimize(env, obj_expr));
  obj_expr.end();

  if(export_model)
  {
    AddNamesToDualVariables(dual_vars,inst);
    cplex.extract(model);
    cplex.exportModel("dual_continuous_model_compact_baseline.lp");
  }

  auto result = optimizeLP(cplex,env,model,inst,time_limit, R0, Rn);

  // // print solution.
  // IloNumArray a_solution(env);
  // cplex.getValues(a_solution,a);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;

  cplex.end();
  env.end();
  return result;
}

Solution<int> * CompactSingleCommodity(Instance& inst, double * R0, double * Rn, double time_limit, bool solve_relax, bool callback, bool find_root_cuts, std::list<UserCutGeneral*> * initial_cuts, HeuristicSolution * initial_sol, bool force_use_all_vehicles)
{
  const Graph * graph = inst.graph();
  int num_vertices = graph->num_vertices();
  int num_arcs = graph->num_arcs();
  int num_routes = inst.num_vehicles();

  IloEnv env;
  IloCplex cplex(env);
  cplex.setOut(env.getNullStream());
  IloModel model(env);

  IloNumVar slack(env, 0, num_routes, ILOFLOAT);
  IloNumVarArray x(env);
  IloNumVarArray y(env);
  IloNumVarArray f(env, num_arcs, 0, IloInfinity, ILOFLOAT);

  if(solve_relax)
  {
        x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOFLOAT);
        y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOFLOAT);
    }
  else
  {
        x = IloNumVarArray(env, num_arcs, 0.0, 1.0, ILOINT);
        y = IloNumVarArray(env, num_vertices, 0.0, 1.0, ILOINT);
    }

	PopulateByRowCompactSingleCommodity(env,model,slack,y,x,f,inst,R0,Rn,force_use_all_vehicles);

  auto result = optimize(cplex,env,model,y,x,std::nullopt,inst,time_limit,solve_relax,callback,find_root_cuts, R0, Rn, initial_cuts, initial_sol);

  cplex.end();
  env.end();
  return result;
}

static void PopulateByRowCompactSingleCommodity(IloEnv& env, IloModel& model, IloNumVar& slack, IloNumVarArray & y, IloNumVarArray & x, IloNumVarArray & f, Instance& instance, double * R0, double* Rn, bool force_use_all_vehicles)
{
  int num_vehicles = instance.num_vehicles();
  int num_vertices = (instance.graph())->num_vertices();
  int num_mandatory = instance.num_mandatory();
  GArc * curr_arc = nullptr, *curr_inv_arc = nullptr;
  const Graph* graph = instance.graph();

  PopulateByRowCommon(env,model,slack,y,x,instance,force_use_all_vehicles);
 
/*IloExpr expv(env);
for(int i = 0; i < num_vertices; i++)
{
if((*(instance.graph()))[i][num_vertices-1] != nullptr) expv +=f[graph->pos(i,num_vertices-1)];
if((*(instance.graph()))[i][0] != nullptr) expv +=f[graph->pos(i,0)];
}
model.add(expv == 0);
expv.end();*/

/*for(int i = 0; i < num_vertices; i++)
{
if((*(instance.graph()))[i][num_vertices-1] != nullptr)  f[graph->pos(i,num_vertices-1)].setUB(0.0);
}*/

for(int i = 1; i < num_vertices; i++)
{
  curr_arc = (*(instance.graph()))[0][i];

  if(curr_arc != nullptr)
  {
    IloExpr exp(env);
    exp += operator*(instance.limit() - curr_arc->distance(),x[graph->pos(0,i)]);
    model.add(f[graph->pos(0,i)] == exp);
    exp.end();
  }
}

/*IloExpr exp2(env);
IloExpr exp1(env);
for(int i = 1; i < num_vertices-1; i++)
{
for(int j = 0; j < num_vertices; j++)
{
curr_arc = (*(instance.graph()))[i][j];
if(curr_arc != nullptr)
{
exp2 += operator*(curr_arc->dist(),x[graph->pos(i,j)]);
}
}

if((*(instance.graph()))[0][i] != nullptr) exp1 += f[graph->pos(0,i)];
}

if((*(instance.graph()))[0][num_vertices-1] != nullptr) exp1 += f[graph->pos(0,num_vertices-1)];

model.add(exp1 == exp2);

exp1.end();
exp2.end();*/

for(int i = 1; i < num_vertices-1; i++)
{
  IloExpr exp1(env);
  IloExpr exp2(env);

  for(int j = 0; j < num_vertices; j++)
  {
    if((*(instance.graph()))[j][i] != nullptr)
    {
      exp1 += f[graph->pos(j,i)];
    }

    curr_arc = (*(instance.graph()))[i][j];
    if(curr_arc != nullptr)
    {
      exp1 -= f[graph->pos(i,j)];
      exp2 += operator*(curr_arc->distance(), x[graph->pos(i,j)]);
    }
  }

  //if(i <= num_mandatory) model.add(exp == 1);
  model.add(exp1 == exp2);

  exp1.end();
  exp2.end();
}

for(int i = 0; i < num_vertices; i++)
{
  for(int j = 0; j < num_vertices; j++)
  {
    curr_arc = (*(instance.graph()))[i][j];
    if(curr_arc != nullptr)
    {
      IloExpr exp(env);
      if(R0 != nullptr) exp += operator*(instance.limit() - R0[i] - curr_arc->distance(),x[graph->pos(i,j)]);
      else  exp += operator*(instance.limit() - curr_arc->distance(),x[graph->pos(i,j)]);
      model.add(f[graph->pos(i,j)] <= exp);
      exp.end();
    }
  }
}

IloExpr obj(env);

for(int j = num_mandatory + 1; j < num_vertices-1; j++)
{
  /*IloExpr exp(env);
  for(int i = 0; i < num_vertices; i++)
  {
    if((*(instance.graph()))[i][j] != nullptr)
    {
      exp += x[graph->pos(i,j)];
    }
  }*/

  //obj += operator*(((instance.graph())->profits())[j],exp);
  obj += operator*((instance.graph()->vertices_info())[j].profit_,y[j]);
  //exp.end();
}

model.add(IloMaximize(env, obj));
obj.end();

for(int i = 0; i < num_vertices; i++)
{
  char strnum3[15];
  sprintf(strnum3,"y(%d)",i);
  y[i].setName(strnum3);
  for(int j = 0; j < num_vertices; j++)
  {
    curr_arc = (*(instance.graph()))[i][j];
    if(curr_arc != nullptr)
    {
      char strnum[26], strnum2[26];
      sprintf(strnum,"f(%d)(%d)",i,j);
      f[graph->pos(i,j)].setName(strnum);
      sprintf(strnum2,"x(%d)(%d)",i,j);
      x[graph->pos(i,j)].setName(strnum2);
    }
  }
}

slack.setName("slack");
}

static void PopulateByRowCommon(IloEnv& env, IloModel& model, IloNumVar& slack, IloNumVarArray & y, IloNumVarArray & x, Instance& instance, bool force_use_all_vehicles)
{
  int num_vehicles = instance.num_vehicles();
  int num_vertices = (instance.graph())->num_vertices();
  int num_mandatory = instance.num_mandatory();
  GArc * curr_arc = nullptr, *curr_inv_arc = nullptr;
  const Graph* graph = instance.graph();

  for(int i = 1; i < num_vertices; ++i)
  {
    IloExpr exp(env);
    auto adj_vertices = graph->AdjVerticesOut(i);
    size_t num_adj_arc = adj_vertices.size();

    for (const auto j: adj_vertices)
      exp += x[graph->pos(i,j)];

    if (i <= num_mandatory)
    {
      if(num_adj_arc == 0) model.add(y[i] == 0); // creates infeasibility if a mandatory node is disconnected from graph!
      else model.add(exp == 1);
      model.add(y[i] == 1);
    }
    else
    {
      if(num_adj_arc == 0) y[i].setUB(0);
      else model.add(exp == y[i]);
    }
    exp.end();
  }

  model.add(y[0] == 1);

  IloExpr expi(env);
  for(const auto& j: graph->AdjVerticesOut(0))
    expi += x[graph->pos(0,j)];

  expi += slack;

  model.add(expi == num_vehicles);
  expi.end();

  if(force_use_all_vehicles)
    slack.setUB(0);

  for(int i = 0; i < num_vertices; ++i)
  {
    IloExpr exp(env);
    for(int j = 0; j < num_vertices; ++j)
    {
      if((*(instance.graph()))[i][j]) exp += x[graph->pos(i,j)];
      if((*(instance.graph()))[j][i]) exp -= x[graph->pos(j,i)];
    }
    model.add(exp == 0);
    exp.end();
  }
}

static void PopulateByRowDualCompactBaselineContinuousSpace(IloEnv& env, IloModel& model, DualVariables &dual_vars, Instance& instance)
{
  const Graph* graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_mandatory = instance.num_mandatory();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto* vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();

  IloExpr exp(env);

  exp += dual_vars.u_1[u_1_var_to_index(0,0,num_vertices)];
  for (auto j: graph->AdjVerticesIn(0))
    exp -= dual_vars.u_2[u_2_var_to_index(graph->pos(j,0),0, num_arcs_from_origin, num_arcs)];

  model.add(exp >= 0);
  exp.end();

  for(int budget_iter = 1; budget_iter <= budget; ++budget_iter)
  {
    IloExpr exp(env);
    exp += dual_vars.u_1[u_1_var_to_index(0,budget_iter,num_vertices)];

    for (auto j: graph->AdjVerticesIn(0))
    {
      exp -= dual_vars.u_2[u_2_var_to_index(graph->pos(j,0),budget_iter, num_arcs_from_origin, num_arcs)];
      exp -= dual_vars.u_3[u_3_var_to_index(graph->pos(j,0),budget_iter, num_arcs_from_origin, num_arcs)];
    }
    model.add(exp >= 0);
    exp.end();
  }

  for(int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for (int i = 1; i < num_vertices; ++i)
    {
      IloExpr exp(env);

      exp -= dual_vars.u_0[u_0_var_to_index(i,budget_iter,num_vertices)];
      exp += dual_vars.u_1[u_1_var_to_index(i,budget_iter,num_vertices)];

      for(auto j: graph->AdjVerticesOut(i))
      {
        exp += dual_vars.u_2[u_2_var_to_index(graph->pos(i,j),budget_iter,num_arcs_from_origin,num_arcs)];
        if(budget_iter < budget)
          exp += dual_vars.u_3[u_3_var_to_index(graph->pos(i,j),budget_iter+1,num_arcs_from_origin,num_arcs)];
      }

      for(auto j: graph->AdjVerticesIn(i))
      {
        if (j != 0)
        {
          exp -= dual_vars.u_2[u_2_var_to_index(graph->pos(j,i),budget_iter,num_arcs_from_origin,num_arcs)];
          if (budget_iter > 0)
            exp -= dual_vars.u_3[u_3_var_to_index(graph->pos(j,i),budget_iter,num_arcs_from_origin,num_arcs)];
        }
      }

      const auto coef = (budget_iter == budget && i > num_mandatory)? - vertices_info[i].decay_ratio_: 0.0;
      model.add(exp >= coef);
      exp.end();


      if(!(*graph)[0][i])
        dual_vars.u_0[u_0_var_to_index(i,budget_iter,num_vertices)].setUB(0.0);
    }
  }
}

static void PopulateByRowCompactBaselineContinuousSpace(IloEnv& env, IloModel& model, IloNumVarArray & x, IloNumVarArray & y, IloNumVarArray &a, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, Instance& instance)
{
  const int num_vertices = (instance.graph())->num_vertices();
  const int num_mandatory = instance.num_mandatory();
  GArc * curr_arc = nullptr;
  const Graph* graph = instance.graph();
  const auto* vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();
  const double route_time_limit = instance.limit();

  for(int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for (const auto &vertex : graph->AdjVerticesOut(0))
    {
      curr_arc = *graph[0][vertex];
      assert(curr_arc);
      y_values.has_value()?
        model.add(a[a_var_to_index(vertex,budget_iter,num_vertices)] >= curr_arc->distance()*(y_values->get())[vertex])
        : model.add(a[a_var_to_index(vertex,budget_iter,num_vertices)] >= operator*(curr_arc->distance(), y[vertex]));
    }

    // i is zero or mandatory.
    for (int i = 0; i <= num_mandatory; ++i)
      model.add(a[a_var_to_index(i,budget_iter,num_vertices)] <= route_time_limit);

    // i is profitable.
    for (int i = num_mandatory+1; i < num_vertices; ++i)
    {
      const auto & vertex_info = vertices_info[i];
      double vertex_deadline = round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
      y_values.has_value()?
        model.add(a[a_var_to_index(i,budget_iter,num_vertices)] <= vertex_deadline * (y_values->get())[i])
        : model.add(a[a_var_to_index(i,budget_iter,num_vertices)] <= operator*(vertex_deadline,y[i]));
    }

    // robust constraints.
    for (int i = 1; i < num_vertices; ++i)
    {
      const auto& vertex_info = vertices_info[i];
      // if mandatory, the deadline of the vertex D_i is the time limit T.
      const double vertex_deadline = (i <= num_mandatory)? route_time_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
        
      for (const int& j: graph->AdjVerticesOut(i))
      {
        IloExpr exp(env);
        curr_arc = (*graph)[i][j];
        assert(curr_arc);
        exp += (a[a_var_to_index(i,budget_iter,num_vertices)] + vertex_info.nominal_service_time_
            + curr_arc->distance() - a[a_var_to_index(j,budget_iter,num_vertices)]);
        
        x_values.has_value()?
          exp -= (vertex_deadline + vertex_info.nominal_service_time_ + curr_arc->distance()) * (1 - ((x_values->get())[graph->pos(i,j)]))
          : exp -= operator*(vertex_deadline + vertex_info.nominal_service_time_ + curr_arc->distance(),1-x[graph->pos(i,j)]);
        
        model.add(exp <= 0);
        exp.end();

        if(budget_iter > 0)
        {
          IloExpr exp(env);
          curr_arc = (*graph)[i][j];
          assert(curr_arc);
          exp += (a[a_var_to_index(i,budget_iter-1,num_vertices)] + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_
            + curr_arc->distance() - a[a_var_to_index(j,budget_iter,num_vertices)]);
          
          x_values.has_value()?
            exp -= (vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + curr_arc->distance()) * (1-(x_values->get())[graph->pos(i,j)])
            : exp -= operator*(vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + curr_arc->distance(),1-x[graph->pos(i,j)]);
          
          model.add(exp <= 0);
          exp.end();
        }
      }
    }
  }
}

static void PopulateByRowCompactBaseline(IloCplex& cplex, IloEnv& env, IloModel& model, IloNumVar& slack, IloNumVarArray & y, IloNumVarArray & x, IloNumVarArray & a, Instance& instance, double * R0, double* Rn, bool force_use_all_vehicles, bool export_model)
{
  int num_vertices = (instance.graph())->num_vertices();
  int num_mandatory = instance.num_mandatory();
  GArc * curr_arc = nullptr;
  const Graph* graph = instance.graph();
  const auto* vertices_info = graph->vertices_info();

  PopulateByRowCommon(env,model,slack,y,x,instance,force_use_all_vehicles);

  PopulateByRowCompactBaselineContinuousSpace(env,model,x,y,a,std::nullopt,std::nullopt,instance);

  // add objective function.
  IloExpr obj(env);
  const int budget = instance.uncertainty_budget();

  for(int i = num_mandatory + 1; i < num_vertices; ++i)
  {
    const auto& vertex_info = vertices_info[i];
    obj += operator*(vertex_info.profit_,y[i]);
    obj -= operator*(vertex_info.decay_ratio_,a[a_var_to_index(i,budget,num_vertices)]);
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  if(export_model)
  {
    // add name to variables.
    for(int i = 0; i < num_vertices; ++i)
    {
      char strnum3[15];
      sprintf(strnum3,"y(%d)",i);
      y[i].setName(strnum3);

      for(const auto j: graph->AdjVerticesOut(i))
      {
        char strnum[26];
        sprintf(strnum,"x(%d)(%d)",i,j);
        x[graph->pos(i,j)].setName(strnum);
      }

      for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
      {
        char strnum[26];
        sprintf(strnum,"a(%d)(%d)",i,budget_iter);
        a[a_var_to_index(i,budget_iter,num_vertices)].setName(strnum);
      }
    }

    slack.setName("slack");

    cplex.extract(model);
    cplex.exportModel("model_compact_baseline.lp");
  }
}