#include "src/benders_generic_callback.h"
#include "src/solution.hpp"

typedef lemon::ListDigraph LemonGraph;
typedef int LimitValueType;

ILOSTLBEGIN

void FillObjectiveExpressionDualCompactBaselineContinuousSpace(IloExpr& obj, DualVariablesBaseline& dual_vars, IloNumArray& x_values, IloNumArray& y_values, std::optional<IloNum> dual_bound_value, const Instance& instance, bool combine_feas_op_cuts)
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
      obj -= operator*(dual_vars.u_0_[dual_vars.u_0_var_to_index(j,budget_iter,num_vertices)], (*graph)[0][j]->distance())*(y_values[j]);

    for (int i = 0; i < num_vertices; ++i)
    {
      const auto vertex_info = vertices_info[i];
      double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
      
      if(i <= num_mandatory)
        obj += operator*(dual_vars.u_1_[dual_vars.u_1_var_to_index(i,budget_iter,num_vertices)], route_limit);
      else
        obj += operator*(dual_vars.u_1_[dual_vars.u_1_var_to_index(i,budget_iter,num_vertices)], vertex_deadline * y_values[i]);
    
      if (i > 0) // only profitable and mandatory.
      {
        for(auto j: graph->AdjVerticesOut(i))
        {
          GArc * arc = (*graph)[i][j];
          const int arc_pos = graph->pos(i,j);
          double coef1 = (vertex_deadline + vertex_info.nominal_service_time_ + arc->distance())*(1-x_values[arc_pos]) - vertex_info.nominal_service_time_ - arc->distance();
          obj += operator*(dual_vars.u_2_[dual_vars.u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef1);

          if(budget_iter > 0)
          {
            double coef2 = (vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance())*(1-x_values[arc_pos]) - vertex_info.nominal_service_time_  - vertex_info.dev_service_time_ - arc->distance();
            obj += operator*(dual_vars.u_3_[dual_vars.u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef2);
          }
        }
      }
    }
  }

  if(combine_feas_op_cuts)
    obj -= operator*(*dual_bound_value,dual_vars.u_dual_bound_);
}

void FillObjectiveExpressionDualCompactSingleCommodityContinuousSpace(IloExpr& obj, DualVariablesSingleCommodity& dual_vars, IloNumArray& x_values, IloNumArray& y_values, std::optional<IloNum> dual_bound_value, const double* R0, const double* Rn, const Instance& instance, bool combine_feas_op_cuts)
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
    for (int i = 0; i < num_vertices; ++i)
    {
      const auto vertex_info = vertices_info[i];
      double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);

      for (int j : graph->AdjVerticesOut(i))
      {
        const GArc * arc = (*graph)[i][j];
        const int arc_pos = graph->pos(i,j);
        double coef1 = (route_limit - R0[i] - vertex_info.nominal_service_time_ - arc->distance())*x_values[arc_pos];
        obj += operator*(dual_vars.v_0_[dual_vars.v_0_var_to_index(arc_pos,budget_iter,num_arcs)], coef1);

        if(i > 0)
        {
          double coef2 = max(route_limit - vertex_deadline,vertex_info.nominal_service_time_ + Rn[i])*x_values[arc_pos];
          obj -= operator*(dual_vars.v_1_[dual_vars.v_1_var_to_index(arc_pos,budget_iter,num_arcs)], coef2);
        }

        double coef3 = (vertex_info.nominal_service_time_ + arc->distance())*x_values[arc_pos];
        obj -= operator*(dual_vars.v_2_[dual_vars.v_2_var_to_index(arc_pos,budget_iter,num_arcs)], coef3);

        if(budget_iter > 0)
        {
          double coef4 = (vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance())*x_values[arc_pos];
          obj -= operator*(dual_vars.v_3_[dual_vars.v_3_var_to_index(arc_pos,budget_iter,num_arcs)], coef4);
        }
      }
    }
  }

  if(combine_feas_op_cuts)
    obj -= operator*(*dual_bound_value,dual_vars.v_dual_bound_);
}

// This routine separates Benders' cuts violated by the current x solution.
// Violated cuts are found by solving the worker LP
//
IloBool SeparateBendersCutBaseline(IloEnv& master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar& dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex& worker_cplex, DualVariablesBaseline* dual_vars, const Instance & instance, IloObjective& worker_obj, IloExpr& cut_expr, bool combine_feas_op_cuts, Solution<double>& solution)
{
  // std::cout << "try to separate benders cut" << std::endl; 
  IloBool violatedCutFound = IloFalse;
  IloEnv worker_env = worker_cplex.getEnv();
  IloModel worker_model = worker_cplex.getModel();

  // Update the objective function in the worker LP.
  worker_model.remove(worker_obj);
  IloExpr obj_expr = worker_obj.getExpr();
  obj_expr.clear();
  FillObjectiveExpressionDualCompactBaselineContinuousSpace(obj_expr,*dual_vars,x_values,y_values,dual_bound_value,instance, combine_feas_op_cuts);
  worker_obj.setExpr(obj_expr);
  worker_model.add(worker_obj);
  obj_expr.end(); 

  // for(IloNum i = 0; i < y.getSize(); ++i)
  // {
  //   std::cout << y[i].getName() << " " << y_values[i] << std::endl;
  // }

  // for(IloNum i = 0; i < x.getSize(); ++i)
  // {
  //   std::cout << x[i].getName() << " " << x_values[i] << std::endl;
  // }
  
  // worker_cplex.exportModel("worker_model_Benders_compact_baseline.lp");
  // getchar(); getchar();

  // Solve the worker LP
  worker_cplex.solve();

  auto status = worker_cplex.getStatus();
  //std::cout << status << std::endl;

  if (status == IloAlgorithm::InfeasibleOrUnbounded || status == IloAlgorithm::Infeasible )
    throw "Infeasible Benders' subproblem";

  // Get the violated cut as an unbounded ray of the worker LP
  if ( status == IloAlgorithm::Unbounded)
  {
    ++solution.num_benders_feas_cuts_;
    //std::cout << "found new feasibility cut " << instance.num_feas_cuts_ << std::endl;
    IloNumVarArray var(worker_env);
    IloNumArray val(worker_env);

    worker_cplex.getRay(val, var);
    IloInt var_id = 0;

    for(IloInt i = 0; i < val.getSize(); ++i)
    {
      var_id = var[i].getId();
      cut_expr += operator*(val[i],dual_vars->ext_ray_coef_[var_id]);
      //std::cout << var[i].getName() << " " << val[i] << " " << dual_vars.ext_ray_coef[var_id] << std::endl;
    }
    // std::cout << cut_expr << std::endl;
    // getchar();getchar();

    var.end();
    val.end();

    return IloTrue;
  }

  // if the problem is solved to optimality and the current master solution is violated, a new optimality cut is found.
  if (status == IloAlgorithm::Optimal)
  {
    double new_dual_bound = worker_cplex.getObjValue();
    std::cout << status << " " << new_dual_bound << " x " << dual_bound_value << std::endl;
    // if combining feas and opt cuts, it suffices to have a bound greater than zero to imply a violated cut.
    if ( !combine_feas_op_cuts && double_greater(dual_bound_value,new_dual_bound) || combine_feas_op_cuts && !double_equals(new_dual_bound,0))
    {
      ++solution.num_benders_opt_cuts_;
      //std::cout << status << " " << new_dual_bound << " x " << dual_bound_value << std::endl;
      //std::cout << "found new optimality cut " << instance.num_opt_cuts_ << std::endl;
      IloNumArray u_0_values(worker_env);
      worker_cplex.getValues(u_0_values, dual_vars->u_0_);
      IloNumArray u_1_values(worker_env);
      worker_cplex.getValues(u_1_values, dual_vars->u_1_);
      IloNumArray u_2_values(worker_env);
      worker_cplex.getValues(u_2_values, dual_vars->u_2_);
      IloNumArray u_3_values(worker_env);
      worker_cplex.getValues(u_3_values, dual_vars->u_3_);

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
          cut_expr -= operator*(u_0_values[dual_vars->u_0_var_to_index(j,budget_iter,num_vertices)]*(*graph)[0][j]->distance(),y[j]);

        for (int i = 0; i < num_vertices; ++i)
        {
          const auto vertex_info = vertices_info[i];
          double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
          if(i <= num_mandatory)
            cut_expr += u_1_values[dual_vars->u_1_var_to_index(i,budget_iter,num_vertices)]*route_limit;
          else
            cut_expr += operator*(u_1_values[dual_vars->u_1_var_to_index(i,budget_iter,num_vertices)]*vertex_deadline,y[i]);
        
          if (i > 0) // only profitable and mandatory.
          {
            for(auto j: graph->AdjVerticesOut(i))
            {
              GArc * arc = (*graph)[i][j];
              const int arc_pos = graph->pos(i,j);
              IloExpr coef1(master_env);
              coef1 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + arc->distance(),1-x[arc_pos]) - vertex_info.nominal_service_time_ - arc->distance());
              cut_expr += operator*(u_2_values[dual_vars->u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef1);
              coef1.end();
              if(budget_iter > 0)
              {
                IloExpr coef2(master_env);
                coef2 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance(),1-x[arc_pos]) - vertex_info.nominal_service_time_  - vertex_info.dev_service_time_ - arc->distance());
                cut_expr += operator*(u_3_values[dual_vars->u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)], coef2);
                coef2.end();
              }
            }
          }
        }
      }

      if(combine_feas_op_cuts)
        cut_expr -= operator*(worker_cplex.getValue(dual_vars->u_dual_bound_),dual_bound);
      else
        cut_expr -= dual_bound;

      // for(IloNum i = 0; i < dual_vars->u_0_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->u_0_[i].getName() << " " << u_0_values[i] << std::endl;
      // }

      // for(IloNum i = 0; i < dual_vars->u_1_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->u_1_[i].getName() << " " << u_1_values[i] << std::endl;
      // }

      // for(IloNum i = 0; i < dual_vars->u_2_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->u_2_[i].getName() << " " << u_2_values[i] << std::endl;
      // }

      // for(IloNum i = 0; i < dual_vars->u_3_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->u_3_[i].getName() << " " << u_3_values[i] << std::endl;
      // }

      //std::cout << cut_expr << std::endl;
      // getchar();getchar();

      u_0_values.end();
      u_1_values.end();
      u_2_values.end();
      u_3_values.end();

      return IloTrue;
    }
  }

  return violatedCutFound;

} // END separate

// This routine separates Benders' cuts violated by the current x solution.
// Violated cuts are found by solving the worker LP
//
IloBool SeparateBendersCutSingleCommodity(IloEnv& master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar& dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex& worker_cplex, DualVariablesSingleCommodity* dual_vars, const double* R0, const double* Rn, const Instance & instance,  IloObjective& worker_obj, IloExpr& cut_expr, bool combine_feas_op_cuts, Solution<double>& solution)
{
  //std::cout << "try to separate benders cut" << std::endl; 
  IloBool violatedCutFound = IloFalse;
  IloEnv worker_env = worker_cplex.getEnv();
  IloModel worker_model = worker_cplex.getModel();

  // Update the objective function in the worker LP.
  worker_model.remove(worker_obj);
  IloExpr obj_expr = worker_obj.getExpr();
  obj_expr.clear();
  FillObjectiveExpressionDualCompactSingleCommodityContinuousSpace(obj_expr,*dual_vars,x_values,y_values,dual_bound_value,R0,Rn,instance,combine_feas_op_cuts);
  worker_obj.setExpr(obj_expr);
  worker_model.add(worker_obj);
  obj_expr.end(); 

  // for(IloNum i = 0; i < y.getSize(); ++i)
  // {
  //   std::cout << y[i].getName() << " " << y_values[i] << std::endl;
  // }

  // for(IloNum i = 0; i < x.getSize(); ++i)
  // {
  //   std::cout << x[i].getName() << " " << x_values[i] << std::endl;
  // }
  
  //worker_cplex.exportModel("worker_model_Benders_compact_single_commodity.lp");
  // getchar(); getchar();

  // Solve the worker LP
  worker_cplex.solve();

  auto status = worker_cplex.getStatus();
  //std::cout << status << std::endl;

  if (status == IloAlgorithm::InfeasibleOrUnbounded || status == IloAlgorithm::Infeasible )
    throw "Infeasible Benders' subproblem";

  // Get the violated cut as an unbounded ray of the worker LP
  if ( status == IloAlgorithm::Unbounded)
  {
    ++solution.num_benders_feas_cuts_;
    //std::cout << "found new feasibility cut " << solution.num_benders_feas_cuts_ << std::endl;
    IloNumVarArray var(worker_env);
    IloNumArray val(worker_env);

    worker_cplex.getRay(val, var);
    IloInt var_id = 0;

    for(IloInt i = 0; i < val.getSize(); ++i)
    {
      var_id = var[i].getId();
      cut_expr += operator*(val[i],dual_vars->ext_ray_coef_[var_id]);
      //std::cout << var[i].getName() << " " << val[i] << " " << dual_vars.ext_ray_coef[var_id] << std::endl;
    }

    // std::cout << cut_expr << std::endl;
    // getchar();getchar();

    var.end();
    val.end();

    return IloTrue;
  }

  // if the problem is solved to optimality and the current master solution is violated, a new optimality cut is found.
  if (status == IloAlgorithm::Optimal)
  {
    double new_dual_bound = worker_cplex.getObjValue();
    //std::cout << status << " " << new_dual_bound << " x " << dual_bound_value << std::endl;
    // if combining feas and opt cuts, it suffices to have a bound greater than zero to imply a violated cut.
    if ( !combine_feas_op_cuts && double_greater(dual_bound_value,new_dual_bound) || combine_feas_op_cuts && !double_equals(new_dual_bound,0))
    {
      ++solution.num_benders_opt_cuts_;
      //std::cout << status << " " << new_dual_bound << " x " << dual_bound_value << std::endl;
      //std::cout << "found new optimality cut " << solution.num_benders_opt_cuts_ << std::endl;
      IloNumArray v_0_values(worker_env);
      worker_cplex.getValues(v_0_values, dual_vars->v_0_);
      IloNumArray v_1_values(worker_env);
      worker_cplex.getValues(v_1_values, dual_vars->v_1_);
      IloNumArray v_2_values(worker_env);
      worker_cplex.getValues(v_2_values, dual_vars->v_2_);
      IloNumArray v_3_values(worker_env);
      worker_cplex.getValues(v_3_values, dual_vars->v_3_);

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
        for (int i = 0; i < num_vertices; ++i)
        {
          const auto vertex_info = vertices_info[i];
          double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);

          for (int j : graph->AdjVerticesOut(i))
          {
            const GArc * arc = (*graph)[i][j];
            const int arc_pos = graph->pos(i,j);
            double coef1 = (route_limit - R0[i] - vertex_info.nominal_service_time_ - arc->distance())*v_0_values[dual_vars->v_0_var_to_index(arc_pos,budget_iter,num_arcs)];
            cut_expr += operator*(coef1, x[arc_pos]);

            if(i > 0)
            {
              double coef2 = max(route_limit - vertex_deadline,vertex_info.nominal_service_time_ + Rn[i])* v_1_values[dual_vars->v_1_var_to_index(arc_pos,budget_iter,num_arcs)];
              cut_expr -= operator*(coef2, x[arc_pos]);
            }

            double coef3 = (vertex_info.nominal_service_time_ + arc->distance())*v_2_values[dual_vars->v_2_var_to_index(arc_pos,budget_iter,num_arcs)];
            cut_expr -= operator*(coef3, x[arc_pos]);

            if(budget_iter > 0)
            {
              double coef4 = (vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance())* v_3_values[dual_vars->v_3_var_to_index(arc_pos,budget_iter,num_arcs)];
              cut_expr -= operator*(coef4, x[arc_pos]);
            }
          }
        }
      }

      if(combine_feas_op_cuts)
      {
        double v_dual_bound = worker_cplex.getValue(dual_vars->v_dual_bound_);
        cut_expr -= operator*(v_dual_bound,dual_bound);
      }else
        cut_expr -= dual_bound;

      // for(IloNum i = 0; i < dual_vars->v_0_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->v_0_[i].getName() << " " << v_0_values[i] << std::endl;
      // }

      // for(IloNum i = 0; i < dual_vars->v_1_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->v_1_[i].getName() << " " << v_1_values[i] << std::endl;
      // }

      // for(IloNum i = 0; i < dual_vars->v_2_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->v_2_[i].getName() << " " << v_2_values[i] << std::endl;
      // }

      // for(IloNum i = 0; i < dual_vars->v_3_.getSize(); ++i)
      // {
      //   std::cout << dual_vars->v_3_[i].getName() << " " << v_3_values[i] << std::endl;
      // }

      // std::cout << cut_expr << std::endl;
      // getchar();getchar();

      v_0_values.end();
      v_1_values.end();
      v_2_values.end();
      v_3_values.end();

      return IloTrue;
    }
  }

  return violatedCutFound;

} // END separate


DualVariables::~DualVariables()
{
  for (auto& it: ext_ray_coef_)
    it.second.end();

  ext_ray_coef_.clear();
}

DualVariablesBaseline::DualVariablesBaseline(IloEnv& env,const Instance & instance, bool combine_feas_op_cuts):DualVariables(instance,combine_feas_op_cuts)
{
  const Graph * graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = instance.uncertainty_budget();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));

  // budget + 1 to consider level 0 of budget! 0,..., budget
  // u_0 defined for all i \in S U P, 0...budget
  u_0_ = IloNumVarArray(env, (budget+1)*(num_vertices-1), 0, IloInfinity, ILOFLOAT);
  // u_1 defined for all i \in N, 0...budget
  u_1_ = IloNumVarArray(env, (budget+1)*num_vertices, 0, IloInfinity, ILOFLOAT);
  // u_2 defined for all arcs (i,j) i != 0, 0...budget
  u_2_ = IloNumVarArray(env, (budget+1)*(num_arcs-num_arcs_from_origin), 0, IloInfinity, ILOFLOAT);
  // u_3 defined for all arcs (i,j) i != 0, 1...budget (considering 0...budget -1)
  u_3_ = IloNumVarArray(env, budget*(num_arcs-num_arcs_from_origin), 0, IloInfinity, ILOFLOAT);
  if(combine_feas_op_cuts)
    u_dual_bound_ = IloNumVar(env,0,IloInfinity,ILOFLOAT);
}

DualVariablesBaseline::~DualVariablesBaseline()
{
  u_0_.end();
  u_1_.end();
  u_2_.end();
  u_3_.end();

  if(combine_feas_op_cuts_)
    u_dual_bound_.end();
}

void DualVariablesBaseline::AddNamesToDualVariables()
{
  const Graph* graph = instance_.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = instance_.uncertainty_budget();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));

  // add name to variables.
  // u_0.
  for(int i = 1; i < num_vertices; ++i)
  {
    for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
    {
      char strnum[28];
      sprintf(strnum,"u_0(%d)(%d)",i,budget_iter);
      u_0_[u_0_var_to_index(i,budget_iter,num_vertices)].setName(strnum);
    }
  }

  // u_1.
  for(int i = 0; i < num_vertices; ++i)
  {
    for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
    {
      char strnum[28];
      sprintf(strnum,"u_1(%d)(%d)",i,budget_iter);
      u_1_[u_1_var_to_index(i,budget_iter,num_vertices)].setName(strnum);
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
        u_2_[u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].setName(strnum);
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
        u_3_[u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].setName(strnum);
      }
    }
  }

  if(combine_feas_op_cuts_)
    u_dual_bound_.setName("u_dual_bound");
}

void DualVariablesBaseline::SetDualVariablesProperties(IloEnv& master_env, MasterVariables& master_vars, const double* R0, const double* Rn)
{
  assert(R0 == nullptr);
  assert(Rn == nullptr);
  const Graph * graph = instance_.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto * vertices_info = graph->vertices_info();
  const int budget = instance_.uncertainty_budget();
  const double route_limit = instance_.limit();
  const int num_mandatory = instance_.num_mandatory();
  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for(auto j: graph->AdjVerticesOut(0))
    {
      IloExpr exp1(master_env);
      exp1 -= operator*((*graph)[0][j]->distance(),master_vars.y[j]);
      ext_ray_coef_[u_0_[u_0_var_to_index(j,budget_iter,num_vertices)].getId()] = exp1;
    }

    for (int i = 0; i < num_vertices; ++i)
    {
      const auto vertex_info = vertices_info[i];
      double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);
      
      IloExpr exp2(master_env);
      if(i <= num_mandatory)
      {
        exp2 += route_limit;
        ext_ray_coef_[u_1_[u_1_var_to_index(i,budget_iter,num_vertices)].getId()] = exp2;
      }else
      {
        exp2 +=operator*(vertex_deadline,master_vars.y[i]);
        ext_ray_coef_[u_1_[u_1_var_to_index(i,budget_iter,num_vertices)].getId()] = exp2;
      }
      if (i > 0) // only profitable and mandatory.
      {
        for(auto j: graph->AdjVerticesOut(i))
        {
          GArc * arc = (*graph)[i][j];
          const int arc_pos = graph->pos(i,j);
          IloExpr exp3(master_env);
          exp3 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + arc->distance(),1-master_vars.x[arc_pos]) - vertex_info.nominal_service_time_ - arc->distance());
          ext_ray_coef_[u_2_[u_2_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].getId()] = exp3;
          if(budget_iter > 0)
          {
            IloExpr exp4(master_env);
            exp4 += (operator*(vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance(),1-master_vars.x[arc_pos]) - vertex_info.nominal_service_time_  - vertex_info.dev_service_time_ - arc->distance());
            ext_ray_coef_[u_3_[u_3_var_to_index(arc_pos,budget_iter,num_arcs_from_origin,num_arcs)].getId()] = exp4;
          }
        }
      }
    }
  }

  // won´t be necessary, because, when combine_feas_op_cuts_, the subproblem always has a bounded optimal solution, so no feasibility cut is separated.
  // if(combine_feas_op_cuts_)
  // {
  //   IloExpr exp5(master_env);
  //   exp5 -= master_vars.dual_bound;
  //   ext_ray_coef_[u_dual_bound_.getId()] = exp5;
  // }
}

int DualVariablesBaseline::u_0_var_to_index(int vertex, int budget, int num_vertices)
{
  assert(num_vertices > 1);
  return budget * (num_vertices-1) + vertex - 1;
} 

int DualVariablesBaseline::u_1_var_to_index(int vertex, int budget, int num_vertices)
{
  return budget * num_vertices + vertex;
}

int DualVariablesBaseline::u_2_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs)
{
  assert(num_arcs > num_arcs_from_zero);
  return budget * (num_arcs-num_arcs_from_zero) + arc_pos - num_arcs_from_zero;
}

int DualVariablesBaseline::u_3_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs)
{
  assert(num_arcs > num_arcs_from_zero);
  assert(budget > 0);
  return (budget-1) * (num_arcs-num_arcs_from_zero) + arc_pos - num_arcs_from_zero;
}

DualVariablesSingleCommodity::DualVariablesSingleCommodity(IloEnv& env,const Instance & instance, bool combine_feas_op_cuts):DualVariables(instance,combine_feas_op_cuts)
{
  const Graph * graph = instance.graph();
  const int num_arcs = graph->num_arcs();
  const int budget = instance.uncertainty_budget();

  // budget + 1 to consider level 0 of budget! 0,..., budget.
  // u_0 defined for all (i,j) \in A, 0...budget.
  v_0_ = IloNumVarArray(env, (budget+1)*num_arcs, 0, IloInfinity, ILOFLOAT);
  // u_1 defined for all (i,j) \in A, 0...budget.
  v_1_ = IloNumVarArray(env, (budget+1)*num_arcs, 0, IloInfinity, ILOFLOAT);
  // u_2 defined for all (i,j) \in A, 0...budget.
  v_2_ = IloNumVarArray(env, (budget+1)*num_arcs, 0, IloInfinity, ILOFLOAT);
  // u_3 defined for all (i,j) \in A, 1...budget.
  v_3_ = IloNumVarArray(env, budget*num_arcs, 0, IloInfinity, ILOFLOAT);
  if(combine_feas_op_cuts)
    v_dual_bound_ = IloNumVar(env,0,IloInfinity,ILOFLOAT);
}

DualVariablesSingleCommodity::~DualVariablesSingleCommodity()
{
  v_0_.end();
  v_1_.end();
  v_2_.end();
  v_3_.end();

  if(combine_feas_op_cuts_)
    v_dual_bound_.end();
}

void DualVariablesSingleCommodity::AddNamesToDualVariables()
{
  const Graph* graph = instance_.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = instance_.uncertainty_budget();

  // add name to variables.
  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for(int i = 0; i < num_vertices; ++i)
    {
      for(int j = i+1; j < num_vertices; ++j)
      {
          GArc* arc = (*graph)[i][j];
          GArc* arc_inv = (*graph)[j][i];
          if(arc != nullptr)
          {
            char strnum[41], strnum1[41], strnum2[41];
            sprintf(strnum,"v_0(%d)(%d)(%d)",i,j,budget_iter);
            sprintf(strnum1,"v_1(%d)(%d)(%d)",i,j,budget_iter);
            sprintf(strnum2,"v_2(%d)(%d)(%d)",i,j,budget_iter);
            
            v_0_[v_0_var_to_index(graph->pos(i,j),budget_iter,num_arcs)].setName(strnum);
            v_1_[v_1_var_to_index(graph->pos(i,j),budget_iter,num_arcs)].setName(strnum1);
            v_2_[v_2_var_to_index(graph->pos(i,j),budget_iter,num_arcs)].setName(strnum2);

            if(budget_iter > 0)
            {
              char strnum3[41];
              sprintf(strnum3,"v_3(%d)(%d)(%d)",i,j,budget_iter);
              v_3_[v_3_var_to_index(graph->pos(i,j),budget_iter,num_arcs)].setName(strnum3);
            }
          }

          if(arc_inv != nullptr)
          {
            char strnum[41], strnum1[41], strnum2[41];
            sprintf(strnum,"v_0(%d)(%d)(%d)",j,i,budget_iter);
            sprintf(strnum1,"v_1(%d)(%d)(%d)",j,i,budget_iter);
            sprintf(strnum2,"v_2(%d)(%d)(%d)",j,i,budget_iter);
            v_0_[v_0_var_to_index(graph->pos(j,i),budget_iter,num_arcs)].setName(strnum);
            v_1_[v_1_var_to_index(graph->pos(j,i),budget_iter,num_arcs)].setName(strnum1);
            v_2_[v_2_var_to_index(graph->pos(j,i),budget_iter,num_arcs)].setName(strnum2);

            if(budget_iter > 0)
            {
              char strnum3[41];
              sprintf(strnum3,"v_3(%d)(%d)(%d)",j,i,budget_iter);
              v_3_[v_3_var_to_index(graph->pos(j,i),budget_iter,num_arcs)].setName(strnum3);
            }
          }
      }
    }
  }

  if(combine_feas_op_cuts_)
    v_dual_bound_.setName("v_dual_bound");
}

void DualVariablesSingleCommodity::SetDualVariablesProperties(IloEnv& master_env, MasterVariables& master_vars, const double* R0, const double* Rn)
{
  const Graph * graph = instance_.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto * vertices_info = graph->vertices_info();
  const int budget = instance_.uncertainty_budget();
  const double route_limit = instance_.limit();
  const int num_mandatory = instance_.num_mandatory();

  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for (int i = 0; i < num_vertices; ++i)
    {
      const auto vertex_info = vertices_info[i];
      double vertex_deadline = (i <= num_mandatory)? route_limit: round_decimals(vertex_info.profit_/vertex_info.decay_ratio_,2);

      for (int j : graph->AdjVerticesOut(i))
      {
        const GArc * arc = (*graph)[i][j];
        const int arc_pos = graph->pos(i,j);
        IloExpr exp1(master_env);
        exp1  += operator*(master_vars.x[arc_pos],route_limit - R0[i] - vertex_info.nominal_service_time_ - arc->distance());
        ext_ray_coef_[v_0_[v_0_var_to_index(arc_pos,budget_iter,num_arcs)].getId()] = exp1;

        if(i > 0)
        {
          IloExpr exp2(master_env);
          exp2 -= operator*(master_vars.x[arc_pos],max(route_limit - vertex_deadline, vertex_info.nominal_service_time_ + Rn[i]));
          ext_ray_coef_[v_1_[v_1_var_to_index(arc_pos,budget_iter,num_arcs)].getId()] = exp2;
        }

        IloExpr exp3(master_env);
        exp3 -= operator*(master_vars.x[arc_pos], vertex_info.nominal_service_time_ + arc->distance());
        ext_ray_coef_[v_2_[v_2_var_to_index(arc_pos,budget_iter,num_arcs)].getId()] = exp3;

        if(budget_iter > 0)
        {
          IloExpr exp4(master_env);
          exp4 -= operator*(master_vars.x[arc_pos], vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + arc->distance());
          ext_ray_coef_[v_3_[v_3_var_to_index(arc_pos,budget_iter,num_arcs)].getId()] = exp4;
        }
      }
    }
  }

  // won´t be necessary, because, when combine_feas_op_cuts_, the subproblem always has a bounded optimal solution, so no feasibility cut is separated.
  // if(combine_feas_op_cuts_)
  // {
  //   IloExpr exp5(master_env);
  //   exp5 -= master_vars.dual_bound;
  //   ext_ray_coef_[this->v_dual_bound_.getId()] = exp5;
  // }
}

int DualVariablesSingleCommodity::v_0_var_to_index(int arc_pos, int budget, int num_arcs)
{
  return budget * num_arcs + arc_pos;
} 

int DualVariablesSingleCommodity::v_1_var_to_index(int arc_pos, int budget, int num_arcs)
{
  return budget * num_arcs + arc_pos;
}

int DualVariablesSingleCommodity::v_2_var_to_index(int arc_pos, int budget, int num_arcs)
{
  return budget * num_arcs + arc_pos;
}

int DualVariablesSingleCommodity::v_3_var_to_index(int arc_pos, int budget, int num_arcs)
{
  assert(budget > 0);
  return (budget-1) * num_arcs + arc_pos;
}


void PopulateByRowDualCompactSingleCommodityContinuousSpaceBaseline(IloEnv& env, IloModel& model, DualVariablesSingleCommodity * dual_vars, const Instance& instance, bool combine_feas_op_cuts)
{
  const Graph* graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_mandatory = instance.num_mandatory();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto* vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();

  for (int i = 0; i < num_vertices; ++i)
  {
    for (int j : graph->AdjVerticesOut(i))
    {
      assert((*graph)[i][j] != nullptr);
      const int arc_pos = graph->pos(i,j);
      for(int budget_iter = 0; budget_iter <= budget; ++budget_iter)
      {
        IloExpr exp(env);

        exp += dual_vars->v_0_[dual_vars->v_0_var_to_index(arc_pos,budget_iter,num_arcs)];
        exp -= dual_vars->v_1_[dual_vars->v_1_var_to_index(arc_pos,budget_iter,num_arcs)];
        exp += dual_vars->v_2_[dual_vars->v_2_var_to_index(arc_pos,budget_iter,num_arcs)];

        for (int k: graph->AdjVerticesOut(j))
          exp -= dual_vars->v_2_[dual_vars->v_2_var_to_index(graph->pos(j,k),budget_iter,num_arcs)];
        
        if (budget_iter > 0)
          exp += dual_vars->v_3_[dual_vars->v_3_var_to_index(arc_pos,budget_iter,num_arcs)];
        
        if (budget_iter < budget)
        {
          for (int k: graph->AdjVerticesOut(j))
            exp -= dual_vars->v_3_[dual_vars->v_3_var_to_index(graph->pos(j,k),budget_iter+1,num_arcs)];
        }

        const auto coef = (budget_iter == budget && j > num_mandatory)? vertices_info[j].decay_ratio_: 0.0;

        // if combine optimality and feasibility cuts, add contibution of additional dual var.
        combine_feas_op_cuts? exp -= operator*(coef,dual_vars->v_dual_bound_) : exp-= coef;
        model.add(exp >= 0);
        exp.end();

        if(j == 0)
          dual_vars->v_1_[dual_vars->v_1_var_to_index(arc_pos,budget_iter,num_arcs)].setUB(0.0);

        if(i == 0)
        {
          dual_vars->v_2_[dual_vars->v_2_var_to_index(arc_pos,budget_iter,num_arcs)].setUB(0.0);
          if(budget_iter > 0)
            dual_vars->v_3_[dual_vars->v_3_var_to_index(arc_pos,budget_iter,num_arcs)].setUB(0.0);
        }
      }
    }
  }

  // add extra constraint to scale/normalize the values of the dual variables.
  if(combine_feas_op_cuts)
  {
    double num_dual_vars = 1.0 + dual_vars->v_0_.getSize() + dual_vars->v_1_.getSize() + dual_vars->v_2_.getSize()+ dual_vars->v_3_.getSize();
    double normalized_coef = 1.0/sqrt(1.0*num_dual_vars); // all of these variables have coefficient equal to 1, so the sum of the powers of 2 is equal to their sum,
    
    IloExpr exp_normalize(env);
    // add all dual variables to the expression.
    for(IloInt var_index = 0; var_index < dual_vars->v_0_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->v_0_[var_index]);
    for(IloInt var_index = 0; var_index < dual_vars->v_1_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->v_1_[var_index]);
    for(IloInt var_index = 0; var_index < dual_vars->v_2_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->v_2_[var_index]);
    for(IloInt var_index = 0; var_index < dual_vars->v_3_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->v_3_[var_index]);
    exp_normalize += operator*(normalized_coef,dual_vars->v_dual_bound_);

  
    //std::cout << num_dual_vars << " " <<  normalized_coef << std::endl;

    model.add(exp_normalize == 1);
    exp_normalize.end();
  }
}

void PopulateByRowDualCompactBaselineContinuousSpaceBaseline(IloEnv& env, IloModel& model, DualVariablesBaseline * dual_vars, const Instance& instance, bool combine_feas_op_cuts)
{
  const Graph* graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_mandatory = instance.num_mandatory();
  const int num_arcs = graph->num_arcs();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));
  const auto* vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();

  IloExpr exp(env);

  exp += dual_vars->u_1_[dual_vars->u_1_var_to_index(0,0,num_vertices)];
  for (auto j: graph->AdjVerticesIn(0))
    exp -= dual_vars->u_2_[dual_vars->u_2_var_to_index(graph->pos(j,0),0, num_arcs_from_origin, num_arcs)];

  model.add(exp >= 0);
  exp.end();

  for(int budget_iter = 1; budget_iter <= budget; ++budget_iter)
  {
    IloExpr exp(env);
    exp += dual_vars->u_1_[dual_vars->u_1_var_to_index(0,budget_iter,num_vertices)];

    for (auto j: graph->AdjVerticesIn(0))
    {
      exp -= dual_vars->u_2_[dual_vars->u_2_var_to_index(graph->pos(j,0),budget_iter, num_arcs_from_origin, num_arcs)];
      exp -= dual_vars->u_3_[dual_vars->u_3_var_to_index(graph->pos(j,0),budget_iter, num_arcs_from_origin, num_arcs)];
    }
    model.add(exp >= 0);
    exp.end();
  }

  for(int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for (int i = 1; i < num_vertices; ++i)
    {
      IloExpr exp(env);

      exp -= dual_vars->u_0_[dual_vars->u_0_var_to_index(i,budget_iter,num_vertices)];
      exp += dual_vars->u_1_[dual_vars->u_1_var_to_index(i,budget_iter,num_vertices)];

      for(auto j: graph->AdjVerticesOut(i))
      {
        exp += dual_vars->u_2_[dual_vars->u_2_var_to_index(graph->pos(i,j),budget_iter,num_arcs_from_origin,num_arcs)];
        if(budget_iter < budget)
          exp += dual_vars->u_3_[dual_vars->u_3_var_to_index(graph->pos(i,j),budget_iter+1,num_arcs_from_origin,num_arcs)];
      }

      for(auto j: graph->AdjVerticesIn(i))
      {
        if (j != 0)
        {
          exp -= dual_vars->u_2_[dual_vars->u_2_var_to_index(graph->pos(j,i),budget_iter,num_arcs_from_origin,num_arcs)];
          if (budget_iter > 0)
            exp -= dual_vars->u_3_[dual_vars->u_3_var_to_index(graph->pos(j,i),budget_iter,num_arcs_from_origin,num_arcs)];
        }
      }

      const auto coef = (budget_iter == budget && i > num_mandatory)? - vertices_info[i].decay_ratio_: 0.0;
      
      // if combine optimality and feasibility cuts, add contibution of additional dual var.
      combine_feas_op_cuts? exp -= operator*(coef,dual_vars->u_dual_bound_) : exp-= coef;
      model.add(exp >= 0);
      exp.end();


      if(!(*graph)[0][i])
        dual_vars->u_0_[dual_vars->u_0_var_to_index(i,budget_iter,num_vertices)].setUB(0.0);
    }
  }

  // add extra constraint to scale/normalize the values of the dual variables.
  if(combine_feas_op_cuts)
  {
    double num_dual_vars = 1.0 + dual_vars->u_0_.getSize() + dual_vars->u_1_.getSize() + dual_vars->u_2_.getSize()+ dual_vars->u_3_.getSize();
    double normalized_coef = 1.0/sqrt(1.0*num_dual_vars); // all of these variables have coefficient equal to 1, so the sum of the powers of 2 is equal to their sum,
    
    IloExpr exp_normalize(env);
    // add all dual variables to the expression.
    for(IloInt var_index = 0; var_index < dual_vars->u_0_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->u_0_[var_index]);
    for(IloInt var_index = 0; var_index < dual_vars->u_1_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->u_1_[var_index]);
    for(IloInt var_index = 0; var_index < dual_vars->u_2_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->u_2_[var_index]);
    for(IloInt var_index = 0; var_index < dual_vars->u_3_.getSize(); ++var_index)
      exp_normalize += operator*(normalized_coef,dual_vars->u_3_[var_index]);
    exp_normalize += operator*(normalized_coef,dual_vars->u_dual_bound_);

  
    //std::cout << num_dual_vars << " " <<  normalized_coef << std::endl;

    model.add(exp_normalize == 1);
    exp_normalize.end();
  }
}

WorkerI::WorkerI(IloEnv& master_env, const Instance& instance, MasterVariables& master_vars, bool combine_feas_op_cuts, bool export_model)
    : worker_vars_(nullptr), worker_cplex_(worker_env_), worker_obj_(worker_env_), combine_feas_op_cuts_(combine_feas_op_cuts) 
{}

WorkerBaseline::WorkerBaseline(IloEnv& master_env, const Instance& instance, MasterVariables& master_vars, bool combine_feas_op_cuts, bool export_model)
    : WorkerI(master_env,instance,master_vars,combine_feas_op_cuts,export_model)
    {

    worker_vars_ = new DualVariablesBaseline(worker_env_,instance,combine_feas_op_cuts);
    // create the worker subproblem!
    // The subproblem will be always the same, except for the objective function, which relies
    // on the values of X and y variables. Then, we can create a unique model beforehand.
    worker_cplex_.setOut(worker_env_.getNullStream());

    // Turn off the presolve reductions and set the CPLEX optimizer
    // to solve the worker LP with primal simplex method.
    worker_cplex_.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
    worker_cplex_.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

    IloModel worker_model(worker_env_);
    worker_cplex_.extract(worker_model);

    // Set empty initial objective function.
    worker_obj_ = IloObjective(worker_env_);
    worker_obj_.setName("obj1");

    worker_obj_.setSense(IloObjective::Minimize);

    worker_model.add(worker_obj_);

    worker_vars_->SetDualVariablesProperties(master_env,master_vars,nullptr,nullptr);

    PopulateByRowDualCompactBaselineContinuousSpaceBaseline(worker_env_,worker_model,static_cast<DualVariablesBaseline*>(worker_vars_),instance,combine_feas_op_cuts);

    if(export_model)
      worker_vars_->AddNamesToDualVariables();
}

WorkerSingleCommodity::WorkerSingleCommodity(IloEnv& master_env, const Instance& instance, MasterVariables& master_vars, const double* R0, const double* Rn, bool combine_feas_op_cuts, bool export_model)
    : WorkerI(master_env,instance,master_vars,combine_feas_op_cuts,export_model)
    {

    worker_vars_ = new DualVariablesSingleCommodity(worker_env_,instance,combine_feas_op_cuts);
    // create the worker subproblem!
    // The subproblem will be always the same, except for the objective function, which relies
    // on the values of X and y variables. Then, we can create a unique model beforehand.
    worker_cplex_.setOut(worker_env_.getNullStream());

    // Turn off the presolve reductions and set the CPLEX optimizer
    // to solve the worker LP with primal simplex method.
    worker_cplex_.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
    worker_cplex_.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

    IloModel worker_model(worker_env_);
    worker_cplex_.extract(worker_model);

    // Set empty initial objective function.
    worker_obj_ = IloObjective(worker_env_);
    worker_obj_.setName("obj1");

    worker_obj_.setSense(IloObjective::Minimize);

    worker_model.add(worker_obj_);

    worker_vars_->SetDualVariablesProperties(master_env,master_vars,R0,Rn);

    PopulateByRowDualCompactSingleCommodityContinuousSpaceBaseline(worker_env_,worker_model,static_cast<DualVariablesSingleCommodity*>(worker_vars_),instance,combine_feas_op_cuts);

    if(export_model)
      worker_vars_->AddNamesToDualVariables();
}