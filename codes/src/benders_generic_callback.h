#pragma once

#include <optional>
#include <ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include "src/instance.h"
#include "src/general.h"
#include "src/solution.hpp"

struct MasterVariables
{
   IloNumVarArray x;
   IloNumVarArray y;
   IloNumVar slack;
   IloNumVar dual_bound;
};

class DualVariables
{
public:
   explicit DualVariables(const Instance &instance, bool combine_feas_op_cuts) : instance_(instance), combine_feas_op_cuts_(combine_feas_op_cuts) {}
   virtual ~DualVariables();

   virtual void AddNamesToDualVariables() = 0;
   virtual void SetDualVariablesProperties(IloEnv &master_env, MasterVariables &master_vars, const double *R0, const double *Rn) = 0;
   std::unordered_map<IloInt, IloExpr> ext_ray_coef_;
   bool combine_feas_op_cuts_;
   const Instance &instance_;
};

class DualVariablesBaseline : public DualVariables
{
public:
   explicit DualVariablesBaseline(IloEnv &env, const Instance &instance, bool combine_feas_op_cuts);
   virtual ~DualVariablesBaseline();
   void AddNamesToDualVariables() final;
   void SetDualVariablesProperties(IloEnv &master_env, MasterVariables &master_vars, const double *R0, const double *Rn) final;
   int u_0_var_to_index(int vertex, int budget, int num_vertices);
   int u_1_var_to_index(int vertex, int budget, int num_vertices);
   int u_2_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs);
   int u_3_var_to_index(int arc_pos, int budget, int num_arcs_from_zero, int num_arcs);

   IloNumVarArray u_0_;
   IloNumVarArray u_1_;
   IloNumVarArray u_2_;
   IloNumVarArray u_3_;
   IloNumVar u_dual_bound_;
};

class DualVariablesSingleCommodity : public DualVariables
{
public:
   explicit DualVariablesSingleCommodity(IloEnv &env, const Instance &instance, bool combine_feas_op_cuts);
   virtual ~DualVariablesSingleCommodity();
   void AddNamesToDualVariables() final;
   void SetDualVariablesProperties(IloEnv &master_env, MasterVariables &master_vars, const double *R0, const double *Rn) final;
   int v_0_var_to_index(int arc_pos, int budget, int num_arcs);
   int v_1_var_to_index(int arc_pos, int budget, int num_arcs);
   int v_2_var_to_index(int arc_pos, int budget, int num_arcs);
   int v_3_var_to_index(int arc_pos, int budget, int num_arcs);

   IloNumVarArray v_0_;
   IloNumVarArray v_1_;
   IloNumVarArray v_2_;
   IloNumVarArray v_3_;
   IloNumVar v_dual_bound_;
};

ILOSTLBEGIN

void FillObjectiveExpressionDualCompactBaselineContinuousSpace(IloExpr &obj, DualVariablesBaseline &dual_vars, IloNumArray &x_values, IloNumArray &y_values, std::optional<IloNum> dual_bound_value, const Instance &instance, bool combine_feas_op_cuts);
IloBool SeparateBendersCutBaseline(IloEnv &master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar &dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex &worker_cplex, DualVariablesBaseline *dual_vars, const Instance &instance, IloObjective &worker_obj, IloExpr &cut_expr, bool combine_feas_op_cuts, Solution<double> &solution);
void PopulateByRowDualCompactBaselineContinuousSpace(IloEnv &env, IloModel &model, DualVariablesBaseline *dual_vars, const Instance &instance, bool combine_feas_op_cuts);

void FillObjectiveExpressionDualCompactSingleCommodityContinuousSpace(IloExpr &obj, DualVariablesSingleCommodity &dual_vars, IloNumArray &x_values, IloNumArray &y_values, std::optional<IloNum> dual_bound_value, const double *R0, const double *Rn, const Instance &instance, bool combine_feas_op_cuts);
IloBool SeparateBendersCutSingleCommodity(IloEnv &master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar &dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex &worker_cplex, DualVariablesSingleCommodity *dual_vars, const double *R0, const double *Rn, const Instance &instance, IloObjective &worker_obj, IloExpr &cut_expr, bool combine_feas_op_cuts, Solution<double> &solution);
void PopulateByRowDualCompactSingleCommodityContinuousSpace(IloEnv &env, IloModel &model, DualVariablesSingleCommodity *dual_vars, const Instance &instance, bool combine_feas_op_cuts);

/** The Benders' worker thread-local class. */
class WorkerI
{
public:
   // The constructor sets up the IloCplex algorithm to solve the worker LP, and
   // creates the worker LP.
   explicit WorkerI(IloEnv &master_env, const Instance &instance, MasterVariables &master_vars, bool combine_feas_op_cuts, bool export_model);

   virtual ~WorkerI()
   {
      worker_env_.end();
   }

   IloEnv &worker_env()
   {
      return worker_env_;
   }
   IloCplex &worker_cplex()
   {
      return worker_cplex_;
   }
   DualVariables *worker_vars()
   {
      return worker_vars_;
   }
   IloObjective &worker_obj()
   {
      return worker_obj_;
   }

protected:
   IloEnv worker_env_;
   IloCplex worker_cplex_;
   DualVariables *worker_vars_ = nullptr;
   IloObjective worker_obj_;
   bool combine_feas_op_cuts_;
};

/** The Benders' worker thread-local class. */
class WorkerBaseline : public WorkerI
{
public:
   // The constructor sets up the IloCplex algorithm to solve the worker LP, and
   // creates the worker LP.
   explicit WorkerBaseline(IloEnv &master_env, const Instance &instance, MasterVariables &master_vars, bool combine_feas_op_cuts, bool export_model);

   ~WorkerBaseline()
   {
      if (worker_vars_ != nullptr)
      {
         delete worker_vars_;
         worker_vars_ = nullptr;
      }
   }
};

/** The Benders' worker thread-local class. */
class WorkerSingleCommodity : public WorkerI
{
public:
   // The constructor sets up the IloCplex algorithm to solve the worker LP, and
   // creates the worker LP.
   explicit WorkerSingleCommodity(IloEnv &master_env, const Instance &instance, MasterVariables &master_vars, const double *R0, const double *Rn, bool combine_feas_op_cuts, bool export_model);

   ~WorkerSingleCommodity()
   {
      if (worker_vars_ != nullptr)
      {
         delete worker_vars_;
         worker_vars_ = nullptr;
      }
   }
};

// Implementation class for the user-defined user cut callback.
// The function BendersUserCallback allows to add Benders' cuts as user cuts.
//
class BendersGenericCallbackI : public IloCplex::Callback::Function
{
public:
   BendersGenericCallbackI(IloEnv &master_env, Instance &instance, MasterVariables &master_vars, bool combine_feas_op_cuts, bool export_model, Solution<double> &solution, IloInt num_workers = 1)
       : master_env_(master_env), instance_(instance), master_vars_(master_vars), combine_feas_op_cuts_(combine_feas_op_cuts), export_model_(export_model), solution_(solution), workers_(num_workers, nullptr) {}

   ~BendersGenericCallbackI()
   {
      IloInt num_workers = workers_.size();
      for (IloInt w = 0; w < num_workers; ++w)
      {
         if (workers_[w] != nullptr)
         {
            delete workers_[w];
            workers_[w] = nullptr;
         }
      }
      workers_.clear();
   }

   virtual WorkerI *AllocateNewWorker() = 0;
   virtual IloBool SeparateBendersCut(IloEnv &master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar &dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex &worker_cplex, DualVariables *dual_vars, const Instance &instance, IloObjective &worker_obj, IloExpr &cut_expr, bool combine_feas_op_cuts, Solution<double> &solution) = 0;

   void invoke(const IloCplex::Callback::Context &context) ILO_OVERRIDE
   {
      int const thread_num = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

      // setup
      if (context.inThreadUp())
      {
         delete (workers_[thread_num]);
         workers_[thread_num] = AllocateNewWorker();
         return;
      }

      // teardown
      if (context.inThreadDown())
      {
         delete workers_[thread_num];
         workers_[thread_num] = nullptr;
         return;
      }

      IloEnv master_env = context.getEnv();
      IloNumArray x_values(master_env);
      IloNumArray y_values(master_env);
      IloNum dual_bound_value;

      // Get the current master solution
      switch (context.getId())
      {
      case IloCplex::Callback::Context::Id::Candidate:
      {
         if (!context.isCandidatePoint()) // The master model is always bounded
            throw IloCplex::Exception(-1, "Unbounded solution");
         context.getCandidatePoint(master_vars_.x, x_values);
         context.getCandidatePoint(master_vars_.y, y_values);
         dual_bound_value = context.getCandidatePoint(master_vars_.dual_bound);
         break;
      }
      case IloCplex::Callback::Context::Id::Relaxation:
      {
         // const int num_nodes = context.getIntInfo(IloCplex::Callback::Context::Info::NodeCount);
         const int id_node = context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID);
         // only separate benders cuts related to the linear relaxation in the root node.
         if (id_node == 0)
         {
            context.getRelaxationPoint(master_vars_.x, x_values);
            context.getRelaxationPoint(master_vars_.y, y_values);
            dual_bound_value = context.getRelaxationPoint(master_vars_.dual_bound);
            break;
         }
         else
         {
            // Free memory.
            x_values.end();
            y_values.end();
            return;
         }
      }
      default:
      {
         // Free memory.
         x_values.end();
         y_values.end();
         throw IloCplex::Exception(-1, "Unexpected contextID");
      }
      }

      auto curr_worker = workers_[thread_num];
      // std::cout << "current master obj " << context.getCandidateObjective() << std::endl;
      // std::cout << "current dual bound " << dual_bound_value << std::endl;
      // Benders' cut separation.
      IloExpr cut_expr(master_env);
      IloBool sep_status = SeparateBendersCut(master_env, master_vars_.x, master_vars_.y, master_vars_.dual_bound, x_values, y_values, dual_bound_value, curr_worker->worker_cplex(), curr_worker->worker_vars(), instance_, curr_worker->worker_obj(), cut_expr, combine_feas_op_cuts_, solution_);
      if (sep_status)
      {
         // std::cout << "Adicionou corte " << solution_.num_benders_feas_cuts_ << std::endl;
         // getchar(); getchar();
         // Add the cut
         IloRange cut(master_env, 0, cut_expr, IloInfinity);

         switch (context.getId())
         {
         case IloCplex::Callback::Context::Id::Candidate:
            context.rejectCandidate(cut);
            break;
         case IloCplex::Callback::Context::Id::Relaxation:
            context.addUserCut(cut,
                               IloCplex::UseCutPurge,
                               IloFalse);
            break;
         default:
            cut.end();
            throw IloCplex::Exception(-1, "Unexpected contextID");
         }
         cut.end();
      }
      else
         cut_expr.end();
   }

protected:
   IloEnv &master_env_;
   const Instance &instance_;
   Solution<double> &solution_;
   MasterVariables &master_vars_;
   bool combine_feas_op_cuts_;
   bool export_model_;
   std::vector<WorkerI *> workers_;
};

class BendersGenericCallbackBaseline : public BendersGenericCallbackI
{
public:
   BendersGenericCallbackBaseline(IloEnv &master_env, Instance &instance, MasterVariables &master_vars, bool combine_feas_op_cuts, bool export_model, Solution<double> &solution, IloInt num_workers = 1)
       : BendersGenericCallbackI(master_env, instance, master_vars, combine_feas_op_cuts, export_model, solution, num_workers) {}

   virtual ~BendersGenericCallbackBaseline() = default;

   virtual WorkerI *AllocateNewWorker() final
   {
      return new WorkerBaseline(master_env_, instance_, master_vars_, combine_feas_op_cuts_, export_model_);
   }
   virtual IloBool SeparateBendersCut(IloEnv &master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar &dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex &worker_cplex, DualVariables *dual_vars, const Instance &instance, IloObjective &worker_obj, IloExpr &cut_expr, bool combine_feas_op_cuts, Solution<double> &solution) final
   {
      return SeparateBendersCutBaseline(master_env, master_vars_.x, master_vars_.y, master_vars_.dual_bound, x_values, y_values, dual_bound_value, worker_cplex, static_cast<DualVariablesBaseline *>(dual_vars), instance_, worker_obj, cut_expr, combine_feas_op_cuts_, solution_);
   }
};

class BendersGenericCallbackSingleCommodity : public BendersGenericCallbackI
{
public:
   BendersGenericCallbackSingleCommodity(IloEnv &master_env, Instance &instance, MasterVariables &master_vars, bool combine_feas_op_cuts, bool export_model, Solution<double> &solution, const double *R0, const double *Rn, IloInt num_workers = 1)
       : BendersGenericCallbackI(master_env, instance, master_vars, combine_feas_op_cuts, export_model, solution, num_workers), R0_(R0), Rn_(Rn) {}

   virtual ~BendersGenericCallbackSingleCommodity() = default;

   virtual WorkerI *AllocateNewWorker() final
   {
      return new WorkerSingleCommodity(master_env_, instance_, master_vars_, R0_, Rn_, combine_feas_op_cuts_, export_model_);
   }
   virtual IloBool SeparateBendersCut(IloEnv &master_env, IloNumVarArray &x, IloNumVarArray &y, IloNumVar &dual_bound, IloNumArray &x_values, IloNumArray &y_values, IloNum dual_bound_value, IloCplex &worker_cplex, DualVariables *dual_vars, const Instance &instance, IloObjective &worker_obj, IloExpr &cut_expr, bool combine_feas_op_cuts, Solution<double> &solution) final
   {
      return SeparateBendersCutSingleCommodity(master_env, master_vars_.x, master_vars_.y, master_vars_.dual_bound, x_values, y_values, dual_bound_value, worker_cplex, static_cast<DualVariablesSingleCommodity *>(dual_vars), R0_, Rn_, instance_, worker_obj, cut_expr, combine_feas_op_cuts_, solution_);
   }

   const double *R0_;
   const double *Rn_;
};