  #include <cmath>
  #include <algorithm>
  #include "local_branching.h"
  #include "CommodityFormulationsForm.h"
  #include "graph_algorithms.h"
  #include "general.h"

  LocalBranching::LocalBranching()
  {
  }

  LocalBranching::~LocalBranching()
  {
    this->Reset();
  }

  void LocalBranching::Reset()
  {
  	(this->env_).end();
  }

  void LocalBranching::Init(Instance& instance)
  {
  	this->curr_instance_ = &instance;
  	Graph * graph = instance.graph();
    int num_vertices = graph->num_vertices();
    int num_arcs = graph->num_arcs();
    int num_routes = instance.num_vehicles();

    (this->solution_).Reset(num_vertices,num_arcs,num_routes);

  	this->env_ = IloEnv();
  	this->cplex_ = IloCplex(this->env_);
  	this->cplex_.setOut(this->env_.getNullStream());
    this->model_ = IloModel(this->env_);

    this->slack_ = IloNumVar(this->env_, 0, num_routes, ILOFLOAT);
  	f = IloNumVarArray(this->env_, num_arcs, 0, IloInfinity, ILOFLOAT);
  	x = IloNumVarArray(this->env_, num_arcs, 0.0, 1.0, ILOINT);
    y = IloNumVarArray(this->env_, num_vertices, 0.0, 1.0, ILOINT);

    this->lb_constraint_ =  IloRange(this->env_,0, y[0], 10);
    (this->model_).add(this->lb_constraint_);

    double * R0 = Dijkstra(graph,false,false);
    double * Rn = Dijkstra(graph,true,false);

    std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

    if(K_LOCAL_BRANCHING_ADD_CUTS)
    {
      std::list<UserCutGeneral*> * cuts = NULL;
      Solution<int> sol;

      //(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
      //(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;

      cuts = CapacitatedSingleCommodity(*(this->curr_instance_),R0,Rn,-1,&sol,true,false,true,NULL,NULL);

      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = false;
      (*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = false;
      (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = false;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = false;

      if (!(sol.is_feasible_))
      {
        (this->solution_).is_infeasible_ = true;
      }else
      {
        CapacitatedSingleCommodityForm(this->env_, this->model_, this->slack_, this->y, this->x, this->f, *(this->curr_instance_), R0, Rn, false);

        // add user cuts as restrictions here!
        AddCutsAsConstraints(cuts);
      }

      DeleteCuts(cuts);
      delete cuts;
      cuts = NULL;
    }else
    {
      // to explude from model these inequalities
      //(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;

      CapacitatedSingleCommodityForm(this->env_, this->model_, this->slack_, this->y, this->x, this->f, *(this->curr_instance_), R0, Rn, false);

      //(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = false;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = false;
    }

    //cplex_.extract(model_);
    //cplex_.exportModel("teste.lp");
    //getchar(); getchar();

    delete [] R0;
    delete [] Rn;
    R0 = NULL;
    Rn = NULL;
  }

  void LocalBranching::BuildHeuristicSolution()
  {
  /*  Graph * graph = (this->curr_instance_)->graph();
    int num_vertices = graph->num_vertices(), v1 = 0, v2 = 0, curr_route_index = 0;
    int last_vertex_added = 0;
    std::vector<std::list<int>> * arcs = graph->adj_lists_;
    VertexStatus * status = NULL;
    Route * curr_route = NULL;
    GArc * curr_arc = NULL;
    size_t cont = 0;

    std::list<int> q;
    q.push_back(0);

    do
    {
      v1 = q.front();
      q.pop_front();

      // end current route
      if(v1 == num_vertices - 1)
      {
          curr_route = &(((this->solution_).routes_vec_)[curr_route_index]);
          curr_arc = (*graph)[last_vertex_added][v1];
          if(curr_arc == NULL) throw "Inappropriate addition of vertex to route";
          else
          {
            (curr_route->time_) += curr_arc->dist();
          }

          last_vertex_added = 0;
          ++curr_route_index;
      }else if(v1 != 0)
      {
        // add vertex to current route
        curr_route = &(((this->solution_).routes_vec_)[curr_route_index]);
        status = &(((this->solution_).vertex_status_vec_)[v1]);

        // remove from list of unvisited_vertices
        ((this->solution_).unvisited_vertices_).erase(status->pos_);

        // adds vertex to route
        status->selected_ = true;
        status->route_ = curr_route_index;
        status->pos_  = (curr_route->vertices_).insert((curr_route->vertices_).end(), v1);

        ((this->solution_).profits_sum_) += (graph->profits())[v1];
        (curr_route->sum_profits_) += (graph->profits())[v1];

        curr_arc = (*graph)[last_vertex_added][v1];
        if(curr_arc == NULL) throw "Inappropriate addition of vertex to route";
        else
        {
          (curr_route->time_) += curr_arc->dist();
        }
        last_vertex_added = v1;
      }

      for(std::list<int>::iterator it = ((*arcs)[v1]).begin(); it != ((*arcs)[v1]).end(); it++)
      {
        v2 = *it;

        if((this->curr_int_x_)[graph->pos(v1,v2)])
        {
            cont++;
            q.push_front(v2);
        }
      }
    }while(!(q.empty()));

    if( ((this->solution_).unvisited_vertices_).empty() ) (this->solution_).is_optimal_ = true;

    if((cont != (this->curr_int_x_).count()) || (!((this->solution_).CheckCorrectness(*(this->curr_instance_)))) )throw "Error Building FP solution";
    */
  }

boost::dynamic_bitset<> * LocalBranching::Run(double time_limit, HeuristicSolution * initial_sol)
  {
    Timestamp * ti = NewTimestamp();
    Timer * timer = GetTimer();
    timer->Clock(ti);

    if(!(initial_sol->is_feasible_) || ((this->solution_).is_infeasible_))
    {
      (this->solution_).is_infeasible_ = true;
      (this->solution_).time_ += timer->CurrentElapsedTime(ti);
      return NULL;
    }

  	std::cout << std::setprecision(2) << std::fixed;

    Graph * graph = (this->curr_instance_)->graph();
    int num_mandatory = graph->num_mandatory();
    int num_vertices = graph->num_vertices();
    int num_arcs = graph->num_arcs();
    boost::dynamic_bitset<> * curr_bitset = NULL;

    Solution<int> sol;
    // main algorithm
    std::cout << initial_sol->profits_sum_ << std::endl;
    //std::cout << (initial_sol->bitset_arcs_).count() << std::endl;

    BuildLeftLBConstraintOnX(initial_sol->bitset_arcs_,K_LOCAL_BRANCHING_K);

    //cplex_.extract(model_);
    //cplex_.exportModel("teste.lp");
    //getchar(); getchar();

    optimizeCapacitatedSingleCommodity(this->cplex_,this->env_,this->model_,y,x,f,*(this->curr_instance_),sol,NULL,time_limit,false,false,false, NULL, NULL,NULL,initial_sol);

    std::cout << sol.lb_ << std::endl;

    if( (sol.is_feasible_) && double_greater(sol.lb_, 1.0*(initial_sol->profits_sum_)) )
    {
      // BUILD BITSET OF CURRENT SOLUTION
      curr_bitset = new boost::dynamic_bitset<>(num_arcs,0);
      RetrieveAndRoundArcValuesBasic(*curr_bitset);
      std::cout << "* " << sol.lb_ << std::endl;
    }

    std::cout << " - " << timer->CurrentElapsedTime(ti) << "s" << std::endl;

    return curr_bitset;

    /*bool can_stop = false;
    do
    {
      //IloModel model(this->env_);
      //sol.reset();
      //std::cout << "iter" << std::endl;
      // generate and solve left branch
      BuildLeftBranchOnX(model,R0,Rn,K_LOCAL_BRANCHING_K);

      //cplex_.extract(model);
      //cplex_.exportModel("teste.lp");
      //getchar(); getchar();

      if(!double_equals(time_limit,-1.0)) time_limit = std::max(0.0, time_limit - timer->CurrentElapsedTime(ti));

      optimizeCapacitatedSingleCommodity(this->cplex_,this->env_,model,y,x,f,*(this->curr_instance_),sol,NULL,30,false,false,false, R0, Rn,cuts,&(this->initial_sol_));

      std::cout << sol.lb_ << std::endl;
      can_stop = true;

      if( (!(sol.is_feasible_)) || (!(double_greater(sol.lb_, 1.0*(((this->solutions_vec_).back()).second)))) ) can_stop = true;
      else
      {
        // BUILD BITSET OF CURRENT SOLUTION AND UPDATE SOLUTIONS_VEC
        RetrieveAndRoundArcValuesBasic(curr_bitset);
        (this->solutions_vec_).push_back(std::pair<boost::dynamic_bitset<>,int>(curr_bitset,sol.lb_));
        std::cout << "* " << sol.lb_ << std::endl;
      }
      model.end();
    }while(!can_stop);

    if(!((this->solution_).is_infeasible_)) BuildHeuristicSolution();

    delete [] R0;
    delete [] Rn;
    R0 = NULL;
    Rn = NULL;*/
  }

  /*void LocalBranching::BuildLeftBranchOnX(IloModel& model, double * R0, double * Rn, int k)
  {
    CapacitatedSingleCommodityForm(this->env_, model, y, x, f, *(this->curr_instance_), R0, Rn, false);

    // add local branching constraints
    for(size_t i = 0; i < (this->solutions_vec_).size() -1 ; ++i)
    {
      BuildRightLBConstraintOnX(model,((this->solutions_vec_)[i]).first,k);
    }

    BuildLeftLBConstraintOnX(model,((this->solutions_vec_).back()).first,k);
  }*/

  void LocalBranching::BuildLeftLBConstraintOnX(boost::dynamic_bitset<>& bitset, int k)
  {
    int num_arcs = ((this->curr_instance_)->graph())->num_arcs();

    IloExpr exp(this->env_);
    for(int i = 0; i < num_arcs; ++i)
    {
        if(bitset[i] == 1) exp += (1 - x[i]);
        else exp += x[i];
    }

    (this->lb_constraint_).setExpr(exp);
    (this->lb_constraint_).setUB(k);

    exp.end();
  }

  /*void LocalBranching::BuildRightLBConstraintOnX(IloModel& model,boost::dynamic_bitset<>& bitset, int k)
  {
    int num_arcs = ((this->curr_instance_)->graph())->num_arcs();

    IloExpr exp(this->env_);
    for(int i = 0; i < num_arcs; ++i)
    {
        if(bitset[i] == 1) exp += (1 - x[i]);
        else exp += x[i];
    }

    model.add(exp >= k +1);
    exp.end();
  }

  void LocalBranching::BuildLeftLBConstraintOnY(IloModel& model,boost::dynamic_bitset<>& bitset, int k)
  {
    int num_vertices = ((this->curr_instance_)->graph())->num_vertices();

    IloExpr exp(this->env_);
    for(int i = 0; i < num_vertices; ++i)
    {
        if(bitset[i] == 1) exp += (1 - y[i]);
        else exp += y[i];
    }

    model.add(exp <= k);
    exp.end();
  }

  void LocalBranching::BuildRightLBConstraintOnY(IloModel& model,boost::dynamic_bitset<>& bitset, int k)
  {
    int num_vertices = ((this->curr_instance_)->graph())->num_vertices();

    IloExpr exp(this->env_);
    for(int i = 0; i < num_vertices; ++i)
    {
        if(bitset[i] == 1) exp += (1 - y[i]);
        else exp += y[i];
    }

    model.add(exp >= k +1);
    exp.end();
  }*/

  void LocalBranching::RetrieveAndRoundArcValuesBasic(boost::dynamic_bitset<> & bitset)
  {
    IloNumArray curr_values_x(this->env_);
    (this->cplex_).getValues(curr_values_x,x);
    Graph * graph = (this->curr_instance_)->graph();
    int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
    std::vector<std::list<int>> * arcs = graph->adj_lists_;
    int v1 = 0, v2 = 0;

    std::list<int> q;
    q.push_front(0);
    std::vector<bool> visited_nodes(num_vertices,false);
    visited_nodes[0] = true;
    int curr_pos = 0;

    bitset.reset();

    // we don't have to check every arc: due to the flow conservation constraints, only perform width first search
    do
    {
      v1 = q.front();
      q.pop_front();

      for(std::list<int>::iterator it = ((*arcs)[v1]).begin(); it != ((*arcs)[v1]).end(); it++)
      {
          v2 = *it;
          curr_pos = graph->pos(v1,v2);

          if(!double_equals(curr_values_x[curr_pos],0.0))
          {
            if(!(visited_nodes[v2]))
            {
                q.push_front(v2);
                visited_nodes[v2] = true;
            }

            //if(double_greater((this->curr_relax_x_)[curr_pos],0.0)) (this->curr_int_x_)[curr_pos] = 1;
            if(round((curr_values_x)[curr_pos])) bitset[curr_pos] = 1;
          }
      }
    }while(!q.empty());
  }

void LocalBranching::AddCutsAsConstraints(std::list<UserCutGeneral*> * initial_cuts)
{
  Graph* graph = (this->curr_instance_)->graph();
  if((initial_cuts != NULL)&&(!(initial_cuts->empty())))
  {
    //IloRangeArray root_cuts(env);
    for(std::list<UserCutGeneral*>::iterator it = initial_cuts->begin(); it != initial_cuts->end(); it++)
    {
      switch((*it)->type_)
      {
        case K_TYPE_GCC_CUT:
        {
          UserCut * curr_user_cut = static_cast<UserCut*>((*it));
          IloExpr exp(this->env_);
          for(std::list<std::pair<int,int>>::iterator it2 = curr_user_cut->lhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->lhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            int v1 = (*it2).first;
            int v2 = (*it2).second;
            //if(!solve_relax) std::cout << "(" << v1 << "," << v2 << ")" << std::endl;
            exp += operator*(x[graph->pos(v1,v2)],(curr_user_cut->lhs_)[graph->pos(v1,v2)]);
            //if(!solve_relax) std::cout << "2" << std::endl;
          }

          for(std::list<int>::iterator it2 = curr_user_cut->rhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->rhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            exp -= y[*it2];
          }


          (this->model_).add(exp >= 0);
          //solution.set_cut_added((*it)->type_,false);

          exp.end();
          break;
        }
        case K_TYPE_CONFLICT_CUT:
        {
          UserCut * curr_user_cut = static_cast<UserCut*>((*it));
          IloExpr exp(this->env_);
          for(std::list<std::pair<int,int>>::iterator it2 = curr_user_cut->lhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->lhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            int v1 = (*it2).first;
            int v2 = (*it2).second;
            exp += operator*(x[graph->pos(v1,v2)],(curr_user_cut->lhs_)[graph->pos(v1,v2)]);
          }

          for(std::list<int>::iterator it2 = curr_user_cut->rhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->rhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            exp -= y[*it2];
          }

          (this->model_).add(exp >= 0);
          //solution.set_cut_added((*it)->type_,false);

          exp.end();
          break;
        }
        case K_TYPE_COVER_BOUND_CUT:
        {
          CoverBoundCut * curr_cover_bound_cut = static_cast<CoverBoundCut*>((*it));

          IloExpr exp(this->env_);
          for(std::list<int>::iterator it2 = curr_cover_bound_cut->cut_elements_.begin(); it2 != curr_cover_bound_cut->cut_elements_.end(); it2++)
          {
            exp += operator*((curr_cover_bound_cut->cut_exp_coefs_)[*it2],y[*it2]);
          }

          (this->model_).add(exp <= curr_cover_bound_cut->rhs_);

          //solution.set_cut_added((*it)->type_,false);

          exp.end();
          break;
        }

        case K_TYPE_CLIQUE_CONFLICT_CUT:
        {
          UserCut * curr_user_cut = static_cast<UserCut*>((*it));
          IloExpr exp(this->env_);
          for(std::list<std::pair<int,int>>::iterator it2 = curr_user_cut->lhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->lhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            int v1 = (*it2).first;
            int v2 = (*it2).second;
            exp += operator*(x[graph->pos(v1,v2)],(curr_user_cut->lhs_)[graph->pos(v1,v2)]);
          }

          for(std::list<int>::iterator it2 = curr_user_cut->rhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->rhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            exp -= y[*it2];
          }

          (this->model_).add(exp >= 0);
          //solution.set_cut_added((*it)->type_,false);

          exp.end();
          break;
        }
        default:
        {
          std::cout << (*it)->type_ << std::endl;
          throw 5;
        }
      }

  //delete (*it);
  //*it = NULL;
    }
  }
}
