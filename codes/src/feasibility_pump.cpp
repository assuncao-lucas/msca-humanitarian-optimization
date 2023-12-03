  #include <cmath>
  #include <algorithm>
  #include "feasibility_pump.h"
  #include "CommodityFormulationsForm.h"
  #include "graph_algorithms.h"
  #include "general.h"

  bool compare_func3(const std::pair<int,double>& v1, const std::pair<int,double>& v2)
  {
    return double_greater(v1.second,v2.second);
  }

  FeasibilityPump::FeasibilityPump()
  {
  	this->original_obj_norm_ = 0.0;
  	this->curr_alpha_ = 0.0;
    this->previous_alpha_ = std::numeric_limits<double>::infinity();
  	this->found_int_y_ = false;
  	this->found_int_x_ = false;
  	this->previous_integrality_gap_ = 0.0;
  	this->curr_integrality_gap_ = 0.0;
  }

  FeasibilityPump::~FeasibilityPump()
  {
    this->Reset();
  }

  void FeasibilityPump::Reset()
  {
  	(this->cplex_).end();
  	(this->model_).end();
  	(this->env_).end();
  }

  void FeasibilityPump::Init(Instance& instance)
  {
  	this->curr_instance_ = &instance;
  	const Graph * graph = instance.graph();
    int num_vertices = graph->num_vertices();
    int num_arcs = graph->num_arcs();
    int num_routes = instance.num_vehicles();
    (this->solution_).Reset(num_vertices,num_arcs,instance.num_vehicles());

  	this->curr_alpha_ = 0.0;

  	this->env_ = IloEnv();
  	this->model_ = IloModel(this->env_);
  	this->cplex_ = IloCplex(this->env_);
  	this->cplex_.setOut(this->env_.getNullStream());

  	this->curr_relax_y_ =  IloNumArray(this->env_);
  	this->curr_relax_x_ = IloNumArray(this->env_);

  	this->curr_int_y_ = boost::dynamic_bitset<>(num_vertices,0);
  	this->curr_int_x_ = boost::dynamic_bitset<>(num_arcs,0);

  	double * R0 = Dijkstra(graph,false,false);
    double * Rn = Dijkstra(graph,true,false);

  	f = IloNumVarArray(this->env_, num_arcs, 0, IloInfinity, ILOFLOAT);
  	x = IloNumVarArray(this->env_, num_arcs, 0.0, 1.0, ILOFLOAT);
    y = IloNumVarArray(this->env_, num_vertices, 0.0, 1.0, ILOFLOAT);
    this->slack_ = IloNumVar(this->env_, 0, num_routes, ILOFLOAT);

  	BuildModel(R0,Rn);

  	delete [] R0;
  	delete [] Rn;
  	R0 = NULL;
  	Rn = NULL;
  }

  void FeasibilityPump::RetrieveAndRoundVertexValues()
  {
    cplex_.getValues(curr_relax_y_,y);
  	const Graph * graph = (this->curr_instance_)->graph();
  	int num_vertices = graph->num_vertices(), num_mandatory = (this->curr_instance_)->num_mandatory();
    this->previous_int_y_ = this->curr_int_y_;
    (this->curr_int_y_).reset();
    (this->curr_integrality_gaps_y_).clear();
  	this->found_int_y_ = true;

  	this->previous_integrality_gap_ = this->curr_integrality_gap_;
  	this->curr_integrality_gap_ = 0.0;

  	for(int i = num_mandatory + 1; i < num_vertices-1; i++)
  	{
      if(!double_equals(curr_relax_y_[i],0.0))
      {
        double curr_integrality_gap_y = 0.0;
    		if(round((this->curr_relax_y_)[i])) (this->curr_int_y_)[i] = 1;
    		if(!double_equals((this->curr_relax_y_)[i],1.0*((this->curr_int_y_)[i])))
    		{
    		    curr_integrality_gap_y = fabs((this->curr_relax_y_)[i] - 1.0*((this->curr_int_y_)[i]));
    		    this->found_int_y_ = false;
            (this->curr_integrality_gap_)+= curr_integrality_gap_y;
    		}

        this->curr_integrality_gaps_y_.push_back(std::pair<int,double>(i,curr_integrality_gap_y));
      }
  	}
  }

  bool FeasibilityPump::BuildPathsBFSIter(std::vector<std::list<std::pair<int,double>>> & residual_graph, std::list<std::list<int>>& paths,
  std::list<int> curr_path, std::vector<char>& colors, size_t vertex)
  {
  	if(colors[vertex] != 'w')
    {
  	  paths.push_back(curr_path);
  	 	return true;
    }

  	curr_path.push_back(vertex);
  	if(vertex == residual_graph.size() - 1)
  	{
  		paths.push_back(curr_path);
  		return true;
  	}else
  	{
  		colors[vertex] = 'g';
  		bool found_path =  false;
  		for(auto it = residual_graph[vertex].begin(); it != residual_graph[vertex].end(); ++it)
  		{
  			if(BuildPathsBFSIter(residual_graph,paths,curr_path,colors, it->first))
  			{
  				colors[vertex] = 'b';
  				found_path = true;
  			}
  		}
  		if(found_path) return true;
  	}
  	colors[vertex] = 'w';
  	return false;
  }

  void FeasibilityPump::BuildPathsBFS(std::vector<std::list<std::pair<int,double>>> & residual_graph, std::list<std::list<int>>& paths)
  {
  	std::list<int> curr_path;
  	std::vector<char> colors(residual_graph.size(),'w');
  	BuildPathsBFSIter(residual_graph, paths, curr_path, colors, 0);
  }

  void FeasibilityPump::RetrieveAndRoundArcValuesBFS()
  {
  	cplex_.getValues(curr_relax_x_,x);
  	const Graph * graph = (this->curr_instance_)->graph();
  	int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
  	std::vector<std::list<int>> * arcs = graph->adj_lists_;
  	int v1 = 0, v2 = 0;
  	std::list<std::list<int>> paths;
  	(this->curr_int_x_).reset();

  	std::vector<std::list<std::pair<int,double>>> residual_graph(num_vertices,std::list<std::pair<int,double>>());

  	std::list<int> q;
    q.push_front(0);
  	std::vector<bool> visited_nodes(num_vertices,false);
  	visited_nodes[0] = true;
  	int curr_pos = 0;

  	this->found_int_x_ = true;
  	this->previous_integrality_gap_ = this->curr_integrality_gap_;
  	this->curr_integrality_gap_ = 0.0;

  	// we don't have to check every arc: due to the flow conservation constraints, only perform width first search
  	do
  	{
  		v1 = q.front();
  		q.pop_front();

  		for(std::list<int>::iterator it = ((*arcs)[v1]).begin(); it != ((*arcs)[v1]).end(); it++)
  		{
  		    v2 = *it;
  		    curr_pos = graph->pos(v1,v2);

  		    if(!double_equals(curr_relax_x_[curr_pos],0.0))
  		    {
  					(residual_graph[v1]).push_back(std::pair<int,double>(v2,curr_relax_x_[curr_pos]));
  					if(!(visited_nodes[v2]))
  					{
  					    q.push_front(v2);
  					    visited_nodes[v2] = true;
  					}
  		    }
  		}
  	}while(!q.empty());

  	// sort each list in residual_graph in non-increasing order of relaxation value
  	for(int i = 0; i < num_vertices; i++)
  	{
  		(residual_graph[i]).sort(compare_func3);
  	}

  	// Call BFS to build the paths_list
  	BuildPathsBFS(residual_graph, paths);

  	int count = 0;
  	for(auto it = paths.begin(); it != paths.end(); ++it)
  	{
  		count++;
  		std::cout << "path " << count << std::endl;
  		for(auto it2 = it->begin(); it2 != it->end(); ++it2)
  		{
  			v2 = *it2;
  			if(it2 != it->begin())
  			{
  				curr_pos = graph->pos(v1,v2);
  				(this->curr_int_x_)[curr_pos] = 1;
  			}
  			v1 = v2;
  			std::cout << *it2 << " ";
  		}
  		std::cout << std::endl;
  		getchar();getchar();
  	}


  	this->curr_integrality_gap_ = 0.0;
  	for(curr_pos = 0; curr_pos < num_arcs; curr_pos++)
  	{
  		if(!double_equals((this->curr_relax_x_)[curr_pos],(this->curr_int_x_)[curr_pos]))
  		{
  			this->found_int_x_ = false;
  			(this->curr_integrality_gap_)+= fabs((this->curr_relax_x_)[curr_pos] - 1.0*((this->curr_int_x_)[curr_pos]));
  		}
  	}

  }

  void FeasibilityPump::RetrieveAndRoundArcValuesBasic()
  {
  	cplex_.getValues(curr_relax_x_,x);
  	Graph * graph = (this->curr_instance_)->graph();
  	int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
  	std::vector<std::list<int>> * arcs = graph->adj_lists_;
  	int v1 = 0, v2 = 0;

  	std::list<int> q;
    q.push_front(0);
  	std::vector<bool> visited_nodes(num_vertices,false);
  	visited_nodes[0] = true;
  	int curr_pos = 0;

  	this->found_int_x_ = true;
  	this->previous_integrality_gap_ = this->curr_integrality_gap_;
  	this->curr_integrality_gap_ = 0.0;

  	this->previous_int_x_ = this->curr_int_x_;
  	(this->curr_int_x_).reset();
    (this->curr_integrality_gaps_x_).clear();

  	// we don't have to check every arc: due to the flow conservation constraints, only perform width first search
  	do
  	{
  		v1 = q.front();
  		q.pop_front();

  		for(std::list<int>::iterator it = ((*arcs)[v1]).begin(); it != ((*arcs)[v1]).end(); it++)
  		{
  		    v2 = *it;
  		    curr_pos = graph->pos(v1,v2);

  		    if(!double_equals(curr_relax_x_[curr_pos],0.0))
  		    {
            double curr_integrality_gap_x = 0.0;
  					if(!(visited_nodes[v2]))
  					{
  					    q.push_front(v2);
  					    visited_nodes[v2] = true;
  					}

  					//if(double_greater((this->curr_relax_x_)[curr_pos],0.0)) (this->curr_int_x_)[curr_pos] = 1;
  					if(round((this->curr_relax_x_)[curr_pos])) (this->curr_int_x_)[curr_pos] = 1;
  					if(!double_equals((this->curr_relax_x_)[curr_pos],(this->curr_int_x_)[curr_pos]))
  					{
  						this->found_int_x_ = false;
              curr_integrality_gap_x = fabs((this->curr_relax_x_)[curr_pos] - 1.0*((this->curr_int_x_)[curr_pos]));
  						(this->curr_integrality_gap_)+= curr_integrality_gap_x;
  					}

            this->curr_integrality_gaps_x_.push_back(std::pair<int,double>(curr_pos,curr_integrality_gap_x));
  		    }
  		}
  	}while(!q.empty());
  }

  void FeasibilityPump::SetNewObjStage1()
  {
  	IloExpr new_obj(this->env_);
  	Graph * graph = (this->curr_instance_)->graph();
  	int num_vertices = graph->num_vertices(), num_mandatory = graph->num_mandatory();
  	//int * profits = graph->profits();
  	//this->curr_obj_const_ = 0;

  	for(int i = num_mandatory + 1; i < num_vertices-1; i++)
  	{
  		if((this->curr_int_y_)[i] == 0)
  		 {
  		 	new_obj += y[i];
  		 }
  		 else if((this->curr_int_y_)[i] == 1)
  		 {
  			new_obj += (1 - y[i]);
  		 	//(this->curr_obj_const_)++;
  		 }
  		 else
  		 {
  		    throw std::string("Found a binary var that, once rounded, is neither 0 nor 1. This shouldn't happen!");
  		 }
  	}

  	(this->objective_).setExpr(operator*((1.0-this->curr_alpha_)/sqrt(num_vertices - num_mandatory -2),new_obj) - operator*((this->curr_alpha_)/this->original_obj_norm_,this->original_obj_expr_));
  	new_obj.end();
  }

  void FeasibilityPump::SetNewObjStage2()
  {
  	IloExpr new_obj(this->env_);
  	Graph * graph = (this->curr_instance_)->graph();
  	std::vector<std::list<int>> * arcs = graph->adj_lists_;
  	int num_vertices = graph->num_vertices(), num_mandatory = graph->num_mandatory(), num_arcs = graph->num_arcs();
  	int v1 = 0, v2 = 0, curr_pos = 0;
  	//int * profits = graph->profits();
  	//this->curr_obj_const_ = 0;

  	for(int i = 0; i < num_vertices; i++)
  	{
  		v1 = i;
  		for(std::list<int>::iterator it = ((*arcs)[v1]).begin(); it != ((*arcs)[v1]).end(); it++)
  		{
  			v2 = *it;
  			curr_pos = graph->pos(v1,v2);
  			if((this->curr_int_x_)[curr_pos] == 0)
  			 {
  			 	new_obj += x[curr_pos];
  			 }
  			 else if((this->curr_int_x_)[curr_pos] == 1)
  			 {
  				new_obj += (1 - x[curr_pos]);
  			 	//(this->curr_obj_const_)++;
  			 }
  			 else
  			 {
  			    throw std::string("Found a binary var that, once rounded, is neither 0 nor 1. This shouldn't happen!");
  			 }
  		 }
  	}

  	(this->objective_).setExpr(operator*((1.0-this->curr_alpha_)/sqrt(num_arcs),new_obj) - operator*((this->curr_alpha_)/this->original_obj_norm_,this->original_obj_expr_));
  	new_obj.end();
  }

  void FeasibilityPump::BuildHeuristicSolution()
  {
    Graph * graph = (this->curr_instance_)->graph();
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

    (this->solution_).BuildBitset(*(this->curr_instance_));
    if((cont != (this->curr_int_x_).count()) || (!((this->solution_).CheckCorrectness(*(this->curr_instance_)))) )throw "Error Building FP solution";
  }

void FeasibilityPump::Run()
  {
    Timestamp * ti = NewTimestamp();
    Timer * timer = GetTimer();
    timer->Clock(ti);
  	std::cout << std::setprecision(2) << std::fixed;
  	double init_stalls_cycle_integrality_gap = std::numeric_limits<double>::infinity();
  	double curr_stalls_cycle_improvement = std::numeric_limits<double>::infinity();
  	this->previous_integrality_gap_ = this->curr_integrality_gap_ = std::numeric_limits<double>::infinity();
  	int stage_0_iter = 1, stage_1_iter = 0, stage_2_iter = 0, stalls_cycle_count = 0;
    int num_mandatory = ((this->curr_instance_)->graph())->num_mandatory();
    int num_vertices = ((this->curr_instance_)->graph())->num_vertices();
    int num_arcs = ((this->curr_instance_)->graph())->num_arcs();
    int num_flips_basis = 0, curr_num_flips = 0;
    if(K_FLIP_BASIS) num_flips_basis = K_FLIP_BASIS;
    else num_flips_basis = ceil(perturbation_flip_percentage*num_vertices);
    int curr_prob_of_flipping = 0;
    int num_perturbations_stage1 = 0, num_perturbations_stage2 = 0, num_restarts_stage1 = 0, num_restarts_stage2 = 0;
    (this->curr_alpha_) = initial_alpha_stage1;
  	// Stage 0: solve LP relaxation with original objective

    if(K_FEASIBILITY_PUMP_ADD_CUTS)
    {
      double * R0 = Dijkstra(((this->curr_instance_)->graph()),false,false);
      double * Rn = Dijkstra(((this->curr_instance_)->graph()),true,false);
      Solution<int> sol;

      std::list<UserCutGeneral*> * cuts = NULL;

      // **** these cuts are selected in BuildModel !!!!
      /*std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();

      //(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
      //(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;*/
      cuts = optimizeCapacitatedSingleCommodity(this->cplex_,this->env_,this->model_,this->y,this->x,this->f,*(this->curr_instance_),sol,NULL,-1,true,false,true, R0, Rn, NULL, NULL);
      DeleteCuts(cuts);
      delete cuts;
      cuts = NULL;
      delete [] R0;
      R0 = NULL;
      delete [] Rn;
      Rn = NULL;

      if (!(sol.is_feasible_))
      {
        (this->solution_).is_infeasible_ = true;
        (this->solution_).time_stage1_ = timer->CurrentElapsedTime(ti);
        return;
      }
    }else
    {
    	cplex_.extract(model_);
    	//cplex_.exportModel("FP.lp");
    	//getchar(); getchar();

      if (!(cplex_.solve()))
      {
        if((cplex_.getCplexStatus() == IloCplex::Infeasible)||(cplex_.getCplexStatus() == IloCplex::InfOrUnbd)) (this->solution_).is_infeasible_ = true;
        (this->solution_).time_stage1_ = timer->CurrentElapsedTime(ti);
        return;
      }
    }

  	this->RetrieveAndRoundVertexValues();

  	init_stalls_cycle_integrality_gap = this->curr_integrality_gap_;

  	// skips Stage 1 if all vertices are integer or K_SOLVE_STAGE_1 is false
  	if(!(this->found_int_y_)&& K_SOLVE_STAGE_1)
  	{
  		// Change objective to minimization
  		this->objective_.setSense(IloObjective::Minimize);

  		// Stage 1: only considers y variables in the new objective function
  		do
  		{
        this->previous_alpha_ = this->curr_alpha_;
  			(this->curr_alpha_)*= alpha_decrease_rate;
  			//std::cout << "alpha: " << this->curr_alpha_ << std::endl;
  			// Update objective function according to distance from current rounded (integer) y values
  			SetNewObjStage1();

  			if(double_less(fabs(this->curr_alpha_ - this->previous_alpha_),K_APLHA_DECREMENT_PRECISION)) stalls_cycle_count++;
        //if(stalls_cycle_count > 0) std::cout << stalls_cycle_count << std::endl;
  			stage_1_iter++;

  			// resolve
  			cplex_.extract(model_);
  			//cplex_.exportModel("FP.lp");
  			//getchar(); getchar();
  			if (!(cplex_.solve()))
        {
          if((cplex_.getCplexStatus() == IloCplex::Infeasible)||(cplex_.getCplexStatus() == IloCplex::InfOrUnbd)) (this->solution_).is_infeasible_ = true;

          (this->solution_).time_stage1_ = timer->CurrentElapsedTime(ti);

          (this->solution_).num_iterations_stage1_ = stage_1_iter + stage_0_iter;
          (this->solution_).num_perturbations_stage1_ = num_perturbations_stage1;
          (this->solution_).num_restarts_stage1_ = num_restarts_stage1;
          return;
        }
  			// update y values
  			this->RetrieveAndRoundVertexValues();

  			if(!(this->found_int_y_))
  			{
  				//std::cout << "curr_improv: " << fabs(this->previous_integrality_gap_ - this->curr_integrality_gap_) << std::endl;
  				curr_stalls_cycle_improvement = fabs(init_stalls_cycle_integrality_gap - this->curr_integrality_gap_);
  				//std::cout << "stalls_cycle_improv: " << curr_stalls_cycle_improvement << std::endl;
  				if(double_greater(curr_stalls_cycle_improvement,K_PUMP_IMPROVEMENT_TOLERANCE))
  				{
  					stalls_cycle_count = 0;
  					init_stalls_cycle_integrality_gap = this->curr_integrality_gap_;
  				}

          if (stalls_cycle_count >= max_stalls_stage1)
          {
            num_restarts_stage1++;
            // means that a "cycle" has been heuristically detected: PERFORM RESTART
            //std::cout << "ciclou" << std::endl;
            for(int i  = num_mandatory + 1; i < num_vertices - 1; i++)
            {
                if((this->curr_int_y_)[i] == (this->previous_int_y_)[i])
                {
                  //std::cout << (this->curr_relax_y_)[i] << " " << ((this->curr_int_y_)[i]) << std::endl;
                  //std::cout << fabs((this->curr_relax_y_)[i] - 1.0*((this->curr_int_y_)[i])) << std::endl;
                  curr_prob_of_flipping = ceil((fabs((this->curr_relax_y_)[i] - 1.0*((this->curr_int_y_)[i])) + 0.03)*100.0);
                  //std::cout << "flip prob: " << curr_prob_of_flipping << std::endl;
                  //int curr_rand = rand()%101;
                  //std::cout << "curr rand: " << curr_rand << std::endl;
                  if(rand()%101 <= curr_prob_of_flipping)
                  {
                    //std::cout << "flippou!" << std::endl;
                    ((this->curr_int_y_)[i]).flip();
                  }
                  //getchar();
                }
            }
            stalls_cycle_count = 0;
  					init_stalls_cycle_integrality_gap = this->curr_integrality_gap_;
            // if detected cycle of size one, PERFORM PERTURBATION
          }else if((stage_1_iter > 1) && (double_less(fabs(this->curr_alpha_ - this->previous_alpha_),K_APLHA_DECREMENT_PRECISION)) && (this->curr_int_y_ == this->previous_int_y_))
    			{
              curr_num_flips = num_flips_basis/2 + 1 + rand()%num_flips_basis;
              int cont = 0;
              num_perturbations_stage1++;
              //std::cout << "num_flips: " << curr_num_flips << std::endl;
              //getchar();getchar();
    				 //std::cout << " *** ciclou!" << std::endl;
             (this->curr_integrality_gaps_y_).sort(compare_func3);
             for(auto it = (this->curr_integrality_gaps_y_).begin(); it != (this->curr_integrality_gaps_y_).end(); ++it)
             {
               //std::cout << it->second << std::endl;
               //std::cout << (this->curr_int_y_[it->first]) << " - ";
               if(cont < curr_num_flips) ((this->curr_int_y_)[it->first]).flip();
               else break;
               //std::cout << (this->curr_int_y_[it->first]) << std::endl;
               //getchar();getchar();
               cont++;
             }
             //std::cout << std::endl;
    				 //getchar();getchar();
    			}
  			}

  		}while((!(this->found_int_y_))&&(stage_1_iter < max_iter_stage1));
  	}

  	this->RetrieveAndRoundArcValuesBasic();

    (this->solution_).time_stage1_ = timer->CurrentElapsedTime(ti);
    timer->Clock(ti);
    //std::cout << "iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
  	//if(this->found_int_y_) std::cout << " * Y INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
  	//if(this->found_int_x_) std::cout << " * X INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;

  	init_stalls_cycle_integrality_gap = std::numeric_limits<double>::infinity();
  	curr_stalls_cycle_improvement = std::numeric_limits<double>::infinity();
  	this->previous_integrality_gap_ = this->curr_integrality_gap_ = std::numeric_limits<double>::infinity();
    stalls_cycle_count = 0;
  	stage_2_iter = 0;
  	(this->curr_alpha_) = initial_alpha_stage2;

  	// Stage 2: only considers x variables in the new objective function
  	if(!(this->found_int_x_))
  	{
  		// Change objective to minimization
  		this->objective_.setSense(IloObjective::Minimize);
  		do
  		{
        this->previous_alpha_ = this->curr_alpha_;
  			(this->curr_alpha_)*= alpha_decrease_rate;
  			// Update objective function according to distance from current rounded (integer) y values
  			SetNewObjStage2();
  			if(double_less(fabs(this->curr_alpha_ - this->previous_alpha_),K_APLHA_DECREMENT_PRECISION)) stalls_cycle_count++;
        //if(stalls_cycle_count > 0) std::cout << stalls_cycle_count << std::endl;
  			stage_2_iter++;

  			//std::cout << "alpha: " << this->curr_alpha_ << std::endl;
  			//std::cout << "iter: " << stage_2_iter << std::endl;
  			// resolve
  			cplex_.extract(model_);
  			//cplex_.exportModel("FP.lp");
  			//getchar(); getchar();
        if (!(cplex_.solve()))
        {
          if((cplex_.getCplexStatus() == IloCplex::Infeasible)||(cplex_.getCplexStatus() == IloCplex::InfOrUnbd)) (this->solution_).is_infeasible_ = true;

          (this->solution_).time_stage2_ = timer->CurrentElapsedTime(ti);

          (this->solution_).num_iterations_stage2_ = stage_2_iter;
          (this->solution_).num_perturbations_stage2_ = num_perturbations_stage2;
          (this->solution_).num_restarts_stage2_ = num_restarts_stage2;
          return;
        }
  			// update y values
  			this->RetrieveAndRoundArcValuesBasic();

  			if(!(this->found_int_x_))
  			{
  				//std::cout << "curr_integrality_gap: " << this->curr_integrality_gap_ << std::endl;
  				//std::cout << "curr_improv: " << fabs(this->previous_integrality_gap_ - this->curr_integrality_gap_) << std::endl;
  				curr_stalls_cycle_improvement = fabs(init_stalls_cycle_integrality_gap - this->curr_integrality_gap_);
  				//std::cout << "stalls_cycle_improv: " << curr_stalls_cycle_improvement << std::endl;
  				if(double_greater(curr_stalls_cycle_improvement,K_PUMP_IMPROVEMENT_TOLERANCE))
  				{
  					stalls_cycle_count = 0;
  					init_stalls_cycle_integrality_gap = this->curr_integrality_gap_;
  				}

                /*if (stalls_cycle_count >= max_stalls_stage2)
  				{
                    num_restarts_stage2++;
          			// means that a "cycle" has been heuristically detected: PERFORM RESTART
                    std::cout << "ciclou" << std::endl;

                    for(int i  = 0; i < num_arcs; ++i)
                    {
                        if((this->curr_int_x_)[i] == (this->previous_int_x_)[i])
                        {
                          //std::cout << (this->curr_relax_y_)[i] << " " << ((this->curr_int_y_)[i]) << std::endl;
                          //std::cout << fabs((this->curr_relax_y_)[i] - 1.0*((this->curr_int_y_)[i])) << std::endl;
                          curr_prob_of_flipping = ceil((fabs((this->curr_relax_x_)[i] - 1.0*((this->curr_int_x_)[i])) + 0.03)*100.0);
                          //std::cout << "flip prob: " << curr_prob_of_flipping << std::endl;
                          //int curr_rand = rand()%101;
                          //std::cout << "curr rand: " << curr_rand << std::endl;
                          if(rand()%101 <= curr_prob_of_flipping)
                          {
                            //std::cout << "flippou!" << std::endl;
                            ((this->curr_int_x_)[i]).flip();
                          }
                          //getchar();
                        }
                    }

                    stalls_cycle_count = 0;
  					init_stalls_cycle_integrality_gap = this->curr_integrality_gap_;

                    // if detected cycle of size one, PERFORM PERTURBATION
  				}else */if((stage_2_iter > 1) && (double_less(fabs(this->curr_alpha_ - this->previous_alpha_),K_APLHA_DECREMENT_PRECISION)) && (this->curr_int_x_ == this->previous_int_x_))
    			{
                      num_perturbations_stage2++;
                      curr_num_flips = num_flips_basis/2 + 1 + rand()%num_flips_basis;
                      int cont = 0;
                      //std::cout << "num_flips: " << curr_num_flips << std::endl;
                      //getchar();getchar();
            				 //std::cout << " *** ciclou!" << std::endl;
                     (this->curr_integrality_gaps_x_).sort(compare_func3);
                     for(auto it = (this->curr_integrality_gaps_x_).begin(); it != (this->curr_integrality_gaps_x_).end(); ++it)
                     {
                       //std::cout << (this->curr_int_x_[it->first]) << " - ";
                       if(cont < curr_num_flips) ((this->curr_int_x_)[it->first]).flip();
                       else break;
                       //std::cout << (this->curr_int_x_[it->first]) << std::endl;
                       //getchar();getchar();
                       cont++;
                     }
                     //std::cout << std::endl;
            				 //getchar();getchar();
    		    }
  			}
  		}while((!(this->found_int_x_))&&(stage_2_iter < max_iter_stage2));
  	}

    (this->solution_).time_stage2_ = timer->CurrentElapsedTime(ti);

    (this->solution_).num_iterations_stage1_ = stage_1_iter + stage_0_iter;
    (this->solution_).num_iterations_stage2_ = stage_2_iter;
    (this->solution_).num_perturbations_stage1_ = num_perturbations_stage1;
    (this->solution_).num_perturbations_stage2_ = num_perturbations_stage2;
    (this->solution_).num_restarts_stage1_ = num_restarts_stage1;
    (this->solution_).num_restarts_stage2_ = num_restarts_stage2;

    if(this->found_int_x_) (this->solution_).is_feasible_ = (this->solution_).found_y_integer_ = (this->solution_).found_x_integer_ = true;
    if(this->found_int_y_) (this->solution_).found_y_integer_ = true;

    if(this->found_int_x_) BuildHeuristicSolution();
  	//std::cout << "iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
  	//if(this->found_int_y_) std::cout << " * Y INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
  	//if(this->found_int_x_) std::cout << " * X INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;

    //getchar(); getchar();
  }

  void FeasibilityPump::BuildModel(double * R0, double* Rn)
  {
  	Graph * graph = (this->curr_instance_)->graph();
  	int num_vertices = graph->num_vertices();
  	int num_mandatory = graph->num_mandatory();

    std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();
    (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;

    if(K_FEASIBILITY_PUMP_ADD_CUTS)
    {
      //(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
      //(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
      //(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
      (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;
    }

  	CapacitatedSingleCommodityForm(this->env_, this->model_, this->slack_, y, x, f, *(this->curr_instance_), R0, Rn, true);

  	(this->cplex_).extract(this->model_);
  	this->objective_ = (this->cplex_).getObjective();
  	this->original_obj_expr_ = this->objective_.getExpr();
  	this->original_obj_norm_ = 0.0;

  	for(int i = num_mandatory + 1; i < num_vertices-1; i++) (this->original_obj_norm_)+= pow(1.0*((graph->profits())[i]),2.0);

  	this->original_obj_norm_ = sqrt(this->original_obj_norm_);
  	//cplex_.extract(model_);
  	//cplex_.exportModel("FP.lp");
  }
