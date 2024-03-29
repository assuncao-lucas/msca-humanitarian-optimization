#include <cmath>
#include <algorithm>
#include "src/feasibility_pump.h"
#include "src/graph_algorithms.h"
#include "src/formulations.h"
#include "src/general.h"

bool compare_func3(const std::pair<int, double> &v1, const std::pair<int, double> &v2)
{
	return double_greater(v1.second, v2.second);
}

FeasibilityPump::FeasibilityPump()
{
	original_obj_norm_ = 0.0;
	curr_alpha_ = 0.0;
	previous_alpha_ = std::numeric_limits<double>::infinity();
	found_int_y_ = false;
	found_int_x_ = false;
	previous_integrality_gap_ = 0.0;
	curr_integrality_gap_ = 0.0;
}

FeasibilityPump::~FeasibilityPump()
{
	Reset();
}

void FeasibilityPump::Reset()
{
	cplex_.end();
	model_.end();
	env_.end();
}

void FeasibilityPump::Init(Instance &instance)
{
	curr_instance_ = &instance;
	const Graph *graph = instance.graph();
	int num_vertices = graph->num_vertices();
	int num_arcs = graph->num_arcs();
	int num_routes = instance.num_vehicles();
	int budget = curr_instance_->uncertainty_budget();
	solution_.Reset(num_vertices, num_arcs, instance.num_vehicles());

	curr_alpha_ = 0.0;

	env_ = IloEnv();
	model_ = IloModel(env_);
	cplex_ = IloCplex(env_);
	cplex_.setOut(env_.getNullStream());

	curr_relax_y_ = IloNumArray(env_);
	curr_relax_x_ = IloNumArray(env_);

	curr_int_y_ = boost::dynamic_bitset<>(num_vertices, 0);
	curr_int_x_ = boost::dynamic_bitset<>(num_arcs, 0);

	double *R0 = Dijkstra(graph, false, false);
	double *Rn = Dijkstra(graph, true, false);

	AllocateMasterVariablesSingleCommodity(env_, master_vars_, *curr_instance_, false, true);
	// budget + 1 to consider level 0 of budget! 0,..., budget
	f_ = IloNumVarArray(env_, num_arcs * (budget + 1), 0, IloInfinity, ILOFLOAT);
	slack_ = IloNumVar(env_, 0, num_routes, ILOFLOAT);

	BuildModel(R0, Rn);

	delete[] R0;
	delete[] Rn;
	R0 = nullptr;
	Rn = nullptr;
}

void FeasibilityPump::RetrieveAndRoundVertexValues()
{
	cplex_.getValues(curr_relax_y_, master_vars_.y);
	const Graph *graph = (curr_instance_)->graph();
	int num_vertices = graph->num_vertices(), num_mandatory = (curr_instance_)->num_mandatory();
	previous_int_y_ = curr_int_y_;
	(curr_int_y_).reset();
	(curr_integrality_gaps_y_).clear();
	found_int_y_ = true;

	previous_integrality_gap_ = curr_integrality_gap_;
	curr_integrality_gap_ = 0.0;

	for (int i = num_mandatory + 1; i < num_vertices - 1; i++)
	{
		if (!double_equals(curr_relax_y_[i], 0.0))
		{
			double curr_integrality_gap_y = 0.0;
			if (round((curr_relax_y_)[i]))
				(curr_int_y_)[i] = 1;
			if (!double_equals((curr_relax_y_)[i], 1.0 * ((curr_int_y_)[i])))
			{
				curr_integrality_gap_y = fabs((curr_relax_y_)[i] - 1.0 * ((curr_int_y_)[i]));
				found_int_y_ = false;
				(curr_integrality_gap_) += curr_integrality_gap_y;
			}

			curr_integrality_gaps_y_.push_back(std::pair<int, double>(i, curr_integrality_gap_y));
		}
	}
}

bool FeasibilityPump::BuildPathsBFSIter(std::vector<std::list<std::pair<int, double>>> &residual_graph, std::list<std::list<int>> &paths,
										std::list<int> curr_path, std::vector<char> &colors, size_t vertex)
{
	if (colors[vertex] != 'w')
	{
		paths.push_back(curr_path);
		return true;
	}

	curr_path.push_back(vertex);
	if (vertex == residual_graph.size() - 1)
	{
		paths.push_back(curr_path);
		return true;
	}
	else
	{
		colors[vertex] = 'g';
		bool found_path = false;
		for (auto it = residual_graph[vertex].begin(); it != residual_graph[vertex].end(); ++it)
		{
			if (BuildPathsBFSIter(residual_graph, paths, curr_path, colors, it->first))
			{
				colors[vertex] = 'b';
				found_path = true;
			}
		}
		if (found_path)
			return true;
	}
	colors[vertex] = 'w';
	return false;
}

void FeasibilityPump::BuildPathsBFS(std::vector<std::list<std::pair<int, double>>> &residual_graph, std::list<std::list<int>> &paths)
{
	std::list<int> curr_path;
	std::vector<char> colors(residual_graph.size(), 'w');
	BuildPathsBFSIter(residual_graph, paths, curr_path, colors, 0);
}

void FeasibilityPump::RetrieveAndRoundArcValuesBFS()
{
	cplex_.getValues(curr_relax_x_, master_vars_.x);
	const Graph *graph = (curr_instance_)->graph();
	int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
	int v1 = 0;
	std::list<std::list<int>> paths;
	(curr_int_x_).reset();

	std::vector<std::list<std::pair<int, double>>> residual_graph(num_vertices, std::list<std::pair<int, double>>());

	std::list<int> q;
	q.push_front(0);
	std::vector<bool> visited_nodes(num_vertices, false);
	visited_nodes[0] = true;
	int curr_pos = 0;

	found_int_x_ = true;
	previous_integrality_gap_ = curr_integrality_gap_;
	curr_integrality_gap_ = 0.0;

	// we don't have to check every arc: due to the flow conservation constraints, only perform width first search
	do
	{
		v1 = q.front();
		q.pop_front();

		for (int v2 : graph->AdjVerticesOut(v1))
		{
			curr_pos = graph->pos(v1, v2);

			if (!double_equals(curr_relax_x_[curr_pos], 0.0))
			{
				(residual_graph[v1]).push_back(std::pair<int, double>(v2, curr_relax_x_[curr_pos]));
				if (!(visited_nodes[v2]))
				{
					q.push_front(v2);
					visited_nodes[v2] = true;
				}
			}
		}
	} while (!q.empty());

	// sort each list in residual_graph in non-increasing order of relaxation value
	for (int i = 0; i < num_vertices; i++)
	{
		(residual_graph[i]).sort(compare_func3);
	}

	// Call BFS to build the paths_list
	BuildPathsBFS(residual_graph, paths);

	int count = 0;
	for (auto it = paths.begin(); it != paths.end(); ++it)
	{
		++count;
		std::cout << "path " << count << std::endl;
		for (auto it2 = it->begin(); it2 != it->end(); ++it2)
		{
			int v2 = *it2;
			if (it2 != it->begin())
			{
				curr_pos = graph->pos(v1, v2);
				(curr_int_x_)[curr_pos] = 1;
			}
			v1 = v2;
			std::cout << *it2 << " ";
		}
		std::cout << std::endl;
		getchar();
		getchar();
	}

	curr_integrality_gap_ = 0.0;
	for (curr_pos = 0; curr_pos < num_arcs; curr_pos++)
	{
		if (!double_equals((curr_relax_x_)[curr_pos], (curr_int_x_)[curr_pos]))
		{
			found_int_x_ = false;
			(curr_integrality_gap_) += fabs((curr_relax_x_)[curr_pos] - 1.0 * ((curr_int_x_)[curr_pos]));
		}
	}
}

void FeasibilityPump::RetrieveAndRoundArcValuesBasic()
{
	cplex_.getValues(curr_relax_x_, master_vars_.x);
	const Graph *graph = (curr_instance_)->graph();
	int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
	int v1 = 0;

	std::list<int> q;
	q.push_front(0);
	std::vector<bool> visited_nodes(num_vertices, false);
	visited_nodes[0] = true;
	int curr_pos = 0;

	found_int_x_ = true;
	previous_integrality_gap_ = curr_integrality_gap_;
	curr_integrality_gap_ = 0.0;

	previous_int_x_ = curr_int_x_;
	(curr_int_x_).reset();
	(curr_integrality_gaps_x_).clear();

	// we don't have to check every arc: due to the flow conservation constraints, only perform width first search
	// do
	// {
	// 	v1 = q.front();
	// 	q.pop_front();

	// 	for (int v2 : graph->AdjVerticesOut(v1))
	// 	{
	// 		curr_pos = graph->pos(v1, v2);

	// 		if (!double_equals(curr_relax_x_[curr_pos], 0.0))
	// 		{
	// 			double curr_integrality_gap_x = 0.0;
	// 			if (!(visited_nodes[v2]))
	// 			{
	// 				q.push_front(v2);
	// 				visited_nodes[v2] = true;
	// 			}

	// 			// if(double_greater((curr_relax_x_)[curr_pos],0.0)) (curr_int_x_)[curr_pos] = 1;
	// 			if (round((curr_relax_x_)[curr_pos]))
	// 				(curr_int_x_)[curr_pos] = 1;
	// 			if (!double_equals((curr_relax_x_)[curr_pos], (curr_int_x_)[curr_pos]))
	// 			{
	// 				found_int_x_ = false;
	// 				curr_integrality_gap_x = fabs((curr_relax_x_)[curr_pos] - 1.0 * ((curr_int_x_)[curr_pos]));
	// 				(curr_integrality_gap_) += curr_integrality_gap_x;
	// 			}

	// 			curr_integrality_gaps_x_.push_back(std::pair<int, double>(curr_pos, curr_integrality_gap_x));
	// 		}
	// 	}
	// } while (!q.empty());

	for (int v1 = 0; v1 < num_vertices; ++v1)
	{
		for (auto v2 : graph->AdjVerticesOut(v1))
		{
			curr_pos = curr_pos = graph->pos(v1, v2);
			if (!double_equals(curr_relax_x_[curr_pos], 0.0))
			{
				double curr_integrality_gap_x = 0.0;
				if (round((curr_relax_x_)[curr_pos]))
					(curr_int_x_)[curr_pos] = 1;
				if (!double_equals((curr_relax_x_)[curr_pos], (curr_int_x_)[curr_pos]))
				{
					found_int_x_ = false;
					curr_integrality_gap_x = fabs((curr_relax_x_)[curr_pos] - 1.0 * ((curr_int_x_)[curr_pos]));
					(curr_integrality_gap_) += curr_integrality_gap_x;
				}

				curr_integrality_gaps_x_.push_back(std::pair<int, double>(curr_pos, curr_integrality_gap_x));
			}
		}
	}

	if (found_int_x_)
	{
		for (int i = 0; i < num_vertices; ++i)
		{
			for (int j = 0; j < num_vertices; ++j)
				if ((i != j) && !(double_equals(curr_relax_x_[curr_instance_->graph()->pos(j, i)], 0.0)))
					std::cout << j << ", " << i << " :" << curr_relax_x_[curr_instance_->graph()->pos(j, i)] << std::endl;
		}
	}
}

void FeasibilityPump::SetNewObjStage1()
{
	IloExpr new_obj(env_);
	const Graph *graph = curr_instance_->graph();
	int num_vertices = graph->num_vertices(), num_mandatory = curr_instance_->num_mandatory();
	// int * profits = graph->profits();
	// curr_obj_const_ = 0;

	for (int i = num_mandatory + 1; i < num_vertices; ++i)
	{
		if ((curr_int_y_)[i] == 0)
		{
			new_obj += master_vars_.y[i];
		}
		else if ((curr_int_y_)[i] == 1)
		{
			new_obj += (1 - master_vars_.y[i]);
			//(curr_obj_const_)++;
		}
		else
		{
			throw std::string("Found a binary var that, once rounded, is neither 0 nor 1. This shouldn't happen!");
		}
	}

	(objective_).setExpr(operator*((1.0 - curr_alpha_) / sqrt(num_vertices - num_mandatory - 1), new_obj) - operator*((curr_alpha_) / original_obj_norm_, original_obj_expr_));
	new_obj.end();
}

void FeasibilityPump::SetNewObjStage2()
{
	IloExpr new_obj(env_);
	const Graph *graph = (curr_instance_)->graph();
	int num_vertices = graph->num_vertices(), num_mandatory = curr_instance_->num_mandatory(), num_arcs = graph->num_arcs();
	int v1 = 0, v2 = 0, curr_pos = 0;
	// int * profits = graph->profits();
	// curr_obj_const_ = 0;

	for (int i = 0; i < num_vertices; ++i)
	{
		v1 = i;
		for (int v2 : graph->AdjVerticesOut(v1))
		{
			curr_pos = graph->pos(v1, v2);
			if ((curr_int_x_)[curr_pos] == 0)
			{
				new_obj += master_vars_.x[curr_pos];
			}
			else if ((curr_int_x_)[curr_pos] == 1)
			{
				new_obj += (1 - master_vars_.x[curr_pos]);
				//(curr_obj_const_)++;
			}
			else
			{
				throw std::string("Found a binary var that, once rounded, is neither 0 nor 1. This shouldn't happen!");
			}
		}
	}

	(objective_).setExpr(operator*((1.0 - curr_alpha_) / sqrt(num_arcs), new_obj) - operator*((curr_alpha_) / original_obj_norm_, original_obj_expr_));
	new_obj.end();
}

void FeasibilityPump::BuildHeuristicSolution()
{
	const Graph *graph = curr_instance_->graph();
	int num_vertices = graph->num_vertices(), v1 = 0, curr_route_index = 0;
	int last_vertex_added = 0;
	bool has_added_vertex_to_route = false;
	VertexStatus *status = nullptr;
	Route *curr_route = nullptr;
	GArc *curr_arc = nullptr;
	size_t cont = 0;

	std::list<int> q;
	q.push_back(0);

	do
	{
		v1 = q.front();
		q.pop_front();

		// end current route
		if ((v1 == 0) && has_added_vertex_to_route)
		{
			curr_route = &(((solution_).routes_vec_)[curr_route_index]);
			curr_arc = (*graph)[last_vertex_added][v1];
			if (curr_arc == nullptr)
				throw "Inappropriate addition of vertex to route";
			//   else
			//   {
			//     (curr_route->time_) += curr_arc->dist();
			//   }

			last_vertex_added = 0;

			// compute route's profit sum and max duration.
			auto [route_sum_profits, route_max_duration] = curr_instance_->ComputeRouteCostsRec(*curr_route, true);
			curr_route->time_ = route_max_duration;
			solution_.profits_sum_ += route_sum_profits;
			curr_route->sum_profits_ = route_sum_profits;

			std::cout << *curr_route << std::endl;

			++curr_route_index;
			has_added_vertex_to_route = false;
			continue; // in this case, should avoid the last step of the loop that adds to stack the neighbors of 0.
		}
		else if (v1 != 0)
		{
			// add vertex to current route
			curr_route = &(((solution_).routes_vec_)[curr_route_index]);
			status = &(((solution_).vertex_status_vec_)[v1]);
			has_added_vertex_to_route = true;

			// remove from list of unvisited_vertices
			((solution_).unvisited_vertices_).erase(status->pos_);

			// adds vertex to route
			status->selected_ = true;
			status->route_ = curr_route_index;
			status->pos_ = (curr_route->vertices_).insert((curr_route->vertices_).end(), v1);

			curr_arc = (*graph)[last_vertex_added][v1];
			if (curr_arc == nullptr)
				throw "Inappropriate addition of vertex to route";
			// else
			// {
			//   (curr_route->time_) += curr_arc->dist();
			// }
			last_vertex_added = v1;
		}

		for (int v2 : graph->AdjVerticesOut(v1))
		{
			if ((curr_int_x_)[graph->pos(v1, v2)])
			{
				++cont;
				q.push_front(v2);
			}
		}
	} while (!(q.empty()));

	// this is not true because of the decreasing profits...the order in which vertices appear change the profit collected.
	// if (((solution_).unvisited_vertices_).empty())
	// 	(solution_).is_optimal_ = true;

	(solution_).BuildBitset(*(curr_instance_));
	if ((cont != (curr_int_x_).count())) // || (!((solution_).CheckCorrectness(*(curr_instance_)))))
		throw "Error Building FP solution";
}

void FeasibilityPump::Run()
{
	Timestamp *ti = NewTimestamp();
	Timer *timer = GetTimer();
	timer->Clock(ti);
	std::cout << std::setprecision(2) << std::fixed;
	double init_stalls_cycle_integrality_gap = std::numeric_limits<double>::infinity();
	double curr_stalls_cycle_improvement = std::numeric_limits<double>::infinity();
	previous_integrality_gap_ = curr_integrality_gap_ = std::numeric_limits<double>::infinity();
	int stage_0_iter = 1, stage_1_iter = 0, stage_2_iter = 0, stalls_cycle_count = 0;
	int num_mandatory = curr_instance_->num_mandatory();
	int num_vertices = curr_instance_->graph()->num_vertices();
	int num_arcs = curr_instance_->graph()->num_arcs();
	int num_flips_basis = 0, curr_num_flips = 0;
	if (K_FLIP_BASIS)
		num_flips_basis = K_FLIP_BASIS;
	else
		num_flips_basis = ceil(perturbation_flip_percentage * num_vertices);
	int curr_prob_of_flipping = 0;
	int num_perturbations_stage1 = 0, num_perturbations_stage2 = 0, num_restarts_stage1 = 0, num_restarts_stage2 = 0;
	(curr_alpha_) = initial_alpha_stage1;
	// Stage 0: solve LP relaxation with original objective

	if (K_FEASIBILITY_PUMP_ADD_CUTS)
	{
		double *R0 = Dijkstra(((curr_instance_)->graph()), false, false);
		double *Rn = Dijkstra(((curr_instance_)->graph()), true, false);
		Solution<double> sol;

		auto root_cuts = new std::list<UserCutGeneral *>();

		bool use_valid_inequalities = false;
		bool find_root_cuts = true;

		CompactSingleCommodity(*curr_instance_, R0, Rn, -1, true, use_valid_inequalities, find_root_cuts, nullptr, nullptr, false, false, root_cuts, sol);

		DeleteCuts(root_cuts);
		if (root_cuts != nullptr)
		{
			delete root_cuts;
			root_cuts = nullptr;
		}

		delete[] R0;
		R0 = nullptr;
		delete[] Rn;
		Rn = nullptr;

		if (!(sol.is_feasible_))
		{
			(solution_).is_infeasible_ = true;
			(solution_).time_stage1_ = timer->CurrentElapsedTime(ti);
			return;
		}
	}
	else
	{
		cplex_.extract(model_);
		// cplex_.exportModel("FP.lp");
		// getchar(); getchar();

		if (!(cplex_.solve()))
		{
			if ((cplex_.getCplexStatus() == IloCplex::Infeasible) || (cplex_.getCplexStatus() == IloCplex::InfOrUnbd))
				(solution_).is_infeasible_ = true;
			(solution_).time_stage1_ = timer->CurrentElapsedTime(ti);
			return;
		}
	}

	RetrieveAndRoundVertexValues();

	init_stalls_cycle_integrality_gap = curr_integrality_gap_;

	// skips Stage 1 if all vertices are integer or K_SOLVE_STAGE_1 is false
	if (!(found_int_y_) && K_SOLVE_STAGE_1)
	{
		// Change objective to minimization
		objective_.setSense(IloObjective::Minimize);

		// Stage 1: only considers y variables in the new objective function
		do
		{
			previous_alpha_ = curr_alpha_;
			(curr_alpha_) *= alpha_decrease_rate;
			// std::cout << "alpha: " << curr_alpha_ << std::endl;
			//  Update objective function according to distance from current rounded (integer) y values
			SetNewObjStage1();

			if (double_less(fabs(curr_alpha_ - previous_alpha_), K_APLHA_DECREMENT_PRECISION))
				stalls_cycle_count++;
			// if(stalls_cycle_count > 0) std::cout << stalls_cycle_count << std::endl;
			stage_1_iter++;

			// resolve
			cplex_.extract(model_);
			// cplex_.exportModel("FP.lp");
			// getchar(); getchar();
			if (!(cplex_.solve()))
			{
				if ((cplex_.getCplexStatus() == IloCplex::Infeasible) || (cplex_.getCplexStatus() == IloCplex::InfOrUnbd))
					(solution_).is_infeasible_ = true;

				(solution_).time_stage1_ = timer->CurrentElapsedTime(ti);

				(solution_).num_iterations_stage1_ = stage_1_iter + stage_0_iter;
				(solution_).num_perturbations_stage1_ = num_perturbations_stage1;
				(solution_).num_restarts_stage1_ = num_restarts_stage1;
				return;
			}
			// update y values
			RetrieveAndRoundVertexValues();

			if (!(found_int_y_))
			{
				// std::cout << "curr_improv: " << fabs(previous_integrality_gap_ - curr_integrality_gap_) << std::endl;
				curr_stalls_cycle_improvement = fabs(init_stalls_cycle_integrality_gap - curr_integrality_gap_);
				// std::cout << "stalls_cycle_improv: " << curr_stalls_cycle_improvement << std::endl;
				if (double_greater(curr_stalls_cycle_improvement, K_PUMP_IMPROVEMENT_TOLERANCE))
				{
					stalls_cycle_count = 0;
					init_stalls_cycle_integrality_gap = curr_integrality_gap_;
				}

				if (stalls_cycle_count >= max_stalls_stage1)
				{
					num_restarts_stage1++;
					// means that a "cycle" has been heuristically detected: PERFORM RESTART
					// std::cout << "ciclou" << std::endl;
					for (int i = num_mandatory + 1; i < num_vertices - 1; i++)
					{
						if ((curr_int_y_)[i] == (previous_int_y_)[i])
						{
							// std::cout << (curr_relax_y_)[i] << " " << ((curr_int_y_)[i]) << std::endl;
							// std::cout << fabs((curr_relax_y_)[i] - 1.0*((curr_int_y_)[i])) << std::endl;
							curr_prob_of_flipping = ceil((fabs((curr_relax_y_)[i] - 1.0 * ((curr_int_y_)[i])) + 0.03) * 100.0);
							// std::cout << "flip prob: " << curr_prob_of_flipping << std::endl;
							// int curr_rand = rand()%101;
							// std::cout << "curr rand: " << curr_rand << std::endl;
							if (rand() % 101 <= curr_prob_of_flipping)
							{
								// std::cout << "flippou!" << std::endl;
								((curr_int_y_)[i]).flip();
							}
							// getchar();
						}
					}
					stalls_cycle_count = 0;
					init_stalls_cycle_integrality_gap = curr_integrality_gap_;
					// if detected cycle of size one, PERFORM PERTURBATION
				}
				else if ((stage_1_iter > 1) && (double_less(fabs(curr_alpha_ - previous_alpha_), K_APLHA_DECREMENT_PRECISION)) && (curr_int_y_ == previous_int_y_))
				{
					curr_num_flips = num_flips_basis / 2 + 1 + rand() % num_flips_basis;
					int cont = 0;
					num_perturbations_stage1++;
					// std::cout << "num_flips: " << curr_num_flips << std::endl;
					// getchar();getchar();
					// std::cout << " *** ciclou!" << std::endl;
					(curr_integrality_gaps_y_).sort(compare_func3);
					for (auto it = (curr_integrality_gaps_y_).begin(); it != (curr_integrality_gaps_y_).end(); ++it)
					{
						// std::cout << it->second << std::endl;
						// std::cout << (curr_int_y_[it->first]) << " - ";
						if (cont < curr_num_flips)
							((curr_int_y_)[it->first]).flip();
						else
							break;
						// std::cout << (curr_int_y_[it->first]) << std::endl;
						// getchar();getchar();
						cont++;
					}
					// std::cout << std::endl;
					// getchar();getchar();
				}
			}

		} while ((!(found_int_y_)) && (stage_1_iter < max_iter_stage1));
	}

	RetrieveAndRoundArcValuesBasic();

	(solution_).time_stage1_ = timer->CurrentElapsedTime(ti);
	timer->Clock(ti);
	// std::cout << "iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
	// if(found_int_y_) std::cout << " * Y INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
	// if(found_int_x_) std::cout << " * X INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;

	init_stalls_cycle_integrality_gap = std::numeric_limits<double>::infinity();
	curr_stalls_cycle_improvement = std::numeric_limits<double>::infinity();
	previous_integrality_gap_ = curr_integrality_gap_ = std::numeric_limits<double>::infinity();
	stalls_cycle_count = 0;
	stage_2_iter = 0;
	(curr_alpha_) = initial_alpha_stage2;

	// Stage 2: only considers x variables in the new objective function
	if (!(found_int_x_))
	{
		// Change objective to minimization
		objective_.setSense(IloObjective::Minimize);
		do
		{
			previous_alpha_ = curr_alpha_;
			(curr_alpha_) *= alpha_decrease_rate;
			// Update objective function according to distance from current rounded (integer) y values
			SetNewObjStage2();
			if (double_less(fabs(curr_alpha_ - previous_alpha_), K_APLHA_DECREMENT_PRECISION))
				stalls_cycle_count++;
			// if(stalls_cycle_count > 0) std::cout << stalls_cycle_count << std::endl;
			stage_2_iter++;

			// std::cout << "alpha: " << curr_alpha_ << std::endl;
			// std::cout << "iter: " << stage_2_iter << std::endl;
			//  resolve
			cplex_.extract(model_);
			// cplex_.exportModel("FP.lp");
			// getchar(); getchar();
			if (!(cplex_.solve()))
			{
				if ((cplex_.getCplexStatus() == IloCplex::Infeasible) || (cplex_.getCplexStatus() == IloCplex::InfOrUnbd))
					(solution_).is_infeasible_ = true;

				(solution_).time_stage2_ = timer->CurrentElapsedTime(ti);

				(solution_).num_iterations_stage2_ = stage_2_iter;
				(solution_).num_perturbations_stage2_ = num_perturbations_stage2;
				(solution_).num_restarts_stage2_ = num_restarts_stage2;
				return;
			}
			// update y values
			RetrieveAndRoundArcValuesBasic();

			if (!(found_int_x_))
			{
				// std::cout << "curr_integrality_gap: " << curr_integrality_gap_ << std::endl;
				// std::cout << "curr_improv: " << fabs(previous_integrality_gap_ - curr_integrality_gap_) << std::endl;
				curr_stalls_cycle_improvement = fabs(init_stalls_cycle_integrality_gap - curr_integrality_gap_);
				// std::cout << "stalls_cycle_improv: " << curr_stalls_cycle_improvement << std::endl;
				if (double_greater(curr_stalls_cycle_improvement, K_PUMP_IMPROVEMENT_TOLERANCE))
				{
					stalls_cycle_count = 0;
					init_stalls_cycle_integrality_gap = curr_integrality_gap_;
				}

				/*if (stalls_cycle_count >= max_stalls_stage2)
				{
					num_restarts_stage2++;
					// means that a "cycle" has been heuristically detected: PERFORM RESTART
					std::cout << "ciclou" << std::endl;

					for(int i  = 0; i < num_arcs; ++i)
					{
						if((curr_int_x_)[i] == (previous_int_x_)[i])
						{
						  //std::cout << (curr_relax_y_)[i] << " " << ((curr_int_y_)[i]) << std::endl;
						  //std::cout << fabs((curr_relax_y_)[i] - 1.0*((curr_int_y_)[i])) << std::endl;
						  curr_prob_of_flipping = ceil((fabs((curr_relax_x_)[i] - 1.0*((curr_int_x_)[i])) + 0.03)*100.0);
						  //std::cout << "flip prob: " << curr_prob_of_flipping << std::endl;
						  //int curr_rand = rand()%101;
						  //std::cout << "curr rand: " << curr_rand << std::endl;
						  if(rand()%101 <= curr_prob_of_flipping)
						  {
							//std::cout << "flippou!" << std::endl;
							((curr_int_x_)[i]).flip();
						  }
						  //getchar();
						}
					}

					stalls_cycle_count = 0;
					init_stalls_cycle_integrality_gap = curr_integrality_gap_;

					// if detected cycle of size one, PERFORM PERTURBATION
				}else */
				if ((stage_2_iter > 1) && (double_less(fabs(curr_alpha_ - previous_alpha_), K_APLHA_DECREMENT_PRECISION)) && (curr_int_x_ == previous_int_x_))
				{
					num_perturbations_stage2++;
					curr_num_flips = num_flips_basis / 2 + 1 + rand() % num_flips_basis;
					int cont = 0;
					// std::cout << "num_flips: " << curr_num_flips << std::endl;
					// getchar();getchar();
					// std::cout << " *** ciclou!" << std::endl;
					(curr_integrality_gaps_x_).sort(compare_func3);
					for (auto it = (curr_integrality_gaps_x_).begin(); it != (curr_integrality_gaps_x_).end(); ++it)
					{
						// std::cout << (curr_int_x_[it->first]) << " - ";
						if (cont < curr_num_flips)
							((curr_int_x_)[it->first]).flip();
						else
							break;
						// std::cout << (curr_int_x_[it->first]) << std::endl;
						// getchar();getchar();
						cont++;
					}
					// std::cout << std::endl;
					// getchar();getchar();
				}
			}
			std::cout << stage_2_iter << std::endl;
		} while ((!(found_int_x_)) && (stage_2_iter < max_iter_stage2));
	}

	(solution_).time_stage2_ = timer->CurrentElapsedTime(ti);

	(solution_).num_iterations_stage1_ = stage_1_iter + (K_SOLVE_STAGE_1 ? stage_0_iter : 0);
	(solution_).num_iterations_stage2_ = stage_2_iter + (!K_SOLVE_STAGE_1 ? stage_0_iter : 0);
	(solution_).num_perturbations_stage1_ = num_perturbations_stage1;
	(solution_).num_perturbations_stage2_ = num_perturbations_stage2;
	(solution_).num_restarts_stage1_ = num_restarts_stage1;
	(solution_).num_restarts_stage2_ = num_restarts_stage2;

	if (found_int_x_)
		(solution_).is_feasible_ = (solution_).found_y_integer_ = (solution_).found_x_integer_ = true;
	if (found_int_y_)
		(solution_).found_y_integer_ = true;

	if (found_int_x_)
		BuildHeuristicSolution();
	// std::cout << "iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
	// if(found_int_y_) std::cout << " * Y INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
	// if(found_int_x_) std::cout << " * X INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;

	// getchar(); getchar();
}

void FeasibilityPump::BuildModel(double *R0, double *Rn)
{
	const Graph *graph = curr_instance_->graph();
	int num_vertices = graph->num_vertices();
	int num_mandatory = curr_instance_->num_mandatory();

	// std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();
	// //(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;

	// if(K_FEASIBILITY_PUMP_ADD_CUTS)
	// {
	//   //(*CALLBACKS_SELECTION)[K_TYPE_GCC_CUT] = true;
	//   //(*CALLBACKS_SELECTION)[K_TYPE_CONFLICT_CUT] = true;
	//   //(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_FLOW_BOUNDS_CUT] = true;
	//   (*CALLBACKS_SELECTION)[K_TYPE_COVER_BOUND_CUT] = true;
	//   (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
	//   (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;
	// }

	cplex_.extract(model_);

	PopulateByRowCompactSingleCommodity(cplex_, env_, model_, master_vars_, f_, *curr_instance_, R0, Rn, false, false);

	cplex_.extract(model_);
	objective_ = cplex_.getObjective();
	original_obj_expr_ = objective_.getExpr();
	original_obj_norm_ = 0.0;
	const auto *vertices_info = graph->vertices_info();
	double route_limit = curr_instance_->limit();

	for (int j = num_mandatory + 1; j < num_vertices; ++j)
	{
		const auto &vertex_info = vertices_info[j];
		double var_coef = vertex_info.profit_ - vertex_info.decay_ratio_ * route_limit;
		original_obj_norm_ += pow(1.0 * var_coef, 2.0);

		for (const int &i : graph->AdjVerticesIn(j))
			original_obj_norm_ += pow(1.0 * vertex_info.decay_ratio_, 2.0);
	}

	original_obj_norm_ = sqrt(original_obj_norm_);
	// cplex_.extract(model_);
	// cplex_.exportModel("FP.lp");
}
