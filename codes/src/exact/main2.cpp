#include <iostream>
#include <string>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include "src/instance.h"
#include "src/graph.h"
#include "src/general.h"
#include "src/formulations.h"
#include "src/timer.h"
#include "src/graph_algorithms.h"

int main()
{
	Instance inst("../instances/R-STOP-DP/","test.txt",1,0.25,3,false);
	

	std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();
	(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
	(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = false;

	std::cout << inst << std::endl;
	const Graph * graph = inst.graph();

	//std::cout << *graph << std::endl;

	double * R0 = Dijkstra(graph,false,false);
	double * Rn = Dijkstra(graph,true,false);

	//inst.WriteToFile("../","teste.txt");
	bool force_use_all_vehicles = false;
	bool export_model = false;
	bool solve_relaxed = true;
	bool use_valid_inequalities = true;
	auto root_cuts = new std::list<UserCutGeneral*>();
	Solution<double> solution(graph->num_vertices());

	// CompactBaseline(inst,R0,Rn,-1,true,false,false,nullptr,nullptr,force_use_all_vehicles,export_model, nullptr,solution);
	// std::cout <<  " LP: " << solution.lp_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();
	// CompactBaseline(inst,R0,Rn,-1,true,use_valid_inequalities,true,nullptr,nullptr,force_use_all_vehicles,export_model, root_cuts,solution);
	// std::cout <<  " LP: " << solution.lp_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();
	// CompactBaseline(inst,R0,Rn,-1,false,use_valid_inequalities,false,root_cuts,nullptr,force_use_all_vehicles,export_model, nullptr,solution);
	
	// if(!solution.is_feasible_)
	// 	std::cout << "Infeasible" << std::endl;
	// solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
	//solution.reset();
	//DeleteCuts(root_cuts);

	CompactSingleCommodity(inst,R0,Rn,-1,true,false,false,nullptr,nullptr,force_use_all_vehicles,export_model, nullptr,solution);

	std::cout <<  " LP: " << solution.lp_ << std::endl;
	std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	solution.reset();
	CompactSingleCommodity(inst,R0,Rn,-1,true,use_valid_inequalities,true,nullptr,nullptr,force_use_all_vehicles,export_model, root_cuts,solution);

	std::cout <<  " LP: " << solution.lp_ << std::endl;
	std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	solution.reset();
	// auto conflicts = inst.conflicts_list();
	// std::cout << "conflicts: " << std::endl;
	// for (auto &conflict: conflicts)
	// {
	// 	for( auto& vertex: conflict)
	// 		std::cout << vertex << " ";
		
	// 	std::cout << std::endl;
	// }

	// CompactSingleCommodity(inst,R0,Rn,-1,false,use_valid_inequalities,false,root_cuts,nullptr,force_use_all_vehicles,export_model,nullptr,solution);

	// if(!solution.is_feasible_)
	// 	std::cout << "Infeasible" << std::endl;
	// solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	solution.reset();
	DeleteCuts(root_cuts);

	BendersCompactBaseline(inst,R0,Rn,-1,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model,solution);
	if (solve_relaxed)
	{
		std::cout <<  " LP: " << solution.lp_ << std::endl;
	}
	else
	{
		solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
		: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	}
	std::cout << "# opt cuts: " << inst.num_opt_cuts_ << std::endl;
	std::cout << "# feas cuts: " << inst.num_feas_cuts_ << std::endl;
	std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
	solution.reset();
	// IloEnv env;
	// IloNumArray y_values(env,inst.graph()->num_vertices());
	// // for(int i = 0; i < inst.graph()->num_vertices(); ++i)
	// // {
	// //   y_values[i] = 0;
	// // }
	// y_values[0] = 1;
	// y_values[2] = 1;
	// y_values[1] = y_values[3] = 1;

	// IloNumArray x_values(env,inst.graph()->num_arcs());
	
	// for(int i  = 0; i < inst.graph()->num_vertices(); ++i)
	// {
	//   for (const auto & j: inst.graph()->AdjVerticesOut(i))
	//   {
	//     x_values[inst.graph()->pos(i,j)] = 0;
	//   }
	// }

	// x_values[inst.graph()->pos(0,1)] = 1;
	// x_values[inst.graph()->pos(1,2)] = 1;
	// x_values[inst.graph()->pos(2,3)] = 1;
	// x_values[inst.graph()->pos(3,0)] = 1;
	
	// auto PrimalSubproblemCompactBaseline(inst,x_values,y_values,nullptr,nullptr,-1,export_model,solution);
	// solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// DualSubproblemCompactBaseline(inst,x_values,y_values,nullptr,nullptr,-1,export_model,solution);
	// solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;

	// env.end();

	delete [] R0; 
	R0 = nullptr;
	delete [] Rn;
	Rn = nullptr;

	DeleteTimer();
	DeleteCallbackSelection();

	return 0;
}
