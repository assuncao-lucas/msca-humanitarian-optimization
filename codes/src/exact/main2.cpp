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

int main()
{
	Instance inst("../instances/R-STOP-DP/","R25_1.txt",1,0.5,1,false);
	//std::cout << inst << std::endl;
	//inst.WriteToFile("../","teste.txt");
	bool force_use_all_vehicles = true;
	bool export_model = false;
	// auto solution = CompactBaseline(inst,nullptr,nullptr,-1,false,false,false,nullptr,nullptr,force_use_all_vehicles,export_model);
	// solution->is_optimal_? std::cout <<  " optimal: " << solution->lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution->lb_ << ", " << solution->ub_ << "]" << std::endl;
	// delete solution;
	auto solution = BendersCompactBaseline(inst,nullptr,nullptr,-1,false,true,false,nullptr,nullptr,force_use_all_vehicles,export_model);
	solution->is_optimal_? std::cout <<  " optimal: " << solution->lb_ << std::endl
	: std::cout <<  " non optimal: [" << solution->lb_ << ", " << solution->ub_ << "]" << std::endl;

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
	
	// auto solution = PrimalSubproblemCompactBaseline(inst,x_values,y_values,nullptr,nullptr,-1,export_model);
	// solution->is_optimal_? std::cout <<  " optimal: " << solution->lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution->lb_ << ", " << solution->ub_ << "]" << std::endl;
	// delete solution;
	// solution = DualSubproblemCompactBaseline(inst,x_values,y_values,nullptr,nullptr,-1,export_model);
	// solution->is_optimal_? std::cout <<  " optimal: " << solution->lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution->lb_ << ", " << solution->ub_ << "]" << std::endl;

	// env.end();
	delete solution;
	solution = nullptr;
	DeleteTimer();
	return 0;
}
